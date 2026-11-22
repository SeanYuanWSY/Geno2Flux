#!/usr/bin/env python3
"""Genome-to-proteome prediction & standardization entry point.

Bridges a raw genome assembly FASTA (plus an external gene predictor or an
already-predicted coding GFF3) into a standardized longest-isoform proteome,
QC summary, and ready-to-run AlphaGEM target/coordinate contracts:

1. **Predict coding GFF3** (or reuse ``--gff`` with ``--tool gff``):
   - ``annevo``: runs ANNEVO fungal nucleotide prediction + HMM decoding.
   - ``helixer``: runs Helixer ``Helixer.py`` or 3-step ``fasta2h5``/``HybridModel``/``helixer_post_bin``.
   - ``tiberius``: runs Tiberius ``tiberius.py`` fungal gene prediction.
   - ``funannotate``: runs RepeatModeler/RepeatMasker soft-masking (or ``--masked-genome``) + ``funannotate predict`` with UniProt/RefSeq protein evidence (``--reference-proteins``).
   - ``custom``: runs a user-supplied command template ``--predictor-cmd`` (replacing ``{genome}`` and ``{gff}``).
   - ``gff``: skips gene calling and standardizes an existing coding GFF3 directly.
2. **Standardize & QC** (``scripts/standardize_annotation.py`` + ``scripts/annotation_qc.py``):
   validates every gene/mRNA/CDS structure, checks table-1 translation, and selects
   one longest representative protein isoform per gene.
3. **Build AlphaGEM contracts** (``scripts/build_coordinate_contract.py``):
   writes ``inputs/proteome.fasta``, ``inputs/genome.fasta``, ``inputs/id_map.tsv``,
   ``inputs/target-input-contract.json``, ``inputs/coordinate-map.tsv``,
   ``inputs/coordinate-segments.tsv``, and ``inputs/gapfill-coordinate-contract.json``.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys


SCRIPTS_DIR = Path(__file__).resolve().parent
ROOT = SCRIPTS_DIR.parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import annotation_qc  # noqa: E402
import build_coordinate_contract  # noqa: E402
import standardize_annotation  # noqa: E402


COMMON_GPU_ENV = {
    "OMP_NUM_THREADS": "4",
    "TF_NUM_INTRAOP_THREADS": "4",
    "TF_NUM_INTEROP_THREADS": "2",
    "TF_FORCE_GPU_ALLOW_GROWTH": "true",
    "PYTHONDONTWRITEBYTECODE": "1",
}


def _run_step(name: str, argv: list[str], *, cwd: Path | None = None, env_overrides: dict[str, str] | None = None) -> dict[str, object]:
    env = dict(os.environ)
    if env_overrides:
        env.update(env_overrides)
    completed = subprocess.run(argv, cwd=cwd, env=env, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"Prediction step {name!r} failed with exit code {completed.returncode}: {argv}")
    return {"step": name, "argv": argv, "cwd": str(cwd or Path.cwd()), "exit_code": 0}


def _predict_gff(args: argparse.Namespace, raw_dir: Path) -> tuple[Path, list[dict[str, object]]]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    out_gff = raw_dir / "prediction.gff3"
    ledger: list[dict[str, object]] = []
    genome = args.genome.resolve()

    if args.tool == "gff":
        if args.gff is None or not args.gff.is_file():
            raise ValueError("--tool gff requires an existing --gff file")
        shutil.copy2(args.gff.resolve(), out_gff)
        return out_gff, [{"step": "reuse_gff", "source_gff": str(args.gff.resolve())}]

    if args.tool == "custom":
        if not args.predictor_cmd:
            raise ValueError("--tool custom requires --predictor-cmd (use {genome} and {gff} placeholders, or accept <genome> <gff> appended)")
        tokens = shlex.split(args.predictor_cmd)
        has_placeholders = any("{genome}" in t or "{gff}" in t for t in tokens)
        if has_placeholders:
            cmd = [t.replace("{genome}", str(genome)).replace("{gff}", str(out_gff)) for t in tokens]
        else:
            cmd = tokens + [str(genome), str(out_gff)]
        ledger.append(_run_step("custom_predict", cmd))
        if not out_gff.is_file():
            raise RuntimeError(f"Custom predictor did not create expected GFF3: {out_gff}")
        return out_gff, ledger

    tool_root = args.tool_root.resolve() if args.tool_root else None

    if args.tool == "annevo":
        annevo_dir = args.annevo_dir.resolve() if args.annevo_dir else (tool_root / "tools/annevo/annevo-current" if tool_root else None)
        python_bin = args.tool_python or (str(tool_root / "envs/annevo/bin/python") if tool_root else sys.executable)
        model_pt = args.model_weights.resolve() if args.model_weights else (tool_root / "weights/annevo/ANNEVO_Fungi.pt" if tool_root else None)
        if annevo_dir is None or model_pt is None:
            raise ValueError("annevo requires --tool-root or explicit --annevo-dir and --model-weights")
        pred_h5 = raw_dir / "prediction.h5"
        ledger.append(
            _run_step(
                "annevo_prediction",
                [
                    python_bin,
                    "prediction.py",
                    "-g",
                    str(genome),
                    "-m",
                    str(model_pt),
                    "-p",
                    str(pred_h5),
                    "-l",
                    args.lineage,
                    "--batch_size",
                    "1",
                    "--num_workers",
                    "2",
                    "--comp",
                ],
                cwd=annevo_dir,
                env_overrides=COMMON_GPU_ENV,
            )
        )
        ledger.append(
            _run_step(
                "annevo_decoding",
                [
                    python_bin,
                    "decoding.py",
                    "-g",
                    str(genome),
                    "-p",
                    str(pred_h5),
                    "-o",
                    str(out_gff),
                    "-t",
                    str(args.threads),
                ],
                cwd=annevo_dir,
                env_overrides=COMMON_GPU_ENV,
            )
        )
        return out_gff, ledger

    if args.tool == "helixer":
        model_h5 = args.model_weights.resolve() if args.model_weights else (tool_root / "weights/helixer/fungi_v0.3_a_0100.h5" if tool_root else None)
        if model_h5 is None:
            raise ValueError("helixer requires --model-weights (e.g. fungi.h5) or --tool-root")
        helixer_cli = args.helixer_bin or shutil.which("Helixer.py")
        if helixer_cli and not tool_root:
            ledger.append(
                _run_step(
                    "helixer_cli",
                    [
                        helixer_cli,
                        "--fasta-path",
                        str(genome),
                        "--gff-output-path",
                        str(out_gff),
                        "--model-filepath",
                        str(model_h5),
                        "--subsequence-length",
                        "21384",
                        "--batch-size",
                        "1",
                    ],
                    env_overrides=COMMON_GPU_ENV,
                )
            )
            return out_gff, ledger
        if tool_root is None:
            raise ValueError("helixer requires Helixer.py in PATH (--helixer-bin) or --tool-root")
        hel_bin = tool_root / "envs/helixer/bin"
        post_bins = sorted((tool_root / "tools/helixerpost").rglob("helixer_post_bin"))
        if not post_bins:
            raise FileNotFoundError(f"helixer_post_bin not found under {tool_root / 'tools/helixerpost'}")
        env = {
            **COMMON_GPU_ENV,
            "TF_USE_LEGACY_KERAS": "1",
            "PATH": str(hel_bin) + os.pathsep + os.environ.get("PATH", ""),
            "LD_LIBRARY_PATH": str(tool_root / "tools/helixerpost-libs"),
        }
        genome_h5 = raw_dir / "genome.h5"
        preds_h5 = raw_dir / "predictions.h5"
        ledger.append(
            _run_step(
                "helixer_encoding",
                [
                    str(hel_bin / "fasta2h5.py"),
                    "--fasta-path",
                    str(genome),
                    "--species",
                    args.case_id,
                    "--h5-output-path",
                    str(genome_h5),
                    "--subsequence-length",
                    "21384",
                ],
                env_overrides=env,
            )
        )
        ledger.append(
            _run_step(
                "helixer_prediction",
                [
                    str(hel_bin / "HybridModel.py"),
                    "--load-model-path",
                    str(model_h5),
                    "--test-data",
                    str(genome_h5),
                    "--prediction-output-path",
                    str(preds_h5),
                    "--val-test-batch-size",
                    "8",
                    "--overlap",
                    "--overlap-offset",
                    "10692",
                    "--core-length",
                    "16038",
                ],
                env_overrides=env,
            )
        )
        ledger.append(
            _run_step(
                "helixer_post",
                [
                    str(post_bins[0]),
                    str(genome_h5),
                    str(preds_h5),
                    "100",
                    "0.1",
                    "0.8",
                    "60",
                    str(out_gff),
                ],
                env_overrides=env,
            )
        )
        return out_gff, ledger

    if args.tool == "tiberius":
        python_bin = args.tool_python or (str(tool_root / "envs/tiberius/bin/python") if tool_root else sys.executable)
        if args.tiberius_dir:
            tib_src = args.tiberius_dir.resolve()
        elif tool_root:
            matches = sorted((tool_root / "tools/tiberius/tiberius-source").glob("Tiberius-*"))
            if not matches:
                raise FileNotFoundError(f"Tiberius-* directory not found under {tool_root}")
            tib_src = matches[0]
        else:
            raise ValueError("tiberius requires --tiberius-dir or --tool-root")
        cfg_yaml = args.model_weights.resolve() if args.model_weights else (tib_src / "model_cfg/fungi.yaml")
        ledger.append(
            _run_step(
                "tiberius_prediction",
                [
                    python_bin,
                    str(tib_src / "tiberius.py"),
                    "--genome",
                    str(genome),
                    "--model_cfg",
                    str(cfg_yaml),
                    "--out",
                    str(out_gff),
                    "--codingseq",
                    str(raw_dir / "cds.fna"),
                    "--protseq",
                    str(raw_dir / "proteins.faa"),
                    "--batch_size",
                    "1",
                    "--min_genome_seqlen",
                    "0",
                    "--group_size_limit",
                    "5000000",
                ],
                env_overrides=COMMON_GPU_ENV,
            )
        )
        return out_gff, ledger

    if args.tool == "funannotate":
        fun_bin = args.funannotate_bin or shutil.which("funannotate") or "funannotate"
        masked = args.masked_genome.resolve() if args.masked_genome else genome
        pred_dir = raw_dir / "funannotate_out"
        cmd = [
            fun_bin,
            "predict",
            "-i",
            str(masked),
            "-o",
            str(pred_dir),
            "--species",
            args.species,
            "--cpus",
            str(args.threads),
        ]
        if args.reference_proteins:
            cmd += ["--protein_evidence", str(args.reference_proteins.resolve())]
        ledger.append(_run_step("funannotate_predict", cmd))
        gffs = sorted((pred_dir / "predict_results").glob("*.gff3"))
        if len(gffs) != 1:
            raise RuntimeError(f"Expected exactly one funannotate GFF3 under {pred_dir / 'predict_results'}, found {len(gffs)}")
        shutil.copy2(gffs[0], out_gff)
        return out_gff, ledger

    raise ValueError(f"Unsupported tool: {args.tool!r}")


def run_predict_proteome(args: argparse.Namespace) -> dict[str, object]:
    out_root = args.output.resolve()
    if out_root.exists() and any(out_root.iterdir()):
        raise ValueError(f"Output directory must not exist or must be empty: {out_root}")
    out_root.mkdir(parents=True, exist_ok=True)

    raw_dir = out_root / "01_raw_prediction"
    std_dir = out_root / "02_standardized"
    qc_dir = out_root / "03_annotation_qc"
    inputs_dir = out_root / "inputs"

    predicted_gff, ledger = _predict_gff(args, raw_dir)
    method_tag = args.method_label or (args.tool if args.tool != "gff" else "gff")
    std_prefix = f"standardized.{method_tag}"

    std_summary = standardize_annotation.standardize(
        args.genome.resolve(),
        predicted_gff,
        method_tag,
        std_dir,
        args.gffread,
        prefix=std_prefix,
    )

    std_gff = std_dir / f"{std_prefix}.standardized.gff3"
    all_faa = std_dir / f"{std_prefix}.proteins.all_isoforms.faa"
    rep_faa = std_dir / f"{std_prefix}.proteins.longest_isoform.faa"
    map_tsv = std_dir / f"{std_prefix}.gene_transcript_protein_map.tsv"
    excl_tsv = std_dir / f"{std_prefix}.excluded_models.tsv"

    qc_summary = annotation_qc.build_qc(
        method=method_tag,
        genome_path=args.genome.resolve(),
        gff_path=std_gff,
        all_proteins_path=all_faa,
        representative_proteins_path=rep_faa,
        mapping_path=map_tsv,
        excluded_path=excl_tsv,
        output_dir=qc_dir,
    )

    contract_receipt = build_coordinate_contract.build(
        genome=args.genome.resolve(),
        gff=std_gff,
        proteins=rep_faa,
        case_id=args.case_id,
        output=inputs_dir,
        prefix=args.prefix,
        source_label=method_tag,
        keep_original_ids=args.keep_original_ids,
    )

    manifest = {
        "status": "completed",
        "tool": args.tool,
        "method_label": method_tag,
        "case_id": args.case_id,
        "prediction_steps": ledger,
        "standardization": std_summary,
        "qc_summary": qc_summary,
        "inputs_dir": str(inputs_dir),
        "proteome_fasta": str(inputs_dir / "proteome.fasta"),
        "proteome_sha256": contract_receipt["files"]["proteome.fasta"]["sha256"],
        "target_proteins": contract_receipt["target_proteins"],
    }
    (out_root / "predict_proteome_receipt.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--tool",
        choices=["annevo", "helixer", "tiberius", "funannotate", "custom", "gff"],
        required=True,
        help="Gene prediction method (or 'gff' to standardize an existing coding GFF3)",
    )
    p.add_argument("--genome", type=Path, required=True, help="Input genome assembly FASTA")
    p.add_argument("--gff", type=Path, help="Existing coding GFF3 (required when --tool gff)")
    p.add_argument("--output", "--output-dir", dest="output", type=Path, required=True, help="New output directory")
    p.add_argument("--case-id", default="sample", help="Sample / strain identifier for AlphaGEM contracts")
    p.add_argument("--prefix", default="PROT", help="Safe protein ID prefix (default: PROT)")
    p.add_argument("--keep-original-ids", action="store_true", help="Keep original transcript/protein IDs if already safe")
    p.add_argument("--method-label", help="Optional label recorded in standardization tables")
    p.add_argument(
        "--gffread",
        default="auto",
        help="Path to gffread binary, 'auto' (use gffread if installed else built-in Python splicer), or 'builtin'",
    )
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--tool-root", type=Path, help="Optional root containing tools/, envs/, and weights/ for annevo/helixer/tiberius")
    p.add_argument("--tool-python", help="Explicit Python interpreter for annevo/tiberius")
    p.add_argument("--annevo-dir", type=Path, help="Path to ANNEVO checkout containing prediction.py & decoding.py")
    p.add_argument("--tiberius-dir", type=Path, help="Path to Tiberius checkout containing tiberius.py")
    p.add_argument("--helixer-bin", help="Path to Helixer.py executable")
    p.add_argument("--funannotate-bin", help="Path to funannotate executable")
    p.add_argument("--masked-genome", type=Path, help="Soft-masked genome FASTA for funannotate")
    p.add_argument("--reference-proteins", type=Path, help="Reference protein FASTA (UniProt/RefSeq) for funannotate --protein_evidence")
    p.add_argument("--species", default="Fungal strain", help="Species label for funannotate")
    p.add_argument("--lineage", default="Fungi", help="Lineage label for ANNEVO (default: Fungi)")
    p.add_argument("--model-weights", type=Path, help="Explicit model weights/config path (.pt / .h5 / .yaml)")
    p.add_argument("--predictor-cmd", help="Command template for --tool custom (supports {genome} and {gff})")
    args = p.parse_args(argv)

    manifest = run_predict_proteome(args)
    print(
        json.dumps(
            {
                "status": manifest["status"],
                "tool": manifest["tool"],
                "target_proteins": manifest["target_proteins"],
                "proteome_fasta": manifest["proteome_fasta"],
                "proteome_sha256": manifest["proteome_sha256"],
                "inputs_dir": manifest["inputs_dir"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
