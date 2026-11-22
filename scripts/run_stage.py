#!/usr/bin/env python3
"""Run one metabolic reconstruction or target-screening stage against prepared inputs."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def _parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--runtime",
        default=os.environ.get("G2F_RUNTIME"),
        help="prepared writable AlphaGEM runtime (or G2F_RUNTIME); required for localize/gapfill",
    )
    common.add_argument(
        "--workdir",
        default=os.environ.get("G2F_STAGE"),
        help="stage inputs/results directory (or G2F_STAGE)",
    )
    common.add_argument(
        "--inputs",
        default="/work/inputs",
        help="prepared input directory",
    )
    common.add_argument(
        "--name",
        default="sample",
        help="AlphaGEM reconstruction name",
    )
    parser = argparse.ArgumentParser(description=__doc__)
    stages = parser.add_subparsers(dest="stage", required=True)
    for stage in ("localize", "gapfill", "target", "fseof", "paralogs", "knockout", "filter"):
        command = stages.add_parser(stage, parents=[common])
        if stage == "localize":
            command.add_argument(
                "--localization-table",
                type=Path,
                help="reuse an explicit DeepLoc CSV or localization TSV",
            )
            command.add_argument(
                "--deeploc-bin",
                default=None,
                help="explicit DeepLoc binary/command (default: auto-detect G2F_DEEPLOC_BIN -> tools/deeploc/ -> PATH, else fallback to default compartments)",
            )
            command.add_argument(
                "--require-deeploc",
                action="store_true",
                help="fail instead of falling back to default compartments when DeepLoc is not installed",
            )
            command.add_argument(
                "--draft-model",
                type=Path,
                help="override draft model path (default: <runtime>/working/<name>/<name>-GEM_withgaps.yml)",
            )
            command.add_argument(
                "--template-model",
                type=Path,
                help="override template model path (default: <runtime>/models/yeast-GEM.xml)",
            )
        if stage == "target":
            command.add_argument(
                "--source-model",
                type=Path,
                help="input metabolic model (default: <workdir>/gapfilled.xml)",
            )
            command.add_argument(
                "--target-spec",
                help="JSON path defining custom metabolites/reactions (omit when using --target-reaction)",
            )
            command.add_argument(
                "--target-reaction",
                help="existing model reaction ID (or override target reaction in --target-spec)",
            )
        if stage in {"fseof", "filter"}:
            command.add_argument("--growth-reaction", default="r_4041")
        if stage in {"fseof", "knockout"}:
            command.add_argument("--target-reaction", default=None)
        if stage == "paralogs":
            command.add_argument("--diamond-bin", default="diamond")
            command.add_argument("--threads", type=int, default=8)
            command.add_argument(
                "--self-hits",
                type=Path,
                help="reuse precomputed 8-column DIAMOND self-alignment TSV",
            )
        if stage == "filter":
            command.add_argument(
                "--essential-genes-tsv",
                type=Path,
                help="precomputed growth-KO TSV or single-column essential gene list; if omitted, computed via --growth-reaction",
            )
            command.add_argument("--essential-threshold", type=float, default=0.5)
    return parser


def _working(args: argparse.Namespace) -> Path:
    return args.runtime / "working" / args.name


def _check_paths(
    parser: argparse.ArgumentParser, inputs: list[Path], outputs: list[Path]
) -> None:
    for path in inputs:
        if not path.is_file():
            parser.error(f"required prepared input is missing: {path}")
    for path in outputs:
        if path.exists() or path.is_symlink():
            parser.error(f"output already exists; choose a fresh stage directory: {path}")


def _write_json(path: Path, value: object) -> None:
    with path.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(value, indent=2, allow_nan=False) + "\n")


def _configure_highs() -> None:
    sys.path.insert(0, str(ROOT / "third_party" / "alphagem" / "src"))
    import cobra
    from optlang_highs_interface import register_highs_interface

    register_highs_interface()
    cobra.Configuration().solver = "highs"


def _product_model(args: argparse.Namespace):
    _configure_highs()
    import cobra

    model = cobra.io.read_sbml_model(str(args.workdir / "product.xml"))
    model.solver.configuration.tolerances.feasibility = 1e-9
    model.solver.configuration.continuous_lp_backend = "linprog"
    return model


def _resolve_stage_target_reaction(args: argparse.Namespace) -> str:
    if getattr(args, "target_reaction", None):
        return str(args.target_reaction)
    receipt_path = args.workdir / "target_receipt.json"
    if receipt_path.is_file():
        data = json.loads(receipt_path.read_text(encoding="utf-8"))
        if data.get("product_reaction"):
            return str(data["product_reaction"])
    raise ValueError(
        "No target reaction is configured; run the target stage first or pass --target-reaction"
    )


def _localize(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    if args.draft_model is None and args.runtime is None:
        parser.error("localize requires --runtime (or G2F_RUNTIME) unless --draft-model is supplied")
    if args.template_model is None and args.runtime is None:
        parser.error("localize requires --runtime (or G2F_RUNTIME) unless --template-model is supplied")
    draft = args.draft_model.resolve() if args.draft_model else (_working(args) / f"{args.name}-GEM_withgaps.yml")
    reference = args.template_model.resolve() if args.template_model else (args.runtime / "models" / "yeast-GEM.xml")
    table = args.localization_table
    import cobra
    from geno2flux.localization import (
        _load_model,
        load_localization_table,
        load_models_for_refinement,
        refine_compartments,
        resolve_deeploc_command,
    )

    resolved_cmd: list[str] | None = None
    loc_source = "localization_table"
    if table is None:
        try:
            resolved_cmd, loc_source = resolve_deeploc_command(args.deeploc_bin)
        except FileNotFoundError as exc:
            parser.error(str(exc))
        if resolved_cmd is None and args.require_deeploc:
            parser.error(
                "DeepLoc not found (checked --deeploc-bin, G2F_DEEPLOC_BIN, tools/deeploc/, and PATH) "
                "and --require-deeploc was set"
            )

    required = [draft, reference]
    if table is not None:
        required.append(table)
    elif resolved_cmd is not None:
        required.append(args.inputs / "proteome.fasta")
    outputs = [args.workdir / "localized.yml", args.workdir / "localization_audit.json"]
    prediction_dir = args.workdir / "deeploc"
    _check_paths(parser, required, outputs + ([prediction_dir] if resolved_cmd is not None else []))

    args.workdir.mkdir(parents=True, exist_ok=True)
    if table is None and resolved_cmd is None:
        print(
            "[geno2flux] DeepLoc not detected; falling back to default model compartments.",
            file=sys.stderr,
        )
        model = _load_model(draft)
        audit = {
            "pulled_reaction_count": 0,
            "pulled_reactions": [],
            "prune_relocated_genes": False,
            "final_reaction_count": len(model.reactions),
            "deeploc_status": "fallback_default_compartments",
            "deeploc_source": "not_found_fallback_default",
            "fallback": "deeploc_unavailable_default_compartments",
        }
        cobra.io.save_yaml_model(model, str(outputs[0]))
        _write_json(outputs[1], audit)
        return

    if table is None:
        assert resolved_cmd is not None
        prediction_dir.mkdir()
        subprocess.run(
            [*resolved_cmd, "-f", str(args.inputs / "proteome.fasta"), "-o", str(prediction_dir)],
            check=True,
        )
        tables = sorted(prediction_dir.glob("results_*.csv"))
        if not tables:
            raise RuntimeError("deeploc2 produced no results_*.csv")
        table = tables[-1]
    model, template = load_models_for_refinement(draft, reference)
    audit = refine_compartments(model, template, load_localization_table(table))
    audit["deeploc_status"] = "reused_table" if args.localization_table is not None else "predicted_with_deeploc"
    audit["deeploc_source"] = loc_source
    cobra.io.save_yaml_model(model, str(outputs[0]))
    _write_json(outputs[1], audit)


def _gapfill(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    if args.runtime is None:
        parser.error("gapfill requires --runtime or G2F_RUNTIME")
    localized = args.workdir / "localized.yml"
    genome = args.inputs / "genome.fasta"
    coordinate_map = args.inputs / "coordinate-map.tsv"
    segments = args.inputs / "coordinate-segments.tsv"
    contract = args.inputs / "gapfill-coordinate-contract.json"
    _check_paths(
        parser,
        [localized, genome, coordinate_map, segments, contract],
        [args.workdir / "gapfilled.xml", args.workdir / "gapfill_audit.json"],
    )
    _configure_highs()
    from geno2flux.highspy_backend import install_highspy_backend

    install_highspy_backend()
    sys.path[:0] = [str(args.runtime), str(args.runtime / "src")]
    from src.gapfilling import gapfill

    shutil.copy2(localized, _working(args) / f"{args.name}-GEM_withgaps.yml")
    previous = Path.cwd()
    try:
        os.chdir(args.runtime)
        from geno2flux.pipeline import _tool_path_scope

        with _tool_path_scope():
            gapfill(
                args.name,
                "yeast",
                "yeast-GEM.xml",
                grothmedium="min",
                search_backend="tblastn_genome",
                gapfill_genome=str(genome),
                gapfill_coordinate_map=str(coordinate_map),
                gapfill_coordinate_segments=str(segments),
                gapfill_coordinate_contract=str(contract),
                strategy="selected_medium",
                bitscore_threshold=200,
                bitscore_operator="gt",
                tie_break="deterministic",
                unmapped_hit_policy="skip_with_audit",
                require_solver_tolerances=True,
                require_optimal_solutions=True,
            )
    finally:
        os.chdir(previous)
    shutil.copy2(_working(args) / f"{args.name}-GEM.xml", args.workdir / "gapfilled.xml")
    shutil.copy2(_working(args) / "gapfill_audit.json", args.workdir / "gapfill_audit.json")


def _target(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    source = args.source_model.resolve() if args.source_model else (args.workdir / "gapfilled.xml")
    destination = args.workdir / "product.xml"
    receipt_path = args.workdir / "target_receipt.json"
    _check_paths(parser, [source], [destination, receipt_path])
    _configure_highs()
    import cobra
    from geno2flux.target_module import mount_target_module

    args.workdir.mkdir(parents=True, exist_ok=True)
    suffix = source.suffix.lower()
    if suffix in {".yml", ".yaml"}:
        model = cobra.io.load_yaml_model(str(source))
    elif suffix == ".json":
        model = cobra.io.load_json_model(str(source))
    else:
        model = cobra.io.read_sbml_model(str(source))
    receipt = mount_target_module(
        model,
        args.target_spec,
        target_reaction=args.target_reaction,
        base_dir=Path.cwd(),
    )
    cobra.io.write_sbml_model(model, str(destination))
    _write_json(receipt_path, receipt)


def _fseof(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    fasta = args.inputs / "proteome.fasta"
    required = [args.workdir / "product.xml"]
    if fasta.is_file():
        required.append(fasta)
    outputs = [
        args.workdir / "fseof.json",
        args.workdir / "all_reactions_fseof.tsv",
    ]
    if fasta.is_file():
        outputs.append(args.workdir / "all_genes_fseof.tsv")
    _check_paths(parser, required, outputs)
    from geno2flux.pipeline import _gene_report, _load_protein_map, _write_tsv
    from geno2flux.published_fseof import run_published_fseof

    model = _product_model(args)
    target_rxn = _resolve_stage_target_reaction(args)
    result = run_published_fseof(
        model,
        args.growth_reaction,
        target_rxn,
        config={"iterations": 10},
    )
    _write_json(args.workdir / "fseof.json", result)
    _write_tsv(
        args.workdir / "all_reactions_fseof.tsv",
        result["reaction_rows"],
        [
            "reaction_id",
            "reaction_name",
            "gpr",
            "gene_ids",
            "boundary",
            "reversible",
            "raven_reported_direction",
            "candidate",
            "direction",
            "raw_candidate",
            "robust_candidate",
            "slope",
            "initial_flux",
            "final_flux",
        ],
    )
    if fasta.is_file():
        id_map = args.inputs / "id_map.tsv"
        proteins = _load_protein_map(id_map if id_map.is_file() else None, fasta)
        gene_rows = _gene_report(
            proteins, result["reaction_rows"], {g.id for g in model.genes}
        )
        _write_tsv(
            args.workdir / "all_genes_fseof.tsv",
            gene_rows,
            [
                "gene_id",
                "in_model",
                "raw_score",
                "robust_score",
                "raw_candidate_reactions",
                "robust_candidate_reactions",
                "all_associated_reactions",
                "raw_rank",
                "robust_rank",
            ],
        )


def _paralogs(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    fasta = args.inputs / "proteome.fasta"
    outputs = [
        args.workdir / "self_hits.tsv",
        args.workdir / "families.json",
        args.workdir / "families.tsv",
    ]
    if args.self_hits is None:
        outputs.append(args.workdir / "diamond")
    required = [fasta, args.workdir / "fseof.json"]
    if args.self_hits is not None:
        required.append(args.self_hits.expanduser().resolve())
    _check_paths(parser, required, outputs)
    from Bio import SeqIO
    from geno2flux.paralog_grouping import (
        TIERS,
        build_families,
        load_self_hits,
        protein_key,
        rank_paralog_families,
        run_diamond_self_alignment,
        write_family_tables,
    )
    from geno2flux.pipeline import resolve_tool_executable

    reactions = json.loads((args.workdir / "fseof.json").read_text(encoding="utf-8"))["reaction_rows"]
    scores = {}
    for protein in SeqIO.parse(fasta, "fasta"):
        gene = protein_key(protein.id)
        slopes = [
            row["slope"]
            for row in reactions
            if row["raw_candidate"] and gene in row["gene_ids"]
        ]
        scores[gene] = sum(slopes) / len(slopes) if slopes else None
    if args.self_hits is not None:
        shutil.copy2(args.self_hits.expanduser().resolve(), args.workdir / "self_hits.tsv")
        hits_path = args.workdir / "self_hits.tsv"
    else:
        hits_path = run_diamond_self_alignment(
            fasta,
            args.workdir / "self_hits.tsv",
            diamond_bin=resolve_tool_executable(args.diamond_bin),
            threads=args.threads,
            workdir=args.workdir / "diamond",
        )
    families = build_families(load_self_hits(hits_path), TIERS["default"])
    ranked = rank_paralog_families(scores, families, all_protein_ids=scores)
    write_family_tables(ranked, args.workdir / "families.json", args.workdir / "families.tsv")


def _knockout(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    _check_paths(
        parser,
        [args.workdir / "product.xml", args.workdir / "families.json"],
        [args.workdir / "knockout"],
    )
    from geno2flux.knockout import run_ko_suite

    target_rxn = _resolve_stage_target_reaction(args)
    rows = json.loads((args.workdir / "families.json").read_text(encoding="utf-8"))["families"]
    families = {row["family_id"]: row["members"] for row in rows}
    summary = run_ko_suite(
        _product_model(args),
        target_rxn,
        families=families,
        out_dir=args.workdir / "knockout",
    )
    _write_json(args.workdir / "knockout" / "ko_summary.json", summary)


def _filter(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    ko_dir = args.workdir / "knockout"
    required = [
        args.workdir / "product.xml",
        args.workdir / "all_genes_fseof.tsv",
        args.workdir / "families.tsv",
        ko_dir / "single_gene_ko.tsv",
        ko_dir / "family_ko.tsv",
    ]
    if args.essential_genes_tsv is not None:
        required.append(args.essential_genes_tsv.expanduser().resolve())
    outputs = [
        args.workdir / "all_genes_fseof_essential_filtered.tsv",
        args.workdir / "families_essential_filtered.tsv",
        ko_dir / "single_gene_ko_essential_filtered.tsv",
        ko_dir / "family_ko_essential_filtered.tsv",
        ko_dir / "essential_filter_detail.json",
    ]
    _check_paths(parser, required, outputs)
    from geno2flux.essential_filter import (
        FAMILY_FIELDS,
        FSEOF_GENE_FIELDS,
        KO_GENE_FIELDS,
        apply_essential_filter,
        essential_family_ids,
        load_essential_gene_ids,
        positional_filtered_rank,
        read_tsv_rows,
        validate_threshold,
        write_filter_manifest_json,
        write_filtered_tsv,
    )
    from geno2flux.knockout import run_single_gene_ko, write_ko_tsv

    threshold = validate_threshold(args.essential_threshold)
    if args.essential_genes_tsv is not None:
        essential_tsv = args.essential_genes_tsv.expanduser().resolve()
    else:
        model = _product_model(args)
        growth_opt, growth_nonopt = run_single_gene_ko(model, args.growth_reaction, ko_dir)
        essential_tsv = ko_dir / "single_gene_ko_growth.tsv"
        write_ko_tsv(essential_tsv, growth_opt, "gene_id")
        write_ko_tsv(
            ko_dir / "single_gene_ko_growth_nonoptimal.tsv",
            growth_nonopt,
            "gene_id",
            nonoptimal=True,
        )

    record = load_essential_gene_ids(essential_tsv, threshold)
    essential_genes = record["ids"]
    _, gene_rows = read_tsv_rows(args.workdir / "all_genes_fseof.tsv")
    _, family_rows = read_tsv_rows(args.workdir / "families.tsv")
    essential_fams = essential_family_ids(family_rows, essential_genes)

    gres = apply_essential_filter(gene_rows, essential_genes, "gene_id", "raw_rank", order="gene_v2")
    write_filtered_tsv(
        args.workdir / "all_genes_fseof_essential_filtered.tsv",
        gres["kept"],
        FSEOF_GENE_FIELDS,
    )
    fres = apply_essential_filter(family_rows, essential_fams, "family_id", "max_score_rank")
    positional_filtered_rank(fres["kept"], "max_score_rank", "family_id", "filtered_rank_max")
    positional_filtered_rank(fres["kept"], "mean_score_rank", "family_id", "filtered_rank_mean")
    write_filtered_tsv(
        args.workdir / "families_essential_filtered.tsv",
        fres["kept"],
        FAMILY_FIELDS + ["filtered_rank_max", "filtered_rank_mean"],
    )
    gene_ko_fields, gene_ko_rows = read_tsv_rows(ko_dir / "single_gene_ko.tsv")
    kg = apply_essential_filter(gene_ko_rows, essential_genes, "gene_id", "ko_rank")
    write_filtered_tsv(
        ko_dir / "single_gene_ko_essential_filtered.tsv",
        kg["kept"],
        gene_ko_fields,
    )
    fam_ko_fields, fam_ko_rows = read_tsv_rows(ko_dir / "family_ko.tsv")
    kf = apply_essential_filter(fam_ko_rows, essential_fams, "family_id", "ko_rank")
    write_filtered_tsv(
        ko_dir / "family_ko_essential_filtered.tsv",
        kf["kept"],
        fam_ko_fields,
    )
    growth_fields, growth_rows = read_tsv_rows(essential_tsv)
    gg = apply_essential_filter(
        [r for r in growth_rows if r.get("status") == "optimal"],
        essential_genes,
        "gene_id",
        "ko_rank",
    )
    write_filtered_tsv(
        ko_dir / "single_gene_ko_growth_essential_filtered.tsv",
        gg["kept"],
        growth_fields if record["mode"] != "explicit_id_list" else KO_GENE_FIELDS,
    )
    write_filter_manifest_json(
        ko_dir / "essential_filter_detail.json",
        {
            "essentiality_source": record["mode"],
            "essentiality_input": str(essential_tsv),
            "threshold": record["threshold"],
            "essential_genes": len(essential_genes),
            "essential_families": len(essential_fams),
            "fseof_genes_kept": len(gres["kept"]),
            "families_kept": len(fres["kept"]),
            "gene_ko_kept": len(kg["kept"]),
            "growth_ko_kept": len(gg["kept"]),
        },
    )


def main(argv: list[str] | None = None) -> None:
    parser = _parser()
    args = parser.parse_args(argv)
    if not args.workdir:
        parser.error("--workdir or G2F_STAGE is required")
    args.workdir = Path(args.workdir).expanduser().resolve()
    args.inputs = Path(args.inputs).expanduser().resolve()
    if args.runtime:
        args.runtime = Path(args.runtime).expanduser().resolve()
    if args.stage == "localize" and args.localization_table is not None:
        args.localization_table = args.localization_table.expanduser().resolve()
    if Path(args.name).name != args.name or args.name in {"", ".", ".."}:
        parser.error("--name must be a reconstruction name, not a path")
    sys.path.insert(0, str(ROOT / "src"))
    handlers = {
        "localize": _localize,
        "gapfill": _gapfill,
        "target": _target,
        "fseof": _fseof,
        "paralogs": _paralogs,
        "knockout": _knockout,
        "filter": _filter,
    }
    handlers[args.stage](args, parser)


if __name__ == "__main__":
    main()
