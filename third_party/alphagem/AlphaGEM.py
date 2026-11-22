#!/usr/bin/env python3
"""Command-line entry point for the repaired AlphaGEM paper workflow."""

import argparse
import sys
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
SRC_DIR = REPO_ROOT / "src"
PLMSEARCH_UTIL = REPO_ROOT / "plmsearchtools" / "plmsearch" / "plmsearch_util"
for module_path in (SRC_DIR, PLMSEARCH_UTIL):
    if str(module_path) not in sys.path:
        sys.path.insert(0, str(module_path))

from src.pipeline_config import (
    PROFILES,
    get_profile,
    override_profile,
    resolve_profile_for_reference,
)
from src.pipeline_runner import PipelineRun, PipelineRunner


def normalize_mode(value: str) -> str:
    normalized = value.strip().lower().replace("-", "_").replace(" ", "_")
    if normalized not in {"structure_alignment", "plmsearch"}:
        raise argparse.ArgumentTypeError(
            "mode must be structure_alignment, 'structure alignment', or plmsearch"
        )
    return normalized


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"invalid boolean value: {value!r}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the complete, stage-based AlphaGEM paper workflow."
    )
    parser.add_argument(
        "--profile",
        choices=sorted(PROFILES),
        default="paper",
        help="Versioned threshold profile. The published paper profile is the default.",
    )
    parser.add_argument(
        "--mode",
        type=normalize_mode,
        default="structure_alignment",
    )
    parser.add_argument(
        "--refname",
        choices=["yeast", "human", "ecoli", "strco", "synechocystis"],
        default="yeast",
    )
    parser.add_argument("--name", required=True, help="Target species/job name.")
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument(
        "--target-input-contract",
        type=Path,
        help=(
            "Versioned JSON pinning the exact target FASTA hashes, record "
            "count, normalized IDs, and provenance. Required for formal "
            "paper-profile runs."
        ),
    )
    parser.add_argument(
        "--gapfill-cds",
        type=Path,
        help=(
            "Legacy nucleotide-CDS proxy input. This is not accepted as the "
            "paper-exact raw-genome gap-fill input."
        ),
    )
    parser.add_argument(
        "--gapfill-genome",
        type=Path,
        help="Raw target-organism genome FASTA for paper tblastn gap filling.",
    )
    parser.add_argument(
        "--gapfill-coordinate-map",
        type=Path,
        help=(
            "Audited TSV mapping every --fasta protein to a raw-genome locus. "
            "Required with --gapfill-genome by the paper profile."
        ),
    )
    parser.add_argument(
        "--gapfill-coordinate-segments",
        type=Path,
        help=(
            "Audited one-row-per-CDS-segment TSV for every target protein. "
            "Required by the paper profile so tblastn HSPs are matched to "
            "exon/CDS unions rather than intron-spanning envelopes."
        ),
    )
    parser.add_argument(
        "--gapfill-coordinate-contract",
        type=Path,
        help=(
            "Versioned JSON pinning the raw genome, coordinate files, "
            "coordinate semantics, official sources, audit command, and "
            "audit-script/resource-manifest hashes."
        ),
    )
    parser.add_argument("--maplist", "--list", dest="maplist", type=Path)
    parser.add_argument("--structure", type=Path)

    parser.add_argument(
        "--clean",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Run CLEAN annotation (enabled for the complete paper workflow).",
    )
    parser.add_argument(
        "--deepec",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Run DeepECtransformer annotation.",
    )
    parser.add_argument(
        "--deepec-python",
        type=Path,
        help=(
            "Python executable from the pinned deepectransformer environment. "
            "Required for a formal paper DeepEC stage."
        ),
    )
    parser.add_argument(
        "--plm-annotation",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Run PLMSearch against Rhea for non-homologous proteins.",
    )
    parser.add_argument(
        "--cleanuse",
        dest="clean",
        type=parse_bool,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )

    parser.add_argument("--grothmedium", choices=["min", "full"], default="min")
    parser.add_argument("--essentialrxns", type=Path)
    parser.add_argument(
        "--solver",
        default="auto",
        help="COBRA solver name, or auto to preserve the environment default.",
    )
    parser.add_argument(
        "--blast-threads",
        type=int,
        default=1,
        help="Threads for each all-vs-all BLASTP direction.",
    )
    parser.add_argument(
        "--esm-model-path",
        default="esm1b_t33_650M_UR50S",
        help="ESM-1b model identifier or exact local .pt file.",
    )
    parser.add_argument(
        "--esm-token-budget",
        type=int,
        default=16384,
        help="Maximum ESM tokens per batch; does not change embedding semantics.",
    )
    parser.add_argument(
        "--prottrans-model-path",
        type=Path,
        help=(
            "Local Rostlab/prot_t5_xl_uniref50 snapshot. Required by the "
            "paper structure-homology stage; remote downloads are disabled."
        ),
    )
    parser.add_argument("--start-at", choices=PipelineRunner.STAGES)
    parser.add_argument("--stop-after", choices=PipelineRunner.STAGES)
    parser.add_argument(
        "--from-checkpoint",
        type=Path,
        help=(
            "Source attempt directory containing imported completed prefix stages. "
            "Requires "
            "--start-at; imported stages are recorded as "
            "imported_from_original_run, never as newly executed."
        ),
    )
    parser.add_argument(
        "--gapfill-tie-break",
        choices=["strict", "deterministic"],
        default="strict",
        help=(
            "Equal-highest tblastn hits mapping to different target proteins: "
            "'strict' refuses (the audited default); 'deterministic' picks "
            "the largest coordinate overlap in bp, then the smallest target "
            "protein ID, and records every broken tie in the gapfill audit."
        ),
    )
    parser.add_argument(
        "--gapfill-unmapped-hit-policy",
        choices=["error", "skip_with_audit"],
        default="error",
        help=(
            "Equal-top tblastn rows overlapping no predicted target "
            "coordinate: 'error' aborts (the audited default, author "
            "behavior); 'skip_with_audit' skips only those rows, records each "
            "skip in the gapfill audit, and a fully unmapped query then "
            "carries no gap evidence."
        ),
    )
    parser.add_argument(
        "--reaction-pool-already-present-policy",
        choices=[
            "author_skip",
            "merge_verified_identical_or",
            "merge_verified_identical_or_multisource",
            "merge_threeway_or_multisource",
        ],
        default="author_skip",
        help=(
            "Handling of candidates whose reaction ID is already in the "
            "model. Non-default values merge the candidate GPR only after "
            "machine-verified identity, and require the strict formal "
            "reaction-pool mode."
        ),
    )
    parser.add_argument(
        "--reaction-pool-conflict-policy",
        choices=["reject", "suffix_new_metabolite", "remap_by_model_identity"],
        default="reject",
        help=(
            "Handling of postmerge same-ID metabolite identity conflicts. "
            "Non-default values require the strict formal reaction-pool mode."
        ),
    )

    # Optional overrides retain CLI compatibility while named profiles remain immutable.
    parser.add_argument("--ss-threshold", type=float)
    parser.add_argument("--plm-coverage", type=float)
    parser.add_argument("--plm-identity", type=int)
    parser.add_argument("--bbh-coverage", type=float)
    parser.add_argument("--bbh-identity", type=int)
    parser.add_argument("--USalignfilter", type=float)
    parser.add_argument("--TMscore", type=float)
    parser.add_argument("--coverage", type=float)
    parser.add_argument("--pLDDT", type=float)
    parser.add_argument("--TMscoretrans", type=float)
    parser.add_argument("--coveragetrans", type=float)
    parser.add_argument("--upTMscore", type=float)
    parser.add_argument("--upcoverage", type=float)
    parser.add_argument("--esp", type=float)
    parser.add_argument("--prottrans_distance", "--prottrans-distance", type=float)
    parser.add_argument("--rhea-threshold", type=float)
    parser.add_argument("--ensemble-threshold", type=float)
    return parser


def load_essential_reactions(path):
    if path is None:
        return None
    import pandas as pd

    table = pd.read_excel(path)
    if "reaction_id" not in table.columns:
        raise ValueError(f"{path} must contain a reaction_id column")
    return table["reaction_id"].dropna().astype(str).tolist()


def build_run(args) -> PipelineRun:
    profile = override_profile(
        resolve_profile_for_reference(
            get_profile(args.profile),
            args.refname,
        ),
        ss_predictor=args.ss_threshold,
        plm_coverage=args.plm_coverage,
        plm_identity=args.plm_identity,
        bbh_coverage=args.bbh_coverage,
        bbh_identity=args.bbh_identity,
        usalign=args.USalignfilter,
        tm_score=args.TMscore,
        structure_coverage=args.coverage,
        plddt=args.pLDDT,
        transporter_tm_score=args.TMscoretrans,
        transporter_coverage=args.coveragetrans,
        high_confidence_tm_score=args.upTMscore,
        high_confidence_coverage=args.upcoverage,
        dbscan_eps=args.esp,
        prottrans_distance=args.prottrans_distance,
        rhea_ss_predictor=args.rhea_threshold,
        ensemble_score=args.ensemble_threshold,
    )
    return PipelineRun(
        name=args.name,
        refname=args.refname,
        fasta=args.fasta,
        mode=args.mode,
        profile=profile,
        target_input_contract=args.target_input_contract,
        gapfill_cds=args.gapfill_cds,
        gapfill_genome=args.gapfill_genome,
        gapfill_coordinate_map=args.gapfill_coordinate_map,
        gapfill_coordinate_segments=args.gapfill_coordinate_segments,
        gapfill_coordinate_contract=args.gapfill_coordinate_contract,
        maplist=args.maplist,
        structure=args.structure,
        clean=args.clean,
        deepec=args.deepec,
        deepec_python=args.deepec_python,
        plm_annotation=args.plm_annotation,
        growth_medium=args.grothmedium,
        essential_reactions=load_essential_reactions(args.essentialrxns),
        solver=args.solver,
        blast_threads=args.blast_threads,
        esm_model_path=args.esm_model_path,
        esm_token_budget=args.esm_token_budget,
        prottrans_model_path=args.prottrans_model_path,
        reaction_pool_already_present_policy=(
            args.reaction_pool_already_present_policy
        ),
        reaction_pool_conflict_policy=args.reaction_pool_conflict_policy,
        from_checkpoint=args.from_checkpoint,
        gapfill_tie_break=args.gapfill_tie_break,
        gapfill_unmapped_hit_policy=args.gapfill_unmapped_hit_policy,
        start_at=args.start_at,
        stop_after=args.stop_after,
    )


def main(argv=None) -> int:
    started = time.time()
    args = build_parser().parse_args(argv)
    run = build_run(args)
    completed = PipelineRunner(REPO_ROOT).run(run)
    elapsed = time.time() - started
    print(f"Completed stages: {', '.join(completed)}")
    print(f"AlphaGEM costs {elapsed:.2f} seconds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
