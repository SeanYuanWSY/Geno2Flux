"""Auditable, stage-based runner for the AlphaGEM paper workflow."""

import ast
import gc
import importlib
import importlib.metadata
import hashlib
import json
import math
import os
import pickle
import re
import sqlite3
import subprocess
import sys
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, List, Optional

from .pipeline_config import (
    PipelineProfile,
    get_profile,
    resolve_profile_for_reference,
)
from .pipeline_contracts import (
    assert_embedding_provenance,
    resolve_eggnog_diamond_db,
    validate_gapfill_cds,
    validate_gapfill_coordinate_contract,
    validate_gapfill_genome,
    validate_job_name,
    validate_prott5_xl_uniref50_snapshot,
    validate_target_input_contract,
)
from .paper_tool_contracts import (
    PAPER_EGGNOG_DATABASE_CONTRACT,
    PAPER_EGGNOG_TAXON_DATABASE_CONTRACTS,
    PAPER_METHOD_EVIDENCE_CONTRACT,
    PAPER_RHEA_ARTIFACT_CONTRACT,
    PAPER_RESOURCE_ARCHIVE_SHA256,
    PAPER_RUNTIME_TOOL_BUNDLE_CONTRACTS,
)


REFERENCE_MODELS = {
    "ecoli": "ecoli-GEM.xml",
    "yeast": "yeast-GEM.xml",
    "strco": "Sco-GEM.xml",
    "human": "Human-GEM.xml",
    "synechocystis": "iSynCJ816.xml",
}

# This is the reaction-pool snapshot shipped with the bundled paper source.
# Version 3 also records pre-addition chemical-identity conflicts separately.
# Version 4 adds the already-present GPR-merge and identity-conflict
# resolution policies and recomputes their recorded effects.
FROZEN_REACTION_POOL_SHA256 = (
    "fbe471dd7d76e947aa2a31af7050b3020a179a3ee64bcd89e21b7a7fe894c611"
)
REACTION_POOL_AUDIT_POLICY_VERSION = 4

# Restated from add_reactions_based_pool so the validator pins the vocabulary
# independently of the builder it audits.
REACTION_POOL_ALREADY_PRESENT_POLICIES = (
    "author_skip",
    "merge_verified_identical_or",
    "merge_verified_identical_or_multisource",
    "merge_threeway_or_multisource",
)
REACTION_POOL_CONFLICT_POLICIES = (
    "reject",
    "suffix_new_metabolite",
    "remap_by_model_identity",
)
REACTION_POOL_MULTISOURCE_MIN_SOURCES = 2
REACTION_POOL_MULTISOURCE_COLUMNS = ("anno_A", "anno_B", "anno_C", "anno_D")


def _rhea_number_key(value):
    """Normalize a spreadsheet rhea cell to the candidate ID's number string."""
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        if not math.isfinite(float(value)) or int(value) != value:
            return None
        return str(int(value))
    text = str(value).strip()
    return text or None


def _author_name_key(name):
    """Mirror of the builder's author metdict name equivalence key.

    Must stay byte-identical in behavior to
    ``add_reactions_based_pool._author_name_key`` (lowercase, parentheses
    stripped); the merge validator and the builder share this relation.
    """
    return (name or "").lower().replace("(", "").replace(")", "")

PROTTRANS_REFERENCE_MODELS = {
    "ecoli": "iML1515.json",
    "yeast": "yeast-GEM.xml",
    "strco": "Sco-GEM.xml",
    "human": "Human-GEM.xml",
    "synechocystis": "iSynCJ816.xml",
}


@dataclass
class PipelineRun:
    name: str
    refname: str
    fasta: Path
    mode: str
    profile: PipelineProfile
    target_input_contract: Optional[Path] = None
    gapfill_cds: Optional[Path] = None
    gapfill_genome: Optional[Path] = None
    gapfill_coordinate_map: Optional[Path] = None
    gapfill_coordinate_segments: Optional[Path] = None
    gapfill_coordinate_contract: Optional[Path] = None
    maplist: Optional[Path] = None
    structure: Optional[Path] = None
    clean: bool = True
    deepec: bool = True
    deepec_python: Optional[Path] = None
    plm_annotation: bool = True
    growth_medium: str = "min"
    essential_reactions: Optional[List[str]] = None
    solver: str = "auto"
    blast_threads: int = 1
    esm_model_path: str = "esm1b_t33_650M_UR50S"
    esm_token_budget: int = 16384
    prottrans_model_path: Optional[Path] = None
    reaction_pool_already_present_policy: str = "author_skip"
    reaction_pool_conflict_policy: str = "reject"
    from_checkpoint: Optional[Path] = None
    gapfill_tie_break: str = "strict"
    gapfill_unmapped_hit_policy: str = "error"
    start_at: Optional[str] = None
    stop_after: Optional[str] = None


class PipelineRunner:
    STAGES = (
        "prepare",
        "homology",
        "draft_model",
        "prepare_nonhomologs",
        "eggnog",
        "clean",
        "deepec",
        "rhea",
        "merge_annotations",
        "reaction_pool",
        "gapfill",
        "validate",
    )
    MANIFEST_SCHEMA_VERSION = 17
    RUNTIME_DISTRIBUTIONS = {
        "biopython": "biopython",
        "cobra": "cobra",
        "fair-esm": "fair-esm",
        "fuzzywuzzy": "fuzzywuzzy",
        "joblib": "joblib",
        "levenshtein": "Levenshtein",
        "numpy": "numpy",
        "openpyxl": "openpyxl",
        "optlang": "optlang",
        "pandas": "pandas",
        "rdflib": "rdflib",
        "requests": "requests",
        "scikit-learn": "scikit-learn",
        "scipy": "scipy",
        "sentencepiece": "sentencepiece",
        "sympy": "sympy",
        "torch": "torch",
        "tqdm": "tqdm",
        "transformers": "transformers",
        "xlsxwriter": "XlsxWriter",
    }

    def __init__(
        self,
        repo_root: Path,
        module_loader: Callable[[str], object] = importlib.import_module,
        verify_outputs: bool = True,
    ):
        self.repo_root = Path(repo_root).resolve()
        self.module_loader = module_loader
        self.verify_outputs = verify_outputs
        self.preflight_by_stage = {}
        self.resource_bundles_by_stage = {}
        self.tool_versions_by_stage = {}
        self._tool_version_cache = {}
        self.solver_record = None
        self._run_source_identity = None
        self._run_runtime_environment = None
        self._checkpoint_import_record = None

    @classmethod
    def selected_stages(
        cls,
        start_at: Optional[str] = None,
        stop_after: Optional[str] = None,
    ) -> List[str]:
        if start_at is not None and start_at not in cls.STAGES:
            raise ValueError(f"Unknown start stage: {start_at}")
        if stop_after is not None and stop_after not in cls.STAGES:
            raise ValueError(f"Unknown stop stage: {stop_after}")
        first = cls.STAGES.index(start_at) if start_at else 0
        last = cls.STAGES.index(stop_after) + 1 if stop_after else len(cls.STAGES)
        if first >= last:
            raise ValueError("start_at must not come after stop_after")
        return list(cls.STAGES[first:last])

    def run(self, run: PipelineRun) -> List[str]:
        with self._formal_no_bytecode(run):
            return self._run_pipeline(run)

    def _run_pipeline(self, run: PipelineRun) -> List[str]:
        stages = self.selected_stages(run.start_at, run.stop_after)
        if run.from_checkpoint is not None and (
            not run.start_at or run.start_at == self.STAGES[0]
        ):
            raise RuntimeError(
                "--from-checkpoint requires --start-at a stage after "
                f"{self.STAGES[0]!r}; a checkpoint only carries completed "
                "prefix stages"
            )
        self._validate_run(run, stages)
        self._assert_formal_run_identity(run)
        if self.verify_outputs:
            self._preflight_resources(run, stages)
        self._assert_formal_run_identity(run)
        resume_manifest = None
        if run.start_at and run.start_at != self.STAGES[0]:
            resume_manifest = self._validate_resume(run)
        self._assert_formal_run_identity(run)
        completed = list(
            resume_manifest.get("completed_stages", [])
            if resume_manifest
            else []
        )
        stage_records = dict(
            resume_manifest.get("stages", {}) if resume_manifest else {}
        )
        with self._working_directory():
            self._configure_solver(run.solver)
            self._assert_formal_run_identity(run)
            for stage in stages:
                self._assert_formal_run_identity(run)
                if self.verify_outputs:
                    self._assert_resources_for_stage(run, stage)
                    self._record_paper_bundles_for_stage(run, stage)
                getattr(self, f"_stage_{stage}")(run)
                self._assert_formal_run_identity(run)
                if self.verify_outputs:
                    self._assert_resources_for_stage(run, stage)
                    self._record_paper_bundles_for_stage(run, stage)
                    artifacts = self._verify_stage_outputs(stage, run)
                else:
                    artifacts = {}
                stage_records[stage] = {
                    "resources": self.preflight_by_stage.get(stage, {}),
                    "resource_bundles": self.resource_bundles_by_stage.get(
                        stage,
                        {},
                    ),
                    "tool_versions": self.tool_versions_by_stage.get(stage, {}),
                    "artifacts": artifacts,
                }
                completed.append(stage)
                self._write_manifest(run, completed, stage_records)
        return completed

    @contextmanager
    def _formal_no_bytecode(self, run: PipelineRun):
        if not (
            self.verify_outputs and self._is_exact_paper_profile(run)
        ):
            yield
            return
        previous_environment = os.environ.get("PYTHONDONTWRITEBYTECODE")
        previous_flag = sys.dont_write_bytecode
        os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
        sys.dont_write_bytecode = True
        try:
            yield
        finally:
            sys.dont_write_bytecode = previous_flag
            if previous_environment is None:
                os.environ.pop("PYTHONDONTWRITEBYTECODE", None)
            else:
                os.environ[
                    "PYTHONDONTWRITEBYTECODE"
                ] = previous_environment

    def _preflight_resources(self, run: PipelineRun, stages: List[str]) -> None:
        self.preflight_by_stage = {}
        self.resource_bundles_by_stage = {}
        self.tool_versions_by_stage = {}
        self.solver_record = None
        for stage in stages:
            self._record_resources_for_stage(run, stage)
            bundle_specs = self._required_resource_bundles_for_stage(
                run,
                stage,
            )
            self.resource_bundles_by_stage[stage] = (
                self._paper_bundle_records(bundle_specs)
            )
            tool_specs = self._required_tools_for_stage(run, stage)
            observed = {}
            mismatches = {}
            for tool, (command, expected) in tool_specs.items():
                cache_key = tuple(command)
                if cache_key not in self._tool_version_cache:
                    self._tool_version_cache[cache_key] = self._tool_version(command)
                version = self._tool_version_cache[cache_key]
                observed[tool] = version
                if not self._tool_version_matches(
                    tool,
                    expected,
                    version,
                ):
                    mismatches[tool] = (expected, version)
            if mismatches:
                raise RuntimeError(
                    f"Tool version preflight failed for stage {stage!r}: "
                    + json.dumps(mismatches, sort_keys=True)
                )
            self.tool_versions_by_stage[stage] = observed

        if any(
            stage in {"draft_model", "reaction_pool", "gapfill", "validate"}
            for stage in stages
        ):
            cobra = self._module("cobra")
            if run.solver != "auto":
                self._apply_solver_selection(cobra, run.solver)
            smoke = cobra.Model("alphagem_solver_preflight")
            interface = smoke.problem
            binary = interface.Variable(
                "alphagem_milp_smoke_binary",
                type="binary",
            )
            smoke.add_cons_vars(
                [
                    binary,
                    interface.Constraint(
                        binary,
                        lb=1,
                        name="alphagem_milp_smoke_constraint",
                    ),
                ]
            )
            smoke.objective = interface.Objective(binary, direction="min")
            solution = smoke.optimize()
            binary_primal = binary.primal
            if (
                solution.status != "optimal"
                or binary_primal is None
                or abs(float(binary_primal) - 1.0) > 1e-7
            ):
                raise RuntimeError(
                    f"Solver preflight failed for {run.solver}: "
                    f"status={solution.status}, binary={binary_primal}"
                )
            interface_name = smoke.solver.interface.__name__
            backend_packages = {
                "optlang.glpk_interface": "swiglpk",
                "optlang.gurobi_interface": "gurobipy",
                "optlang.cplex_interface": "cplex",
                "optlang.scipy_interface": "scipy",
                # The bundled HiGHS interface delegates to scipy.optimize.milp;
                # scipy's version is the backend version evidence (SciPy does
                # not expose the bundled HiGHS build's own version string).
                "optlang_highs_interface": "scipy",
            }
            package_versions = {}
            for package in (
                "cobra",
                "optlang",
                backend_packages.get(interface_name),
            ):
                if package is None:
                    continue
                try:
                    package_versions[package] = importlib.metadata.version(
                        package
                    )
                except importlib.metadata.PackageNotFoundError:
                    package_versions[package] = "unknown"
            solver_tolerances = {}
            for tolerance_name in (
                "feasibility",
                "integrality",
                "optimality",
            ):
                try:
                    solver_tolerances[tolerance_name] = float(
                        getattr(
                            smoke.solver.configuration.tolerances,
                            tolerance_name,
                        )
                    )
                except (AttributeError, TypeError, ValueError):
                    solver_tolerances[tolerance_name] = None
            self.solver_record = {
                "requested": run.solver,
                "interface": interface_name,
                "backend_evidence": "observed_runtime",
                "paper_specified_backend": False,
                "package_versions": package_versions,
                "tolerances": solver_tolerances,
                "milp_smoke": {
                    "status": solution.status,
                    "objective": float(solution.objective_value),
                    "binary_primal": float(binary_primal),
                },
            }
            highs_receipt = getattr(
                getattr(smoke.solver, "problem", None), "last_result_receipt",
                None,
            )
            if highs_receipt is not None:
                self.solver_record["milp_smoke"]["backend_receipt"] = (
                    highs_receipt
                )

    def _required_resources_for_stage(self, run: PipelineRun, stage: str):
        data_dir = self.repo_root / "data_available"
        eggnog_data = self.repo_root / "tools" / "eggnog_mapper" / "data"
        reference_model = self.repo_root / "models" / REFERENCE_MODELS[run.refname]
        resources = {}

        def add(logical_name, path):
            resources[logical_name] = Path(path).resolve()

        if stage == "prepare":
            add("reference_fasta", data_dir / f"{run.refname}.fasta")
            if run.target_input_contract is not None:
                add(
                    "target_input_contract",
                    self._resolve_input(run.target_input_contract),
                )
            model_candidate = self._resolve_input(Path(run.esm_model_path))
            if run.profile.is_paper_method and not model_candidate.is_file():
                raise FileNotFoundError(
                    "The paper profile requires an exact local ESM weight file via "
                    "--esm-model-path; model-name downloads are not reproducible"
                )
            if model_candidate.is_file():
                add("esm_model_weights", model_candidate)
        elif stage == "homology":
            add("reference_fasta", data_dir / f"{run.refname}.fasta")
            add("reference_annotation", data_dir / f"{run.refname}.xlsx")
            if run.mode == "plmsearch":
                add(
                    "plmsearch_model",
                    self.repo_root / "plmsearchtools" / "model" / "plmsearch.sav",
                )
            else:
                add("reference_model", reference_model)
                add(
                    "prottrans_reference_model",
                    self.repo_root
                    / "models"
                    / PROTTRANS_REFERENCE_MODELS[run.refname],
                )
                if run.prottrans_model_path is not None:
                    prottrans_root = self._resolve_input(run.prottrans_model_path)
                    validate_prott5_xl_uniref50_snapshot(prottrans_root)
                    prottrans_files = sorted(
                        path
                        for path in prottrans_root.rglob("*")
                        if path.is_file()
                    )
                    if not prottrans_files:
                        raise FileNotFoundError(
                            "ProtTrans model snapshot is empty: "
                            f"{prottrans_root}"
                        )
                    for path in prottrans_files:
                        relative = path.relative_to(prottrans_root)
                        add(f"prottrans_snapshot/{relative}", path)
                add(
                    "orthofinder_entrypoint",
                    self.repo_root / "tools" / "OrthoFinder" / "orthofinder.py",
                )
                add(
                    "usalign_binary",
                    self.repo_root / "tools" / "USalign" / "USalign",
                )
                add("structure_maplist", self._resolve_input(run.maplist))
                for index, path in enumerate(
                    sorted(self._resolve_input(run.structure).glob("*.pdb"))
                ):
                    add(f"target_structure_{index:06d}_{path.name}", path)
                reference_structure_dir = (
                    self.repo_root / "struct_data" / run.refname
                )
                reference_structures = sorted(
                    reference_structure_dir.glob("*.pdb")
                )
                if not reference_structures:
                    raise FileNotFoundError(
                        "Reference structure directory contains no PDB files: "
                        f"{reference_structure_dir}"
                    )
                for index, path in enumerate(reference_structures):
                    add(f"reference_structure_{index:06d}_{path.name}", path)
        elif stage == "draft_model":
            add("reference_model", reference_model)
        elif stage == "eggnog":
            add(
                "eggnog_mapper_entrypoint",
                self.repo_root / "tools" / "eggnog_mapper" / "emapper.py",
            )
            add("eggnog_database", eggnog_data / "eggnog.db")
            add(
                "eggnog_taxon_diamond_database",
                resolve_eggnog_diamond_db(
                    run.refname,
                    eggnog_data,
                    mode=run.profile.eggnog_database_mode,
                ),
            )
        elif stage == "clean" and run.clean:
            clean_root = self.repo_root / "tools" / "CLEAN" / "app"
            clean_entrypoint = clean_root / "CLEAN_infer_fasta.py"
            if self.verify_outputs and self._is_exact_paper_profile(run):
                self._assert_repo_local_path(
                    clean_root,
                    label="Formal CLEAN application root",
                )
                self._assert_repo_local_path(
                    clean_entrypoint,
                    label="Formal CLEAN entry point",
                )
            add(
                "clean_entrypoint",
                clean_entrypoint,
            )
            if self._is_exact_paper_profile(run):
                clean_esm_weights = self._resolve_input(
                    Path(run.esm_model_path)
                )
                add("clean_esm_model_weights", clean_esm_weights)
                add(
                    "clean_esm_regression_weights",
                    clean_esm_weights.with_name(
                        clean_esm_weights.stem
                        + "-contact-regression.pt"
                    ),
                )
        elif stage == "deepec" and run.deepec:
            tools_root = self.repo_root / "tools"
            deepec_root = (
                tools_root / "DeepProZyme"
            )
            for path in (tools_root, deepec_root):
                if path.is_symlink():
                    raise RuntimeError(
                        "Formal DeepEC bundle path does not permit symlinks: "
                        f"{path}"
                    )
            add(
                "deepec_entrypoint",
                deepec_root / "run_deepectransformer.py",
            )
            if deepec_root.is_dir():
                resolved_root = deepec_root.resolve()
                try:
                    resolved_root.relative_to(self.repo_root)
                except ValueError as error:
                    raise RuntimeError(
                        "DeepEC bundle root escapes the repository: "
                        f"{deepec_root} -> {resolved_root}"
                    ) from error
                for path in sorted(deepec_root.rglob("*")):
                    if path.is_symlink():
                        raise RuntimeError(
                            "Formal DeepEC bundle does not permit symlinks: "
                            f"{path}"
                        )
                    if (
                        not path.is_file()
                        or "__pycache__" in path.parts
                        or path.suffix == ".pyc"
                    ):
                        continue
                    resolved = path.resolve(strict=True)
                    try:
                        resolved.relative_to(resolved_root)
                    except ValueError as error:
                        raise RuntimeError(
                            "DeepEC bundle contains a file escaping its root: "
                            f"{path} -> {resolved}"
                        ) from error
                    relative = path.relative_to(deepec_root)
                    add(f"deepec_bundle/{relative}", path)
            if run.deepec_python is not None:
                add(
                    "deepec_python",
                    self._resolve_input(run.deepec_python),
                )
        elif stage == "rhea" and run.plm_annotation:
            add(
                "plmsearch_model",
                self.repo_root / "plmsearchtools" / "model" / "plmsearch.sav",
            )
            add(
                "rhea_uniprot_mapping",
                self.repo_root / "rhea" / "rhea2uniprot_sprot.tsv",
            )
            for index in range(1, 7):
                shard = (
                    self.repo_root
                    / "rhea"
                    / f"rhea_protembedding{index}.pkl"
                )
                add(f"rhea_embedding_{index}", shard)
        elif stage == "merge_annotations":
            add("rhea_rdf", self.repo_root / "rhea" / "rhea.rdf")
            add(
                "rhea_to_kegg_mapping",
                data_dir / "rhea2kegg_reaction.tsv",
            )
        elif stage == "reaction_pool":
            add("reference_model", reference_model)
            add(
                "general_to_specific_mapping",
                data_dir / "general2specific.pkl",
            )
            add(
                "reaction_pool",
                self.repo_root / "models" / "reaction_pool_final.json",
            )
        elif stage == "gapfill":
            add("reference_model", reference_model)
            add("reference_fasta", data_dir / f"{run.refname}.fasta")
            if run.growth_medium == "full":
                add(
                    "selected_full_medium",
                    data_dir / f"{run.refname}_full_medium.pkl",
                )
            if run.profile.gapfill_search_backend == "tblastn_genome":
                add(
                    "target_gapfill_genome",
                    self._resolve_input(run.gapfill_genome),
                )
                add(
                    "target_gapfill_coordinate_map",
                    self._resolve_input(run.gapfill_coordinate_map),
                )
                add(
                    "target_gapfill_coordinate_segments",
                    self._resolve_input(run.gapfill_coordinate_segments),
                )
                add(
                    "target_gapfill_coordinate_contract",
                    self._resolve_input(run.gapfill_coordinate_contract),
                )
            elif run.profile.gapfill_search_backend == "tblastn_cds":
                add(
                    "target_gapfill_nucleotide",
                    self._resolve_input(run.gapfill_cds),
                )
        return resources

    def _record_resources_for_stage(self, run, stage):
        records = self._resource_records_for_stage(run, stage)
        self.preflight_by_stage[stage] = records
        return records

    def _resource_records_for_stage(self, run, stage):
        resources = self._required_resources_for_stage(run, stage)
        missing = [
            f"{logical_name}: {path}"
            for logical_name, path in resources.items()
            if not path.is_file() or path.stat().st_size == 0
        ]
        if missing:
            raise FileNotFoundError(
                f"Formal AlphaGEM preflight for stage {stage!r} "
                "found missing resources:\n"
                + "\n".join(missing)
            )
        records = {
            logical_name: {
                "observed_path": str(path),
                **self._file_record(path),
            }
            for logical_name, path in resources.items()
        }
        if (
            stage == "eggnog"
            and self._is_exact_paper_profile(run)
            and "eggnog_database" in records
        ):
            records["eggnog_database"] = (
                self._paper_eggnog_database_record(
                    resources["eggnog_database"],
                    records["eggnog_database"],
                )
            )
            taxon_contract = PAPER_EGGNOG_TAXON_DATABASE_CONTRACTS.get(
                run.refname
            )
            if taxon_contract is not None:
                records["eggnog_taxon_diamond_database"] = (
                    self._paper_eggnog_taxon_database_record(
                        resources["eggnog_taxon_diamond_database"],
                        records["eggnog_taxon_diamond_database"],
                        taxon_contract,
                    )
                )
        if stage == "rhea" and self._is_exact_paper_profile(run):
            contracted_resources = PAPER_RHEA_ARTIFACT_CONTRACT["resources"]
            missing_contracts = set(contracted_resources) - set(records)
            if missing_contracts:
                raise RuntimeError(
                    "Formal Rhea resource contract is incomplete: "
                    + ", ".join(sorted(missing_contracts))
                )
            for logical_name, specification in contracted_resources.items():
                records[logical_name] = self._paper_rhea_resource_record(
                    logical_name,
                    resources[logical_name],
                    records[logical_name],
                    specification,
                )
        return records

    def _paper_eggnog_database_record(self, path, observed_record=None):
        expected_path = self._assert_repo_local_path(
            self.repo_root / PAPER_EGGNOG_DATABASE_CONTRACT["path"],
            label="Formal eggNOG database",
            require_file=True,
        ).resolve(strict=True)
        path = Path(path).resolve(strict=True)
        if path != expected_path:
            raise RuntimeError(
                "Formal eggNOG database is not the contracted repository "
                f"file: {path}"
            )
        observed_record = observed_record or {
            "observed_path": str(path),
            **self._file_record(path),
        }
        expected_identity = {
            key: PAPER_EGGNOG_DATABASE_CONTRACT[key]
            for key in ("bytes", "sha256")
        }
        observed_identity = {
            key: observed_record[key]
            for key in ("bytes", "sha256")
        }
        if observed_identity != expected_identity:
            raise RuntimeError(
                "Formal eggNOG database differs from the verified "
                "materialized database: "
                + json.dumps(
                    {
                        "expected": expected_identity,
                        "observed": observed_identity,
                    },
                    sort_keys=True,
                )
            )
        try:
            connection = sqlite3.connect(
                f"{path.as_uri()}?mode=ro",
                uri=True,
            )
            try:
                row = connection.execute(
                    "SELECT * FROM version LIMIT 1"
                ).fetchone()
            finally:
                connection.close()
        except sqlite3.Error as error:
            raise RuntimeError(
                "Formal eggNOG database is not a readable contracted SQLite "
                f"database: {path}"
            ) from error
        expected_version = PAPER_EGGNOG_DATABASE_CONTRACT[
            "expected_database_version"
        ]
        if row is None or str(row[0]) != expected_version:
            raise RuntimeError(
                "Formal eggNOG database has an unexpected internal version: "
                f"expected={expected_version!r}, observed={row!r}"
            )
        return {
            **observed_record,
            "database_version": str(row[0]),
            "identity_source": PAPER_EGGNOG_DATABASE_CONTRACT[
                "identity_source"
            ],
            "official_download_url": PAPER_EGGNOG_DATABASE_CONTRACT[
                "official_download_url"
            ],
            "download_provenance": PAPER_EGGNOG_DATABASE_CONTRACT[
                "download_provenance"
            ],
        }

    def _paper_eggnog_taxon_database_record(
        self,
        path,
        observed_record,
        specification,
    ):
        path = Path(path).resolve(strict=True)
        if path.name != specification["filename"]:
            raise RuntimeError(
                "Formal eggNOG taxonomic DIAMOND database has an unexpected "
                f"filename: expected={specification['filename']!r}, "
                f"observed={path.name!r}"
            )
        expected_identity = {
            key: specification[key] for key in ("bytes", "sha256")
        }
        observed_identity = {
            key: observed_record[key] for key in ("bytes", "sha256")
        }
        if observed_identity != expected_identity:
            raise RuntimeError(
                "Formal eggNOG taxonomic DIAMOND database differs from the "
                "verified Linux x86_64 DIAMOND v2.1.10 artifact derived from "
                "the verified eggNOG Fungi selected FASTA: "
                + json.dumps(
                    {
                        "expected": expected_identity,
                        "observed": observed_identity,
                    },
                    sort_keys=True,
                )
            )
        return {
            **observed_record,
            "database_role": specification["database_role"],
            "diamond_version": specification["diamond_version"],
            "diamond_format_version": specification["diamond_format_version"],
            "database_sequences": specification["sequences"],
            "database_letters": specification["letters"],
            "source_urls": specification["source_urls"],
            "source_sha256": specification["source_sha256"],
            "selection_contract": specification["selection_contract"],
            "identity_source": specification["identity_source"],
        }

    def _paper_rhea_resource_record(
        self,
        logical_name,
        path,
        observed_record,
        specification,
    ):
        expected_path = self._assert_repo_local_path(
            self.repo_root / specification["path"],
            label=f"Formal {logical_name} artifact",
            require_file=True,
        ).resolve(strict=True)
        path = Path(path).resolve(strict=True)
        if path != expected_path:
            raise RuntimeError(
                "Formal Rhea resource is not the contracted repository file: "
                f"{path}"
            )
        expected_identity = {
            key: specification[key]
            for key in ("bytes", "sha256")
        }
        observed_identity = {
            key: observed_record[key]
            for key in ("bytes", "sha256")
        }
        if observed_identity != expected_identity:
            raise RuntimeError(
                f"Formal {logical_name} differs from the verified author "
                "archive artifact: "
                + json.dumps(
                    {
                        "expected": expected_identity,
                        "observed": observed_identity,
                    },
                    sort_keys=True,
                )
            )
        record = {
            **observed_record,
            "identity_source": PAPER_RHEA_ARTIFACT_CONTRACT[
                "identity_source"
            ],
            "source_archive_sha256": PAPER_RHEA_ARTIFACT_CONTRACT[
                "source_archive_sha256"
            ],
            "paper_declared_version": PAPER_RHEA_ARTIFACT_CONTRACT[
                "paper_declared_version"
            ],
            "evidence_boundary": PAPER_RHEA_ARTIFACT_CONTRACT[
                "evidence_boundary"
            ],
            "archive_member_path": specification["archive_member_path"],
            "artifact_role": specification["artifact_role"],
        }
        for key in (
            "audited_key_count",
            "audited_data_row_count",
            "audited_unique_accession_count",
        ):
            if key in specification:
                record[key] = specification[key]
        if logical_name.startswith("rhea_embedding_"):
            for key in (
                "upstream_uniprot_snapshot",
                "embedding_model_identity",
                "pooling_semantics",
                "sharding_rule",
                "exact_regeneration_from_released_inputs",
                "embedding_set_audit",
            ):
                record[key] = PAPER_RHEA_ARTIFACT_CONTRACT[key]
        else:
            record["upstream_uniprot_snapshot"] = (
                PAPER_RHEA_ARTIFACT_CONTRACT[
                    "upstream_uniprot_snapshot"
                ]
            )
            record["exact_regeneration_from_released_inputs"] = (
                PAPER_RHEA_ARTIFACT_CONTRACT[
                    "exact_regeneration_from_released_inputs"
                ]
            )
        return record

    def _assert_resources_for_stage(self, run, stage):
        expected = self.preflight_by_stage.get(stage)
        if expected is None:
            raise RuntimeError(
                f"Formal stage {stage!r} lacks its resource preflight record"
            )
        observed = self._resource_records_for_stage(run, stage)
        if observed != expected:
            raise RuntimeError(
                "Formal stage resource changed after preflight: "
                f"{stage}"
            )
        return observed

    def _required_resource_bundles_for_stage(
        self,
        run: PipelineRun,
        stage: str,
    ):
        if not self._is_exact_paper_profile(run):
            return {}
        names = []
        if stage == "homology":
            names.append(
                "plmsearch"
                if run.mode == "plmsearch"
                else "orthofinder"
            )
            if run.mode == "structure_alignment":
                names.append("usalign")
        elif stage == "eggnog":
            names.extend(("eggnog_mapper", "eggnog_static_data"))
        elif stage == "clean" and run.clean:
            names.append("clean")
        elif stage == "deepec" and run.deepec:
            names.append("deepec_transformer")
        elif stage == "rhea" and run.plm_annotation:
            names.append("plmsearch")
        return {
            name: PAPER_RUNTIME_TOOL_BUNDLE_CONTRACTS[name]
            for name in names
        }

    def _record_paper_bundles_for_stage(self, run, stage):
        records = self._paper_bundle_records(
            self._required_resource_bundles_for_stage(run, stage)
        )
        self.resource_bundles_by_stage[stage] = records
        return records

    def _paper_bundle_records(self, specifications):
        return {
            logical_name: self._paper_bundle_record(
                logical_name,
                specification,
            )
            for logical_name, specification in specifications.items()
        }

    def _paper_bundle_record(self, logical_name, specification):
        root = self.repo_root / specification["path"]
        self._assert_repo_local_path(
            root,
            label=f"Formal {logical_name} bundle",
            require_file=root.is_file(),
        )
        if not root.exists() or not (root.is_file() or root.is_dir()):
            raise FileNotFoundError(
                f"Formal {logical_name} bundle is missing: {root}"
            )
        excluded = tuple(
            Path(relative).parts
            for relative in specification["excluded_relative_roots"]
        )

        def is_excluded(relative):
            return any(
                relative.parts[: len(prefix)] == prefix
                for prefix in excluded
            )

        candidates = [root] if root.is_file() else sorted(root.rglob("*"))
        entries = []
        for path in candidates:
            relative_to_bundle = (
                Path(path.name)
                if root.is_file()
                else path.relative_to(root)
            )
            if (
                "__pycache__" in relative_to_bundle.parts
                or relative_to_bundle.suffix == ".pyc"
            ):
                raise RuntimeError(
                    "Formal paper bundle contains forbidden Python "
                    f"bytecode: {path}"
                )
            if is_excluded(relative_to_bundle):
                continue
            if path.is_symlink():
                raise RuntimeError(
                    f"Formal {logical_name} bundle contains a symlink: {path}"
                )
            if path.is_dir():
                continue
            if not path.is_file():
                raise RuntimeError(
                    "Formal paper bundle contains a non-regular entry: "
                    f"{path}"
                )
            file_record = self._file_record(path)
            repo_relative = path.relative_to(self.repo_root).as_posix()
            entries.append(
                (repo_relative, file_record["sha256"], file_record["bytes"])
            )
        if not entries:
            raise RuntimeError(
                f"Formal {logical_name} bundle contains no regular files"
            )
        entries.sort()
        aggregate = hashlib.sha256()
        for relative, file_sha256, _size in entries:
            aggregate.update(
                f"{file_sha256}  {relative}\n".encode("utf-8")
            )
        observed_identity = {
            "file_count": len(entries),
            "bytes": sum(size for _relative, _sha256, size in entries),
            "sha256": aggregate.hexdigest(),
        }
        expected_identity = {
            key: specification[key]
            for key in ("file_count", "bytes", "sha256")
        }
        identity_scope = specification.get(
            "identity_scope",
            "author_archive_tree",
        )
        evidence_label = (
            "verified post-install runtime tree"
            if identity_scope == "post_install_runtime_tree"
            else "verified author resource archive"
        )
        if observed_identity != expected_identity:
            raise RuntimeError(
                f"Formal {logical_name} bundle differs from the "
                f"{evidence_label}: "
                + json.dumps(
                    {
                        "expected": expected_identity,
                        "observed": observed_identity,
                    },
                    sort_keys=True,
                )
            )
        record = {
            "observed_path": str(root),
            "paper_declared_version": specification[
                "paper_declared_version"
            ],
            "identity_scope": identity_scope,
            "identity_source": (
                "verified_post_install_runtime_tree"
                if identity_scope == "post_install_runtime_tree"
                else "verified_author_data_source_archive"
            ),
            "source_archive_sha256": PAPER_RESOURCE_ARCHIVE_SHA256,
            "excluded_relative_roots": list(
                specification["excluded_relative_roots"]
            ),
            "forbidden_generated_files": [
                "**/__pycache__/**",
                "**/*.pyc",
            ],
            **observed_identity,
        }
        if "source_archive_identity" in specification:
            record["source_archive_identity"] = dict(
                specification["source_archive_identity"]
            )
        if "tool_expected_database_version" in specification:
            record["tool_expected_database_version"] = specification[
                "tool_expected_database_version"
            ]
        return record

    def _required_tools_for_stage(self, run: PipelineRun, stage: str):
        if stage == "homology" and run.mode == "plmsearch":
            return {
                "blastp": (["blastp", "-version"], "2.16.0"),
                "diamond": (["diamond", "version"], "2.1.10"),
            }
        if stage == "homology" and run.mode == "structure_alignment":
            return {
                "blastp": (["blastp", "-version"], "2.16.0"),
                "foldseek": (["foldseek", "version"], "9.427df8a"),
                "orthofinder_declared": (
                    self._python_string_assignment_command(
                        self.repo_root
                        / "tools"
                        / "OrthoFinder"
                        / "scripts_of"
                        / "util.py",
                        "version",
                    ),
                    "3.0.1b1",
                ),
                "usalign": (
                    [
                        str(
                            self.repo_root
                            / "tools"
                            / "USalign"
                            / "USalign"
                        )
                    ],
                    "20241108",
                ),
            }
        if stage == "eggnog":
            return {
                "diamond": (["diamond", "version"], "2.1.10"),
                "eggnog_mapper_declared": (
                    self._python_string_assignment_command(
                        self.repo_root
                        / "tools"
                        / "eggnog_mapper"
                        / "eggnogmapper"
                        / "version.py",
                        "__VERSION__",
                    ),
                    "2.1.12",
                ),
            }
        if stage == "deepec" and run.deepec:
            if run.deepec_python is None:
                return {}
            python = str(self._resolve_input(run.deepec_python))
            return {
                "deepec_python": (
                    [python, "--version"],
                    "Python 3.7.7",
                ),
                "deepec_torch": (
                    [
                        python,
                        "-c",
                        "import torch; print(torch.__version__)",
                    ],
                    "1.7.0",
                ),
                "deepec_transformers": (
                    [
                        python,
                        "-c",
                        "import transformers; print(transformers.__version__)",
                    ],
                    "3.5.1",
                ),
            }
        if stage == "gapfill":
            tools = {
                "cobra": (
                    [
                        sys.executable,
                        "-c",
                        "import cobra; print(cobra.__version__)",
                    ],
                    "0.29.1",
                )
            }
            if run.profile.gapfill_search_backend in {
                "tblastn_genome",
                "tblastn_cds",
            }:
                tools["tblastn"] = (["tblastn", "-version"], "2.16.0")
                tools["makeblastdb"] = (
                    ["makeblastdb", "-version"],
                    "2.16.0",
                )
            return tools
        if stage in {"draft_model", "reaction_pool", "validate"}:
            return {
                "cobra": (
                    [
                        sys.executable,
                        "-c",
                        "import cobra; print(cobra.__version__)",
                    ],
                    "0.29.1",
                )
            }
        return {}

    @staticmethod
    def _python_string_assignment_command(path: Path, variable: str):
        script = (
            "import ast, pathlib, sys; "
            "tree = ast.parse(pathlib.Path(sys.argv[1]).read_text("
            "encoding='utf-8')); "
            "node = next(node for node in tree.body "
            "if isinstance(node, ast.Assign) "
            "and any(isinstance(target, ast.Name) "
            "and target.id == sys.argv[2] for target in node.targets)); "
            "value = ast.literal_eval(node.value); "
            "print(value if isinstance(value, str) else "
            "sys.exit('assignment is not a string'))"
        )
        return [sys.executable, "-c", script, str(path), variable]

    @staticmethod
    def _tool_version(command) -> str:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        return (result.stdout or result.stderr).strip().splitlines()[0]

    @staticmethod
    def _tool_version_matches(tool, expected, observed) -> bool:
        if expected is None:
            return True
        if tool in {
            "deepec_python",
            "deepec_torch",
            "deepec_transformers",
            "orthofinder_declared",
            "eggnog_mapper_declared",
        }:
            return observed.strip() == expected
        return re.search(
            rf"(?<![0-9A-Za-z.]){re.escape(expected)}(?![0-9A-Za-z.])",
            observed,
        ) is not None

    def _validate_resume(self, run: PipelineRun) -> dict:
        self._assert_formal_run_identity(run)
        manifest_path = (
            self.repo_root
            / "working"
            / run.name
            / "alphagem_run_manifest.json"
        )
        if not manifest_path.is_file():
            raise RuntimeError(
                f"Cannot resume without a run manifest: {manifest_path}"
            )
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        import_record = manifest.get("imported_from_checkpoint")
        if import_record is not None or run.from_checkpoint is not None:
            if import_record is None or run.from_checkpoint is None:
                raise RuntimeError(
                    "Checkpoint-imported manifests require --from-checkpoint, "
                    "and --from-checkpoint requires an imported manifest"
                )
            return self._validate_checkpoint_resume(run, manifest)
        expected = {
            "schema_version": self.MANIFEST_SCHEMA_VERSION,
            "name": run.name,
            "refname": run.refname,
            "mode": run.mode,
            "profile": asdict(run.profile),
            "solver": run.solver,
            "target_input_contract": self._target_input_contract_record(run),
            "blast_threads": run.blast_threads,
            "esm_token_budget": run.esm_token_budget,
            "clean": run.clean,
            "deepec": run.deepec,
            "deepec_python": (
                str(self._resolve_input(run.deepec_python))
                if run.deepec_python is not None
                else None
            ),
            "plm_annotation": run.plm_annotation,
            "growth_medium": run.growth_medium,
            "essential_reactions": sorted(run.essential_reactions or []),
            "reaction_pool_policies": self._reaction_pool_policy_record(run),
            "source_code": (
                self._run_source_identity
                if self.verify_outputs and self._is_exact_paper_profile(run)
                else self._source_code_identity()
            ),
            "runtime_environment": (
                self._run_runtime_environment
                if self.verify_outputs and self._is_exact_paper_profile(run)
                else self._runtime_environment()
            ),
            "paper_method_evidence": self._paper_method_evidence(run),
        }
        mismatches = {
            key: (value, manifest.get(key))
            for key, value in expected.items()
            if manifest.get(key) != value
        }
        recorded_solver = manifest.get("resolved_solver")
        if (
            recorded_solver is not None
            and self.solver_record is not None
            and recorded_solver != self.solver_record
        ):
            mismatches["resolved_solver"] = (
                self.solver_record,
                recorded_solver,
            )
        # Older manifests predate this field; the default value is compatible.
        if manifest.get("gapfill_tie_break", "strict") != run.gapfill_tie_break:
            mismatches["gapfill_tie_break"] = (
                run.gapfill_tie_break,
                manifest.get("gapfill_tie_break"),
            )
        if (
            manifest.get("gapfill_unmapped_hit_policy", "error")
            != run.gapfill_unmapped_hit_policy
        ):
            mismatches["gapfill_unmapped_hit_policy"] = (
                run.gapfill_unmapped_hit_policy,
                manifest.get("gapfill_unmapped_hit_policy"),
            )
        source_hash = self._file_record(self._resolve_input(run.fasta))["sha256"]
        if manifest.get("source_fasta", {}).get("sha256") != source_hash:
            mismatches["source_fasta.sha256"] = (
                source_hash,
                manifest.get("source_fasta", {}).get("sha256"),
            )
        current_gapfill_input = self._gapfill_sequence_record(run)
        manifest_gapfill_input = manifest.get("gapfill_sequence_input")
        if (
            "gapfill" in manifest.get("completed_stages", [])
            and not self._same_gapfill_input_identity(
                manifest_gapfill_input,
                current_gapfill_input,
            )
        ):
            mismatches["gapfill_sequence_input"] = (
                current_gapfill_input,
                manifest_gapfill_input,
            )
        start_index = self.STAGES.index(run.start_at)
        expected_completed = list(self.STAGES[:start_index])
        if manifest.get("completed_stages", []) != expected_completed:
            mismatches["completed_stages"] = (
                expected_completed,
                manifest.get("completed_stages", []),
            )
        expected_boundary = self._manifest_gapfill_method_boundary(
            run,
            expected_completed,
        )
        if manifest.get("gapfill_method_boundary") != expected_boundary:
            mismatches["gapfill_method_boundary"] = (
                expected_boundary,
                manifest.get("gapfill_method_boundary"),
            )
        if mismatches:
            raise RuntimeError(
                "Resume manifest is incompatible with this run: "
                + json.dumps(mismatches, sort_keys=True)
            )
        if (
            self.verify_outputs
            and manifest.get("verification_mode")
            != self._verification_mode(run, expected_completed)
        ):
            raise RuntimeError(
                "Cannot use an unverified manifest for a formal resume"
            )
        stage_records = manifest.get("stages", {})
        for stage in expected_completed:
            if stage not in stage_records:
                raise RuntimeError(
                    f"Resume manifest has no record for completed stage: {stage}"
                )
            stage_record = stage_records[stage]
            if self.verify_outputs:
                current_resources = self._resource_records_for_stage(
                    run,
                    stage,
                )
                recorded_resources = stage_record.get("resources", {})
                if set(recorded_resources) != set(current_resources):
                    raise RuntimeError(
                        "Resume resource contract is incomplete for stage "
                        f"{stage}: expected={sorted(current_resources)}, "
                        f"recorded={sorted(recorded_resources)}"
                    )
                for logical_name, current_record in current_resources.items():
                    if current_record != recorded_resources[logical_name]:
                        raise RuntimeError(
                            "Resume resource changed, is missing, or violates "
                            f"its formal contract: {stage}/{logical_name}"
                        )
                expected_tools = self._required_tools_for_stage(run, stage)
                recorded_tools = stage_record.get("tool_versions", {})
                if set(recorded_tools) != set(expected_tools):
                    raise RuntimeError(
                        "Resume tool contract is incomplete for stage "
                        f"{stage}: expected={sorted(expected_tools)}, "
                        f"recorded={sorted(recorded_tools)}"
                    )
                for tool, (_command, expected_version) in expected_tools.items():
                    if (
                        not self._tool_version_matches(
                            tool,
                            expected_version,
                            recorded_tools[tool],
                        )
                    ):
                        raise RuntimeError(
                            f"Resume tool version is incompatible for "
                            f"{stage}/{tool}: {recorded_tools[tool]}"
                        )
                expected_bundles = (
                    self._required_resource_bundles_for_stage(run, stage)
                )
                recorded_bundles = stage_record.get(
                    "resource_bundles",
                    {},
                )
                if set(recorded_bundles) != set(expected_bundles):
                    raise RuntimeError(
                        "Resume paper-tool bundle contract is incomplete "
                        f"for stage {stage}: "
                        f"expected={sorted(expected_bundles)}, "
                        f"recorded={sorted(recorded_bundles)}"
                    )
                for logical_name, specification in expected_bundles.items():
                    current_record = self._paper_bundle_record(
                        logical_name,
                        specification,
                    )
                    if current_record != recorded_bundles[logical_name]:
                        raise RuntimeError(
                            "Resume paper-tool bundle changed or is "
                            f"incompatible: {stage}/{logical_name}"
                        )
                expected_artifacts = {
                    str(path.relative_to(self.repo_root))
                    for path in self._expected_outputs(stage, run)
                }
                recorded_artifacts = stage_record.get("artifacts", {})
                if set(recorded_artifacts) != expected_artifacts:
                    raise RuntimeError(
                        "Resume artifact contract is incomplete for stage "
                        f"{stage}: expected={sorted(expected_artifacts)}, "
                        f"recorded={sorted(recorded_artifacts)}"
                    )
            else:
                expected_resources = {
                    logical_name: Path(record["observed_path"])
                    for logical_name, record in stage_record.get(
                        "resources", {}
                    ).items()
                }
                recorded_resources = stage_record.get("resources", {})
                recorded_artifacts = stage_record.get("artifacts", {})
                for logical_name, path in expected_resources.items():
                    record = recorded_resources[logical_name]
                    if not path.is_file() or not self._same_file_identity(
                        record,
                        self._file_record(path),
                    ):
                        raise RuntimeError(
                            "Resume resource changed or is missing: "
                            f"{stage}/{logical_name}: {path}"
                        )
            for relative_path, record in recorded_artifacts.items():
                path = self.repo_root / relative_path
                if not path.is_file() or not self._same_file_identity(
                    record,
                    self._file_record(path),
                ):
                    raise RuntimeError(
                        f"Resume artifact changed or is missing: {path}"
                    )
        return manifest

    def _validate_checkpoint_resume(self, run: PipelineRun, manifest: dict) -> dict:
        """Validate a checkpoint-imported manifest written by the importer.

        Imported prefix stages are trusted only through the importer's
        per-file SHA-256 records, which are re-verified against the copied
        files now. Their provenance stays imported_from_original_run; this
        path never treats them as executed by the current code.
        """
        import_record = manifest.get("imported_from_checkpoint")
        required_keys = {
            "importer_version",
            "source_run_dir",
            "source_manifest",
            "source_completed_stages",
            "solver_change",
            "code_comparison",
            "copied_files",
            "new_source_identity",
        }
        if not isinstance(import_record, dict) or not required_keys <= set(
            import_record
        ):
            raise RuntimeError(
                "Checkpoint import record is incomplete: missing "
                f"{sorted(required_keys - set(import_record or {}))}"
            )
        declared_source = str(Path(run.from_checkpoint).resolve())
        recorded_source = str(Path(import_record["source_run_dir"]).resolve())
        if declared_source != recorded_source:
            raise RuntimeError(
                "--from-checkpoint does not match the imported manifest's "
                f"source: {declared_source} != {recorded_source}"
            )
        solver_change = import_record["solver_change"]
        if (
            not isinstance(solver_change, dict)
            or solver_change.get("to") != run.solver
            or solver_change.get("declared_solver_backend_change") is not True
        ):
            raise RuntimeError(
                "Checkpoint resume requires the explicitly declared solver "
                "backend change recorded at import time"
            )
        expected = {
            "schema_version": self.MANIFEST_SCHEMA_VERSION,
            "name": run.name,
            "refname": run.refname,
            "mode": run.mode,
            "profile": asdict(run.profile),
            "solver": run.solver,
            "blast_threads": run.blast_threads,
            "esm_token_budget": run.esm_token_budget,
            "clean": run.clean,
            "deepec": run.deepec,
            "deepec_python": (
                str(self._resolve_input(run.deepec_python))
                if run.deepec_python is not None
                else None
            ),
            "plm_annotation": run.plm_annotation,
            "growth_medium": run.growth_medium,
            "essential_reactions": sorted(run.essential_reactions or []),
            "reaction_pool_policies": self._reaction_pool_policy_record(run),
            "paper_method_evidence": self._paper_method_evidence(run),
        }
        mismatches = {
            key: (value, manifest.get(key))
            for key, value in expected.items()
            if manifest.get(key) != value
        }

        def strip_paths(record):
            if isinstance(record, dict):
                return {
                    key: strip_paths(value)
                    for key, value in record.items()
                    if key != "observed_path"
                }
            return record

        content_fields = {
            "target_input_contract": self._target_input_contract_record(run),
            "gapfill_sequence_input": self._gapfill_sequence_record(run),
            "gapfill_coordinate_contract": (
                self._gapfill_coordinate_contract_record(run)
            ),
        }
        for key, current in content_fields.items():
            recorded = manifest.get(key)
            if strip_paths(recorded) != strip_paths(current):
                mismatches[key] = (strip_paths(current), strip_paths(recorded))
        # Older manifests predate this field; the default value is compatible.
        if manifest.get("gapfill_tie_break", "strict") != run.gapfill_tie_break:
            mismatches["gapfill_tie_break"] = (
                run.gapfill_tie_break,
                manifest.get("gapfill_tie_break"),
            )
        if (
            manifest.get("gapfill_unmapped_hit_policy", "error")
            != run.gapfill_unmapped_hit_policy
        ):
            mismatches["gapfill_unmapped_hit_policy"] = (
                run.gapfill_unmapped_hit_policy,
                manifest.get("gapfill_unmapped_hit_policy"),
            )
        source_hash = self._file_record(self._resolve_input(run.fasta))[
            "sha256"
        ]
        if manifest.get("source_fasta", {}).get("sha256") != source_hash:
            mismatches["source_fasta.sha256"] = (
                source_hash,
                manifest.get("source_fasta", {}).get("sha256"),
            )
        if (
            manifest.get("source_code")
            != import_record["new_source_identity"]
            or self._source_code_identity()
            != import_record["new_source_identity"]
        ):
            mismatches["source_code"] = (
                import_record["new_source_identity"],
                self._source_code_identity(),
            )
        if mismatches:
            raise RuntimeError(
                "Checkpoint-imported manifest is incompatible with this run: "
                + json.dumps(mismatches, sort_keys=True)
            )
        start_index = self.STAGES.index(run.start_at)
        expected_completed = list(self.STAGES[:start_index])
        imported_stages = import_record["source_completed_stages"]
        if (
            manifest.get("completed_stages", []) != expected_completed
            or not set(imported_stages) <= set(expected_completed)
        ):
            raise RuntimeError(
                "Checkpoint import does not carry exactly the completed "
                f"prefix stages: {expected_completed}"
            )
        if (
            self.verify_outputs
            and manifest.get("verification_mode")
            != self._verification_mode(run, expected_completed)
        ):
            raise RuntimeError(
                "Cannot use an unverified manifest for a formal resume"
            )
        if (
            manifest.get("gapfill_method_boundary")
            != self._manifest_gapfill_method_boundary(run, expected_completed)
        ):
            raise RuntimeError(
                "Checkpoint-imported manifest has an incompatible gapfill "
                "method boundary"
            )
        stage_records = manifest.get("stages", {})
        for stage in imported_stages:
            stage_record = stage_records.get(stage)
            if not isinstance(stage_record, dict) or not isinstance(
                stage_record.get("imported_from_original_run"), dict
            ):
                raise RuntimeError(
                    "Checkpoint stage lacks its imported_from_original_run "
                    f"provenance: {stage}"
                )
            expected_artifacts = {
                str(path.relative_to(self.repo_root))
                for path in self._expected_outputs(stage, run)
            }
            recorded_artifacts = stage_record.get("artifacts", {})
            if set(recorded_artifacts) != expected_artifacts:
                raise RuntimeError(
                    "Checkpoint stage artifact contract is incomplete for "
                    f"stage {stage!r}: expected={sorted(expected_artifacts)}, "
                    f"recorded={sorted(recorded_artifacts)}"
                )
        copied_files = import_record["copied_files"]
        for relative_path, record in sorted(copied_files.items()):
            path = self.repo_root / relative_path
            if not path.is_file() or not self._same_file_identity(
                record, self._file_record(path)
            ):
                raise RuntimeError(
                    "Checkpoint-imported file changed or is missing: "
                    f"{path}"
                )
        self._checkpoint_import_record = import_record
        return manifest

    def _validate_run(self, run: PipelineRun, stages: List[str]) -> None:
        validate_job_name(run.name)
        if run.gapfill_tie_break not in {"strict", "deterministic"}:
            raise ValueError(
                f"Unsupported gapfill tie_break: {run.gapfill_tie_break!r}"
            )
        if run.gapfill_unmapped_hit_policy not in {"error", "skip_with_audit"}:
            raise ValueError(
                "Unsupported gapfill unmapped_hit_policy: "
                f"{run.gapfill_unmapped_hit_policy!r}"
            )
        if run.refname not in REFERENCE_MODELS:
            choices = ", ".join(sorted(REFERENCE_MODELS))
            raise ValueError(f"Unsupported reference {run.refname!r}; choose: {choices}")
        if run.mode not in {"plmsearch", "structure_alignment"}:
            raise ValueError(f"Unsupported mode: {run.mode}")
        fasta = self._resolve_input(run.fasta)
        if not fasta.is_file():
            raise FileNotFoundError(f"Target FASTA not found: {fasta}")
        if (
            self.verify_outputs
            and run.profile.is_paper_method
            and not self._is_exact_paper_profile(run)
        ):
            raise ValueError(
                "A formal paper reproduction cannot use parameter overrides; "
                "run the derived profile with output verification disabled and "
                "label it as a sensitivity/control run"
            )
        if self.verify_outputs and self._is_exact_paper_profile(run):
            disabled_components = [
                name
                for name, enabled in (
                    ("CLEAN", run.clean),
                    ("DeepECtransformer", run.deepec),
                    ("PLMSearch-Rhea annotation", run.plm_annotation),
                )
                if not enabled
            ]
            if disabled_components:
                raise ValueError(
                    "A formal paper reproduction requires the complete "
                    "annotation ensemble; disabled components="
                    f"{disabled_components}"
                )
        if (
            self.verify_outputs
            and run.profile.is_paper_method
            and run.solver == "auto"
        ):
            raise ValueError(
                "Formal paper runs require an explicit --solver so the "
                "optimization and MILP backend can be versioned and audited"
            )
        target_contract = None
        if run.target_input_contract is not None:
            contract_path = self._resolve_input(run.target_input_contract)
            target_contract = validate_target_input_contract(
                contract_path,
                fasta,
            )
        if self.verify_outputs and run.profile.is_paper_method:
            if target_contract is None:
                raise ValueError(
                    "Formal paper runs require --target-input-contract to pin "
                    "the exact target proteome and provenance"
                )
        if self.verify_outputs and self._is_exact_paper_profile(run):
            self._initialize_formal_run_identity(run)
            source_code = self._run_source_identity
            if (
                not source_code["git_repository"]
                or not source_code["commit"]
            ):
                raise RuntimeError(
                    "Formal paper runs require a clean, committed Git "
                    "source tree; no repository commit could be established"
                )
            if source_code["dirty"]:
                raise RuntimeError(
                    "Formal paper runs require a clean, committed source tree; "
                    f"dirty entries={source_code['status_porcelain']}"
                )
        if run.profile.is_paper_method and run.gapfill_cds is not None:
            raise ValueError(
                "The paper method does not accept --gapfill-cds. Supply the "
                "raw genome and an audited coordinate map instead."
            )
        if (
            "gapfill" in stages
            and run.profile.gapfill_search_backend == "tblastn_genome"
        ):
            if (
                run.gapfill_genome is None
                or run.gapfill_coordinate_map is None
                or run.gapfill_coordinate_segments is None
                or run.gapfill_coordinate_contract is None
            ):
                raise ValueError(
                    "The paper gapfill stage requires --gapfill-genome, "
                    "--gapfill-coordinate-map, and "
                    "--gapfill-coordinate-segments, and "
                    "--gapfill-coordinate-contract"
                )
            genome = self._resolve_input(run.gapfill_genome)
            coordinate_map = self._resolve_input(run.gapfill_coordinate_map)
            coordinate_segments = self._resolve_input(
                run.gapfill_coordinate_segments
            )
            coordinate_contract = self._resolve_input(
                run.gapfill_coordinate_contract
            )
            if not genome.is_file():
                raise FileNotFoundError(
                    f"Gap-fill raw genome FASTA not found: {genome}"
                )
            if not coordinate_map.is_file():
                raise FileNotFoundError(
                    f"Gap-fill coordinate map not found: {coordinate_map}"
                )
            if not coordinate_segments.is_file():
                raise FileNotFoundError(
                    "Gap-fill coordinate-segment map not found: "
                    f"{coordinate_segments}"
                )
            if not coordinate_contract.is_file():
                raise FileNotFoundError(
                    "Gap-fill coordinate contract not found: "
                    f"{coordinate_contract}"
                )
            if target_contract is None:
                raise RuntimeError(
                    "Formal gap-fill validation has no target input contract"
                )
            validated_coordinate_contract = validate_gapfill_coordinate_contract(
                coordinate_contract,
                fasta_path=fasta,
                genome_path=genome,
                coordinate_map_path=coordinate_map,
                coordinate_segments_path=coordinate_segments,
                expected_case_id=target_contract["case_id"],
            )
            validate_gapfill_genome(
                coordinate_map=coordinate_map,
                genome_fasta=genome,
                protein_fasta=fasta,
                coordinate_segments=coordinate_segments,
                translation_exceptions=validated_coordinate_contract[
                    "translation_exceptions"
                ],
            )
        elif (
            "gapfill" in stages
            and run.profile.gapfill_search_backend == "tblastn_cds"
        ):
            if run.gapfill_cds is None:
                raise ValueError(
                    "The paper gapfill stage requires --gapfill-cds; "
                    "a raw contig assembly is insufficient without GFF mapping"
                )
            gapfill_cds = self._resolve_input(run.gapfill_cds)
            if not gapfill_cds.is_file():
                raise FileNotFoundError(
                    f"Gap-fill nucleotide CDS FASTA not found: {gapfill_cds}"
                )
            validate_gapfill_cds(gapfill_cds, fasta)
        if (
            run.mode == "structure_alignment"
            and "homology" in stages
        ):
            if run.maplist is None or run.structure is None:
                raise ValueError(
                    "structure_alignment requires both --list/--maplist and --structure"
                )
            if not self._resolve_input(run.maplist).is_file():
                raise FileNotFoundError(f"Mapping list not found: {run.maplist}")
            if not self._resolve_input(run.structure).is_dir():
                raise FileNotFoundError(f"Structure directory not found: {run.structure}")
            if not any(self._resolve_input(run.structure).glob("*.pdb")):
                raise ValueError(
                    f"Structure directory contains no PDB files: {run.structure}"
                )
            if run.profile.is_paper_method:
                if run.prottrans_model_path is None:
                    raise ValueError(
                        "The paper structure-homology stage requires "
                        "--prottrans-model-path"
                    )
                prottrans_model = self._resolve_input(run.prottrans_model_path)
                if not prottrans_model.is_dir():
                    raise NotADirectoryError(
                        f"ProtTrans model snapshot not found: {prottrans_model}"
                    )
        if (
            self.verify_outputs
            and self._is_exact_paper_profile(run)
            and run.deepec
        ):
            if run.deepec_python is None and "deepec" in stages:
                raise ValueError(
                    "A formal paper DeepECtransformer stage requires "
                    "--deepec-python for the pinned Python 3.7 environment"
                )
            if run.deepec_python is not None:
                deepec_python = self._resolve_input(run.deepec_python)
                if not deepec_python.is_file() or not os.access(
                    deepec_python,
                    os.X_OK,
                ):
                    raise FileNotFoundError(
                        "DeepECtransformer Python is missing or not "
                        f"executable: {deepec_python}"
                    )
                run.deepec_python = deepec_python.resolve(strict=True)

    def _stage_prepare(self, run: PipelineRun) -> None:
        working = self.repo_root / "working" / run.name
        working.mkdir(parents=True, exist_ok=True)
        self._module("fasta_handle").handle(
            str(self._resolve_input(run.fasta)),
            run.name,
            strict=(
                self.verify_outputs
                and self._is_exact_paper_profile(run)
            ),
        )
        target_embedding = working / f"{run.name}_embedding.pkl"
        target_metadata = working / f"{run.name}_embedding.meta.json"
        model_candidate = self._resolve_input(Path(run.esm_model_path))
        model_path = model_candidate if model_candidate.is_file() else None
        pooling = (
            "mean_residue_excluding_bos_eos"
            if run.profile.embedding_semantics == "corrected"
            else "legacy_truncated_slice_includes_eos_for_long_sequences"
        )
        target_fasta = working / f"{run.name}.fasta"
        if not self._embedding_cache_is_reusable(
            target_embedding,
            target_fasta,
            target_metadata,
            model_identifier=run.esm_model_path,
            model_path=model_path,
            pooling=pooling,
        ):
            self._module("plmsearch_embedding").embedding_generate(
                run.name,
                esm_model_path=run.esm_model_path,
                nogpu=False,
                token_budget=run.esm_token_budget,
                checkpoint_every=25,
                embedding_semantics=run.profile.embedding_semantics,
            )
        if target_embedding.is_file():
            assert_embedding_provenance(
                target_embedding,
                target_fasta,
                target_metadata,
                model_identifier=run.esm_model_path,
                model_path=model_path,
                pooling=pooling,
            )
        reference_fasta = (
            self.repo_root / "data_available" / f"{run.refname}.fasta"
        )
        reference_embedding = (
            self.repo_root / "working" / f"{run.refname}_embedding.pkl"
        )
        reference_metadata = (
            self.repo_root / "working" / f"{run.refname}_embedding.meta.json"
        )
        if not reference_fasta.is_file():
            raise FileNotFoundError(
                f"Reference FASTA not found: {reference_fasta}"
            )
        if run.mode == "plmsearch":
            if not self._embedding_cache_is_reusable(
                reference_embedding,
                reference_fasta,
                reference_metadata,
                model_identifier=run.esm_model_path,
                model_path=model_path,
                pooling=pooling,
            ):
                self._module("plmsearch_embedding").embedding_generate_base(
                    run.refname,
                    esm_model_path=run.esm_model_path,
                    nogpu=False,
                    token_budget=run.esm_token_budget,
                    checkpoint_every=25,
                    embedding_semantics=run.profile.embedding_semantics,
                )
            if reference_embedding.is_file():
                assert_embedding_provenance(
                    reference_embedding,
                    reference_fasta,
                    reference_metadata,
                    model_identifier=run.esm_model_path,
                    model_path=model_path,
                    pooling=pooling,
                )

    @staticmethod
    def _embedding_cache_is_reusable(
        embedding_path: Path,
        fasta_path: Path,
        metadata_path: Path,
        *,
        model_identifier: str,
        model_path: Optional[Path],
        pooling: str,
    ) -> bool:
        if not embedding_path.is_file():
            return False
        try:
            assert_embedding_provenance(
                embedding_path,
                fasta_path,
                metadata_path,
                model_identifier=model_identifier,
                model_path=model_path,
                pooling=pooling,
            )
        except (EOFError, OSError, ValueError, pickle.UnpicklingError) as error:
            print(
                "Embedding cache is not reusable and will be regenerated: "
                f"{embedding_path}: {error}"
            )
            return False
        return True

    def _stage_homology(self, run: PipelineRun) -> None:
        if run.mode == "plmsearch":
            self._run_plmsearch_homology(run)
        else:
            self._run_structure_homology(run)

    def _run_plmsearch_homology(self, run: PipelineRun) -> None:
        threshold = run.profile.plmsearch
        self._module("generate_list").generare_list_without_structure(run.name)
        self._module("plmsearch").ss_predictor(
            run.refname,
            run.name,
            threshold.ss_predictor,
        )
        self._module("ss_predict_choose").ss_predict_choose(
            run.name,
            run.refname,
            threshold.ss_predictor,
        )
        cov_filter = self._module("cov_filter")
        cov_filter.result_blast(run.name, run.refname, threshold.ss_predictor)
        cov_filter.cov_filter(
            run.name,
            run.refname,
            threshold.ss_predictor,
            threshold.coverage,
            threshold.identity,
        )
        cov_filter.bbh(
            run.name,
            run.refname,
            threshold.ss_predictor,
            threshold.bbh_coverage,
            threshold.bbh_identity,
            threshold.direction,
            threshold.bbh_backend,
            run.blast_threads,
        )
        cov_filter.merge(
            run.name,
            run.refname,
            threshold.ss_predictor,
            threshold.coverage,
            threshold.identity,
            threshold.bbh_coverage,
            threshold.bbh_identity,
            threshold.direction,
        )
        self._module("mkdir").mvfile(
            run.name,
            run.refname,
            threshold.ss_predictor,
            threshold.coverage,
            threshold.identity,
            threshold.bbh_coverage,
            threshold.bbh_identity,
            threshold.direction,
        )

    def _run_structure_homology(self, run: PipelineRun) -> None:
        threshold = run.profile.structure
        structure_dir = str(self._resolve_input(run.structure))
        maplist = str(self._resolve_input(run.maplist))
        refmodel = REFERENCE_MODELS[run.refname]
        formal_paper_run = (
            self.verify_outputs and self._is_exact_paper_profile(run)
        )

        staged_structure = self._module(
            "generate_list"
        ).generare_list_with_structure(
            run.name,
            maplist,
            structure_dir,
            allow_missing_structures=formal_paper_run,
        )
        staged_path = Path(staged_structure) if staged_structure else None
        if (
            staged_path is None
            or not staged_path.is_dir()
            or not any(
                path.is_file() and path.stat().st_size > 0
                for path in staged_path.glob("*.pdb")
            )
        ):
            raise RuntimeError(
                "Structure staging did not return a directory with non-empty "
                f"PDB files: {staged_structure!r}"
            )
        self._module("orthofinder_datahandle").datahandel(run.name, run.refname)
        self._module("US_align_knock1").US_align_find(
            run.name,
            staged_structure,
            run.refname,
            "",
            write_audit=formal_paper_run,
        )
        self._module("US_align_choose").US_align_choose(
            run.name,
            threshold.usalign,
            threshold.coverage,
            threshold.plddt,
        )
        self._module("foldseek_knock1").foldseekfind(
            staged_structure,
            run.name,
            run.refname,
            allow_score_fallback=not (
                self.verify_outputs
                and self._is_exact_paper_profile(run)
            ),
            require_fresh_output=(
                self.verify_outputs
                and self._is_exact_paper_profile(run)
            ),
        )
        self._module("foldseek_choose").foldseek_choose(
            run.name,
            run.refname,
            staged_structure,
        )
        self._module("foldseek_choose2").foldseek_choose2(
            run.refname,
            refmodel,
            run.name,
            threshold.transporter_tm_score,
            threshold.transporter_coverage,
            threshold.tm_score,
            threshold.coverage,
            threshold.plddt,
        )
        self._module("foldseekcluster").cluster(
            run.name,
            threshold.dbscan_eps,
            uptm=threshold.high_confidence_tm_score,
            upcov=threshold.high_confidence_coverage,
        )
        self._module("homolog_concat2").homo(run.name)
        self._module("filter_by_prottrans").filter_by_prottrans(
            run.name,
            run.refname,
            run.profile.prottrans_distance,
            model_path=(
                str(
                    self._resolve_input(
                        run.prottrans_model_path
                    ).resolve(strict=True)
                )
                if run.prottrans_model_path is not None
                else None
            ),
            local_files_only=run.profile.is_paper_method,
            write_audit=formal_paper_run,
        )

    def _stage_draft_model(self, run: PipelineRun) -> None:
        refmodel = REFERENCE_MODELS[run.refname]
        self._module("model_build1").modelbuild(refmodel, run.name)
        self._module("model_build2").modelbuild(refmodel, run.name)

    def _stage_prepare_nonhomologs(self, run: PipelineRun) -> None:
        self._module("use_clean").get_fasta(
            run.name,
            require_safe_paths=(
                self.verify_outputs
                and self._is_exact_paper_profile(run)
            ),
        )

    def _stage_eggnog(self, run: PipelineRun) -> None:
        self._module("use_eggnog").eggnog(
            run.name,
            run.refname,
            1,
            run.profile.eggnog_database_mode,
        )

    def _release_cuda_cache(self) -> None:
        """Release completed in-process GPU work before a child loads ESM."""

        gc.collect()
        torch = self._module("torch")
        cuda = getattr(torch, "cuda", None)
        is_available = getattr(cuda, "is_available", None)
        empty_cache = getattr(cuda, "empty_cache", None)
        if (
            callable(is_available)
            and callable(empty_cache)
            and is_available()
        ):
            empty_cache()

    def _stage_clean(self, run: PipelineRun) -> None:
        if run.clean:
            self._release_cuda_cache()
            self._module("use_clean").clean(
                run.name,
                require_fresh_output=(
                    self.verify_outputs
                    and self._is_exact_paper_profile(run)
                ),
                esm_model_path=(
                    str(self._resolve_input(Path(run.esm_model_path)))
                    if self.verify_outputs
                    and self._is_exact_paper_profile(run)
                    else None
                ),
                esm_token_budget=run.esm_token_budget,
            )

    def _stage_deepec(self, run: PipelineRun) -> None:
        if run.deepec:
            runtime_contract = (
                self._assert_deepec_runtime_contract(run)
            )
            self._module(
                "use_deepectransformer"
            ).use_deepectransformer(
                run.name,
                python_executable=(
                    str(self._resolve_input(run.deepec_python))
                    if run.deepec_python is not None
                    else None
                ),
                expected_python_sha256=(
                    runtime_contract["python_sha256"]
                ),
                expected_entrypoint_sha256=(
                    runtime_contract["entrypoint_sha256"]
                ),
            )

    def _assert_deepec_runtime_contract(
        self,
        run: PipelineRun,
    ):
        if not (
            self.verify_outputs
            and self._is_exact_paper_profile(run)
        ):
            return {
                "python_sha256": None,
                "entrypoint_sha256": None,
            }
        recorded_resources = self.preflight_by_stage.get("deepec", {})
        recorded_tools = self.tool_versions_by_stage.get("deepec", {})
        python_record = recorded_resources.get("deepec_python")
        if python_record is None or run.deepec_python is None:
            raise RuntimeError(
                "Formal DeepEC stage lacks its preflight interpreter record"
            )
        expected_resources = self._required_resources_for_stage(
            run,
            "deepec",
        )
        if set(recorded_resources) != set(expected_resources):
            raise RuntimeError(
                "DeepEC resource set changed after preflight"
            )
        for logical_name, path in expected_resources.items():
            if not path.is_file() or not self._same_file_identity(
                recorded_resources[logical_name],
                self._file_record(path),
            ):
                raise RuntimeError(
                    "DeepEC resource changed after preflight: "
                    f"{logical_name}={path}"
                )
        tool_specs = self._required_tools_for_stage(run, "deepec")
        if set(recorded_tools) != set(tool_specs):
            raise RuntimeError(
                "Formal DeepEC stage lacks its complete version preflight"
            )
        for tool, (command, expected) in tool_specs.items():
            observed = self._tool_version(command)
            if (
                not self._tool_version_matches(tool, expected, observed)
                or observed != recorded_tools[tool]
            ):
                raise RuntimeError(
                    "DeepEC runtime changed after preflight: "
                    f"{tool} expected={expected!r} "
                    f"preflight={recorded_tools[tool]!r} "
                    f"observed={observed!r}"
                )
        return {
            "python_sha256": python_record["sha256"],
            "entrypoint_sha256": recorded_resources[
                "deepec_entrypoint"
            ]["sha256"],
        }

    def _stage_rhea(self, run: PipelineRun) -> None:
        if run.plm_annotation:
            self._module("use_rhea").rhea(
                run.name,
                run.refname,
                run.profile.rhea_ss_predictor,
            )

    def _stage_merge_annotations(self, run: PipelineRun) -> None:
        self._module("multiple_annotition_for_nonhomogene").nonhome(
            run.name,
            run.clean,
            run.deepec,
            run.plm_annotation,
        )

    def _stage_reaction_pool(self, run: PipelineRun) -> None:
        self._module("general_reactions_to_spe").get_spe_reactions(
            run.name,
            run.profile.ensemble_score,
        )
        self._module("add_reactions_based_pool").model_reaction(
            run.name,
            REFERENCE_MODELS[run.refname],
            strict_identity=(
                self.verify_outputs and self._is_exact_paper_profile(run)
            ),
            already_present_policy=run.reaction_pool_already_present_policy,
            conflict_policy=run.reaction_pool_conflict_policy,
        )

    def _stage_gapfill(self, run: PipelineRun) -> None:
        self._module("gapfilling").gapfill(
            run.name,
            run.refname,
            REFERENCE_MODELS[run.refname],
            grothmedium=run.growth_medium,
            esstinal=run.essential_reactions,
            strategy=run.profile.gapfill_strategy,
            search_backend=run.profile.gapfill_search_backend,
            gapfill_cds=(
                str(self._resolve_input(run.gapfill_cds))
                if run.gapfill_cds is not None
                else None
            ),
            gapfill_genome=(
                str(self._resolve_input(run.gapfill_genome))
                if run.gapfill_genome is not None
                else None
            ),
            gapfill_coordinate_map=(
                str(self._resolve_input(run.gapfill_coordinate_map))
                if run.gapfill_coordinate_map is not None
                else None
            ),
            gapfill_coordinate_segments=(
                str(self._resolve_input(run.gapfill_coordinate_segments))
                if run.gapfill_coordinate_segments is not None
                else None
            ),
            gapfill_coordinate_contract=(
                str(self._resolve_input(run.gapfill_coordinate_contract))
                if run.gapfill_coordinate_contract is not None
                else None
            ),
            bitscore_threshold=run.profile.gapfill_bitscore_threshold,
            bitscore_operator=run.profile.gapfill_bitscore_operator,
            tie_break=run.gapfill_tie_break,
            unmapped_hit_policy=run.gapfill_unmapped_hit_policy,
            require_solver_tolerances=(
                self.verify_outputs and self._is_exact_paper_profile(run)
            ),
            require_optimal_solutions=(
                self.verify_outputs and self._is_exact_paper_profile(run)
            ),
            require_essential_reactions=(
                self.verify_outputs and self._is_exact_paper_profile(run)
            ),
        )

    def _stage_validate(self, run: PipelineRun) -> None:
        working = self.repo_root / "working" / run.name
        xml_path = working / f"{run.name}-GEM.xml"
        json_path = working / f"{run.name}-GEM.json"
        if not xml_path.is_file() or not json_path.is_file():
            raise RuntimeError(
                f"Pipeline ended without both final GEM formats: {working}"
            )
        if not self.verify_outputs:
            return
        cobra = self._module("cobra")
        xml_model = cobra.io.read_sbml_model(str(xml_path))
        json_model = cobra.io.load_json_model(str(json_path))
        differences = self._compare_models(xml_model, json_model)
        xml_solution = xml_model.optimize()
        json_solution = json_model.optimize()
        for label, solution in (
            ("xml", xml_solution),
            ("json", json_solution),
        ):
            if (
                solution.status != "optimal"
                or solution.objective_value is None
                or solution.objective_value <= 1e-5
            ):
                differences.append(
                    f"{label} growth status={solution.status}, "
                    f"objective={solution.objective_value}"
                )
        if differences:
            raise RuntimeError(
                "Final model validation failed: " + "; ".join(differences[:20])
            )
        payload = {
            "xml_sha256": self._file_record(xml_path)["sha256"],
            "json_sha256": self._file_record(json_path)["sha256"],
            "reactions": len(xml_model.reactions),
            "metabolites": len(xml_model.metabolites),
            "genes": len(xml_model.genes),
            "solver": str(xml_model.solver.interface),
            "status": xml_solution.status,
            "objective_value": float(xml_solution.objective_value),
            "formats_equivalent": True,
        }
        (working / "model_validation.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True),
            encoding="utf-8",
        )

    @staticmethod
    def _compare_models(left, right) -> List[str]:
        differences = []
        collections = (
            ("reaction", left.reactions, right.reactions),
            ("metabolite", left.metabolites, right.metabolites),
            ("gene", left.genes, right.genes),
        )
        for label, left_items, right_items in collections:
            left_ids = {item.id for item in left_items}
            right_ids = {item.id for item in right_items}
            if left_ids != right_ids:
                differences.append(
                    f"{label} IDs differ: "
                    f"left_only={len(left_ids - right_ids)}, "
                    f"right_only={len(right_ids - left_ids)}"
                )
        shared_reactions = {
            reaction.id for reaction in left.reactions
        } & {reaction.id for reaction in right.reactions}
        for reaction_id in sorted(shared_reactions):
            left_reaction = left.reactions.get_by_id(reaction_id)
            right_reaction = right.reactions.get_by_id(reaction_id)
            left_stoich = {
                metabolite.id: float(value)
                for metabolite, value in left_reaction.metabolites.items()
            }
            right_stoich = {
                metabolite.id: float(value)
                for metabolite, value in right_reaction.metabolites.items()
            }
            if left_stoich != right_stoich:
                differences.append(f"{reaction_id} stoichiometry differs")
            if (
                abs(left_reaction.lower_bound - right_reaction.lower_bound) > 1e-9
                or abs(left_reaction.upper_bound - right_reaction.upper_bound) > 1e-9
            ):
                differences.append(f"{reaction_id} bounds differ")
            if (
                left_reaction.gene_reaction_rule
                != right_reaction.gene_reaction_rule
            ):
                differences.append(f"{reaction_id} GPR differs")
        return differences

    def _highs_interface(self):
        try:
            return self._module("optlang_highs_interface")
        except ModuleNotFoundError:
            return self._module("src.optlang_highs_interface")

    def _apply_solver_selection(self, cobra, solver: str) -> None:
        """Apply the run solver, registering the bundled HiGHS interface.

        Registration inserts the bundled module into cobra's shared solver
        lookup (cobra.util.solver.solvers, the same object behind
        cobra.core.configuration.SOLVERS), so `--solver highs` flows through
        Configuration and model creation exactly like the shipped backends.
        """
        if solver == "highs":
            self._highs_interface().register_highs_interface()
        cobra.Configuration().solver = solver

    def _configure_solver(self, solver: str) -> None:
        if solver == "auto":
            return
        cobra = self._module("cobra")
        self._apply_solver_selection(cobra, solver)

    def _module(self, name: str):
        return self.module_loader(name)

    def _resolve_input(self, path: Path) -> Path:
        path = Path(path)
        return path if path.is_absolute() else self.repo_root / path

    def _assert_repo_local_path(
        self,
        path: Path,
        *,
        label: str,
        require_file: bool = False,
    ) -> Path:
        path = Path(path)
        if not path.is_absolute():
            path = self.repo_root / path
        try:
            relative = path.relative_to(self.repo_root)
        except ValueError as error:
            raise RuntimeError(
                f"{label} escapes the repository: {path}"
            ) from error
        current = self.repo_root
        for component in relative.parts:
            current = current / component
            if current.is_symlink():
                raise RuntimeError(
                    f"{label} path does not permit symlinks: {current}"
                )
        resolved = path.resolve(strict=require_file)
        try:
            resolved.relative_to(self.repo_root)
        except ValueError as error:
            raise RuntimeError(
                f"{label} resolves outside the repository: "
                f"{path} -> {resolved}"
            ) from error
        if require_file and (
            not path.is_file()
            or path.is_symlink()
            or path.stat().st_size == 0
        ):
            raise RuntimeError(
                f"{label} is missing, empty, or not a regular file: {path}"
            )
        return path

    def _expected_outputs(self, stage: str, run: PipelineRun) -> List[Path]:
        working = self.repo_root / "working" / run.name
        outputs = {
            "prepare": [
                working / f"{run.name}.fasta",
                working / f"{run.name}_embedding.pkl",
                working / f"{run.name}_embedding.meta.json",
            ]
            + (
                [
                    self.repo_root
                    / "working"
                    / f"{run.refname}_embedding.pkl",
                    self.repo_root
                    / "working"
                    / f"{run.refname}_embedding.meta.json",
                ]
                if run.mode == "plmsearch"
                else []
            ),
            "homology": [
                working / f"{run.name}.xlsx",
                working / f"matrix_homolog{run.name}.xlsx",
            ]
            + (
                [
                    working / f"aln{run.name}.csv",
                    working / f"matrix_foldseek_{run.name}.xlsx",
                ]
                + (
                    [
                        working / f"matrix_orthofinder{run.name}.xlsx",
                        working
                        / f"matrix_USalign_filtered{run.name}.xlsx",
                        working / f"matrix_homo_part1{run.name}.xlsx",
                        working
                        / f"matrix_foldseek_filtered1_{run.name}.xlsx",
                        working
                        / f"matrix_foldseek_filtered2_{run.name}.xlsx",
                        working / f"matrix_homo_part2{run.name}.xlsx",
                        working
                        / (
                            f"matrix_homolog{run.name}"
                            "_preforprottrans.xlsx"
                        ),
                        working / f"{run.name}.emb",
                        working / "structure_staging_audit.json",
                        working / "usalign_execution.json",
                        working / "prottrans_execution.json",
                    ]
                    if self.verify_outputs
                    and self._is_exact_paper_profile(run)
                    else []
                )
                + (
                    [working / "foldseek_execution.json"]
                    if self.verify_outputs
                    and self._is_exact_paper_profile(run)
                    else []
                )
                if run.mode == "structure_alignment"
                else []
            ),
            "draft_model": [
                working / f"model1{run.name}.yml",
                working / f"model2{run.name}.yml",
                working / "draft_model_audit.json",
            ],
            "prepare_nonhomologs": [working / f"{run.name}_homoleft.fasta"],
            "eggnog": [working / f"eggec2{run.name}.xlsx"],
            "clean": (
                [working / f"{run.name}_homoleft_maxsep.csv"]
                + (
                    [working / "clean_execution.json"]
                    if self._is_exact_paper_profile(run)
                    and self.verify_outputs
                    else []
                )
                if run.clean
                else []
            ),
            "deepec": (
                [
                    working
                    / f"{run.name}_deepec_result"
                    / "DeepECv2_result.txt"
                ]
                if run.deepec
                else []
            ),
            "rhea": (
                [working / f"{run.name}_rheaid_plmsearch.xlsx"]
                if run.plm_annotation
                else []
            ),
            "merge_annotations": [
                working / f"{run.name}multiple_score.xlsx",
                working / "annotation_merge_audit.json",
            ],
            "reaction_pool": [
                working / f"{run.name}gpr_score.xlsx",
                working / f"{run.name}-GEM_withgaps.yml",
                working / "reaction_pool_audit.json",
                working / "reaction_pool_coverage.json",
            ] + (
                [working / "reaction_pool_validation_v3.json"]
                if self.verify_outputs and self._is_exact_paper_profile(run)
                else []
            ),
            "gapfill": [
                working / f"{run.name}-GEM.xml",
                working / f"{run.name}-GEM.json",
                working / "gapfill_audit.json",
            ],
            "validate": [
                working / f"{run.name}-GEM.xml",
                working / f"{run.name}-GEM.json",
                working / "model_validation.json",
            ],
        }
        return outputs[stage]

    def _verify_stage_outputs(self, stage: str, run: PipelineRun):
        if stage in {
            "draft_model",
            "merge_annotations",
            "reaction_pool",
        }:
            # The formal reaction-pool check writes a versioned receipt, which
            # must exist before artifact hashes enter the run manifest.
            self._validate_stage_audit(stage, run)
        records = {}
        for path in self._expected_outputs(stage, run):
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError(
                    f"Stage {stage!r} did not produce a non-empty artifact: {path}"
                )
            records[str(path.relative_to(self.repo_root))] = self._file_record(path)
        if stage == "clean" and self._is_exact_paper_profile(run):
            self._validate_clean_audit(run)
        if (
            stage == "homology"
            and run.mode == "structure_alignment"
            and self._is_exact_paper_profile(run)
        ):
            self._validate_structure_staging_audit(run)
            self._validate_usalign_audit(run)
            self._validate_prottrans_audit(run)
            self._validate_foldseek_audit(run)
        if stage == "gapfill":
            working = self.repo_root / "working" / run.name
            audit_path = working / "gapfill_audit.json"
            audit = json.loads(audit_path.read_text(encoding="utf-8"))
            if self._is_exact_paper_profile(run):
                self._gapfill_method_boundary_record(audit)
                model_outputs = audit.get("model_outputs")
                output_paths = {
                    "xml": working / f"{run.name}-GEM.xml",
                    "json": working / f"{run.name}-GEM.json",
                }
                if not isinstance(model_outputs, dict):
                    raise RuntimeError(
                        "Formal gapfill audit lacks final model identities"
                    )
                for output_format, output_path in output_paths.items():
                    if not self._same_file_identity(
                        model_outputs.get(output_format),
                        self._file_record(output_path),
                    ):
                        raise RuntimeError(
                            "Formal gapfill audit does not match the final "
                            f"{output_format.upper()} model"
                        )
            if audit.get("sequence_search_performed"):
                raw_output = working / f"outfile{run.name}.csv"
                execution_path = working / "gapfill_sequence_search.json"
                for path in (raw_output, execution_path):
                    if not path.is_file():
                        raise RuntimeError(
                            "Gap-fill sequence search is missing its auditable "
                            f"raw artifact: {path}"
                        )
                    records[str(path.relative_to(self.repo_root))] = (
                        self._file_record(path)
                    )
                execution = json.loads(
                    execution_path.read_text(encoding="utf-8")
                )
                recorded_output = execution.get("raw_output", {})
                observed_output = self._file_record(raw_output)
                if not self._same_file_identity(
                    recorded_output,
                    observed_output,
                ):
                    raise RuntimeError(
                        "Gap-fill raw BLAST output does not match its "
                        "execution record"
                    )
        return records

    def _validate_structure_staging_audit(self, run: PipelineRun) -> None:
        working = self.repo_root / "working" / run.name
        audit_path = working / "structure_staging_audit.json"
        mapping_path = working / f"{run.name}.xlsx"
        stage_dir = working / "structures"
        try:
            audit = json.loads(audit_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise RuntimeError(
                f"Invalid structure staging audit: {audit_path}"
            ) from error
        required = {
            "allow_missing_structures",
            "mapping_row_count",
            "missing_structure_count",
            "missing_structures",
            "staged_structure_count",
            "staged_structures",
        }
        if not required.issubset(audit):
            raise RuntimeError(
                "Structure staging audit is missing required fields: "
                f"{sorted(required - set(audit))}"
            )
        pandas = importlib.import_module("pandas")
        mapping = pandas.read_excel(mapping_path)
        mapping_rows = len(mapping)
        missing = audit["missing_structures"]
        staged = audit["staged_structures"]
        if (
            audit["allow_missing_structures"] is not True
            or audit["mapping_row_count"] != mapping_rows
            or audit["missing_structure_count"] != len(missing)
            or audit["staged_structure_count"] != len(staged)
            or audit["mapping_row_count"]
            != audit["missing_structure_count"]
            + audit["staged_structure_count"]
        ):
            raise RuntimeError(
                "Structure staging audit violates its count conservation "
                "contract"
            )
        structure_columns = [
            column
            for column in ("Structure", "Structure Data")
            if column in mapping.columns
        ]
        if len(structure_columns) != 1:
            raise RuntimeError(
                "Structure staging mapping must have one structure column"
            )
        structure_column = structure_columns[0]
        source_root = self._resolve_input(run.structure)
        expected_missing = []
        expected_staged = []
        expected_sources = {}
        for row_index, row in mapping.iterrows():
            entry = str(row["Entry"]).strip()
            structure_value = str(row[structure_column]).strip()
            filename = f"AF-{entry}-F1-model_v4.pdb"
            basename = Path(structure_value).name
            candidate_names = [basename]
            if not basename.lower().endswith(".pdb"):
                candidate_names.extend(
                    [
                        f"{basename}.pdb",
                        f"AF-{basename}-F1-model_v4.pdb",
                    ]
                )
            source = next(
                (
                    source_root / candidate
                    for candidate in candidate_names
                    if (source_root / candidate).is_file()
                ),
                None,
            )
            if source is not None:
                expected_staged.append(
                    {"entry": entry, "file": filename}
                )
                expected_sources[filename] = source
            else:
                expected_missing.append(
                    {
                        "entry": entry,
                        "row_index": int(row_index),
                        "structure_value": structure_value,
                    }
                )
        observed_missing = [
            {
                "entry": record.get("entry"),
                "row_index": record.get("row_index"),
                "structure_value": record.get("structure_value"),
            }
            for record in missing
        ]
        observed_staged = [
            {
                "entry": record.get("entry"),
                "file": record.get("file"),
            }
            for record in staged
        ]
        if (
            observed_missing != expected_missing
            or observed_staged != expected_staged
        ):
            raise RuntimeError(
                "Structure staging audit does not match the source partition"
            )
        staged_names = [record.get("file") for record in staged]
        if (
            any(not isinstance(name, str) or not name for name in staged_names)
            or len(staged_names) != len(set(staged_names))
            or {
                path.name for path in stage_dir.glob("*.pdb")
            }
            != set(staged_names)
        ):
            raise RuntimeError(
                "Structure staging audit does not bind the staged PDB set"
            )
        for record in staged:
            path = stage_dir / record["file"]
            expected = {
                "bytes": record.get("bytes"),
                "sha256": record.get("sha256"),
            }
            if self._file_record(path) != expected:
                raise RuntimeError(
                    "Structure staging audit has a stale PDB record: "
                    f"{path}"
                )
            if self._file_record(path) != self._file_record(
                expected_sources[record["file"]]
            ):
                raise RuntimeError(
                    "Staged PDB does not match its resolved source bytes: "
                    f"{path}"
                )

    def _validate_usalign_audit(self, run: PipelineRun) -> None:
        working = self.repo_root / "working" / run.name
        audit_path = working / "usalign_execution.json"
        input_path = working / f"matrix_orthofinder{run.name}.xlsx"
        matrix_path = working / f"matrix_USalign_filtered{run.name}.xlsx"
        try:
            audit = json.loads(audit_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise RuntimeError(
                f"Invalid US-align execution audit: {audit_path}"
            ) from error
        required = {
            "failed_pair_count",
            "failures",
            "matrix_output",
            "orthofinder_input",
            "successful_pair_count",
            "total_pair_count",
        }
        if not required.issubset(audit):
            raise RuntimeError(
                "US-align audit is missing required fields: "
                f"{sorted(required - set(audit))}"
            )
        failures = audit["failures"]
        if (
            audit["failed_pair_count"] != len(failures)
            or audit["total_pair_count"]
            != audit["successful_pair_count"] + audit["failed_pair_count"]
            or audit["matrix_output"] != self._file_record(matrix_path)
            or audit["orthofinder_input"] != self._file_record(input_path)
        ):
            raise RuntimeError(
                "US-align audit violates its artifact or count contract"
            )
        pandas = importlib.import_module("pandas")
        matrix = pandas.read_excel(matrix_path)
        orthofinder = pandas.read_excel(input_path)
        reference_table = pandas.read_excel(
            self.repo_root / "data_available" / f"{run.refname}.xlsx"
        )
        reference_map = {
            str(reference_table.iat[index, 2]): str(
                reference_table.iat[index, 0]
            )
            for index in range(len(reference_table.index))
        }
        if (
            len(matrix) != audit["total_pair_count"]
            or len(orthofinder) != audit["total_pair_count"]
            or matrix.shape[1] < 8
        ):
            raise RuntimeError(
                "US-align audit row counts do not match its matrix inputs"
            )
        matrix_pairs = matrix.iloc[:, 1:3].astype(str).to_numpy().tolist()
        orthofinder_pairs = (
            orthofinder.iloc[:, 1:3].astype(str).to_numpy().tolist()
        )
        if matrix_pairs != orthofinder_pairs:
            raise RuntimeError(
                "US-align output does not preserve OrthoFinder pair order"
            )
        metrics = matrix.iloc[:, 3:8].apply(
            pandas.to_numeric,
            errors="coerce",
        )
        any_null = metrics.isna().any(axis=1)
        all_null = metrics.isna().all(axis=1)
        if (any_null & ~all_null).any():
            raise RuntimeError(
                "US-align matrix contains partially-null metric rows"
            )
        null_metric_rows = set(
            int(index)
            for index in metrics.index[all_null]
        )
        failure_rows = [
            failure.get("row_index") for failure in failures
        ]
        if (
            any(
                not isinstance(index, int)
                or index < 0
                or index >= len(matrix)
                for index in failure_rows
            )
            or len(failure_rows) != len(set(failure_rows))
            or set(failure_rows) != null_metric_rows
        ):
            raise RuntimeError(
                "US-align audit failures do not match matrix null-metric rows"
            )
        for failure in failures:
            index = failure["row_index"]
            expected_resolved_reference = reference_map.get(
                str(failure.get("reference_id")),
                str(failure.get("reference_id")),
            )
            if (
                failure.get("reason") != "missing_structure"
                or str(matrix.iat[index, 1])
                != str(failure.get("reference_id"))
                or str(matrix.iat[index, 2])
                != str(failure.get("target_id"))
            ):
                raise RuntimeError(
                    "US-align audit failure metadata does not match its "
                    f"matrix row: {index}"
                )
            if (
                str(failure.get("resolved_reference_id"))
                != expected_resolved_reference
            ):
                raise RuntimeError(
                    "US-align audit has an invalid resolved reference ID "
                    f"at row {index}"
                )
            missing_structures = failure.get("missing_structures")
            if (
                not isinstance(missing_structures, list)
                or not missing_structures
                or not set(missing_structures)
                <= {"reference", "target"}
            ):
                raise RuntimeError(
                    "US-align missing-structure failure has invalid details"
                )
            target_path = (
                working
                / "structures"
                / f"AF-{failure['target_id']}-F1-model_v4.pdb"
            )
            reference_path = (
                self.repo_root
                / "struct_data"
                / run.refname
                / (
                    "AF-{}-F1-model_v4.pdb".format(
                        failure["resolved_reference_id"]
                    )
                )
            )
            actual_missing = []
            if not reference_path.is_file():
                actual_missing.append("reference")
            if not target_path.is_file():
                actual_missing.append("target")
            if actual_missing != missing_structures:
                raise RuntimeError(
                    "US-align missing-structure audit does not match "
                    f"filesystem state at row {index}"
                )

    def _validate_prottrans_audit(self, run: PipelineRun) -> None:
        working = self.repo_root / "working" / run.name
        audit_path = working / "prottrans_execution.json"
        embedding_path = working / f"{run.name}.emb"
        matrix_path = working / f"matrix_homolog{run.name}.xlsx"
        try:
            audit = json.loads(audit_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise RuntimeError(
                f"Invalid ProtTrans execution audit: {audit_path}"
            ) from error
        required = {
            "embedding_count",
            "embedding_output",
            "fallback_reason",
            "final_device",
            "final_dtype",
            "group_count",
            "matrix_output",
            "model_source",
            "requested_device",
            "status",
            "torch_version",
            "transformers_version",
            "uniform_embedding_precision",
        }
        if not required.issubset(audit):
            raise RuntimeError(
                "ProtTrans audit is missing required fields: "
                f"{sorted(required - set(audit))}"
            )
        if (
            audit["status"] != "completed"
            or audit["uniform_embedding_precision"] is not True
            or audit["embedding_count"] < 0
            or audit["group_count"] < 0
            or audit["embedding_output"] != self._file_record(embedding_path)
            or audit["matrix_output"] != self._file_record(matrix_path)
        ):
            raise RuntimeError(
                "ProtTrans audit violates its output contract"
            )
        expected_model_source = str(
            self._resolve_input(run.prottrans_model_path).resolve(
                strict=True
            )
        )
        runtime_environment = (
            self._run_runtime_environment or self._runtime_environment()
        )
        packages = runtime_environment["packages"]
        if audit["model_source"] != expected_model_source:
            raise RuntimeError(
                "ProtTrans audit does not bind the configured model source"
            )
        if (
            audit["torch_version"] != packages["torch"]
            or audit["transformers_version"] != packages["transformers"]
        ):
            raise RuntimeError(
                "ProtTrans audit does not bind the formal runtime versions"
            )
        final_device = audit["final_device"]
        final_dtype = audit["final_dtype"]
        if (final_device, final_dtype) not in {
            ("cuda", "float16"),
            ("cpu", "float32"),
        }:
            raise RuntimeError(
                "ProtTrans audit has an invalid uniform device/dtype pair"
            )
        fallback_reason = audit["fallback_reason"]
        if fallback_reason not in {
            None,
            "cuda_inference_oom_full_cpu_recompute",
            "cuda_model_transfer_oom",
        }:
            raise RuntimeError(
                "ProtTrans audit has an unsupported fallback reason"
            )
        requested_device = audit["requested_device"]
        if requested_device not in {"cpu", "cuda"}:
            raise RuntimeError(
                "ProtTrans audit has an invalid requested device"
            )
        if fallback_reason is not None and final_device != "cpu":
            raise RuntimeError(
                "ProtTrans fallback must finish with uniform CPU embeddings"
            )
        if requested_device == "cpu" and (
            final_device != "cpu" or fallback_reason is not None
        ):
            raise RuntimeError(
                "ProtTrans CPU request has an invalid execution state"
            )
        if requested_device == "cuda" and (
            (fallback_reason is None and final_device != "cuda")
            or (fallback_reason is not None and final_device != "cpu")
        ):
            raise RuntimeError(
                "ProtTrans CUDA request changed device without an audit "
                "fallback reason"
            )

    def _validate_foldseek_audit(self, run: PipelineRun) -> None:
        working = self.repo_root / "working" / run.name
        paths = {
            "audit": working / "foldseek_execution.json",
            "raw_output": working / f"aln{run.name}.csv",
            "matrix_output": working / f"matrix_foldseek_{run.name}.xlsx",
        }
        for logical_name, path in paths.items():
            self._assert_repo_local_path(
                path,
                label=f"Formal Foldseek {logical_name}",
                require_file=True,
            )
        try:
            audit = json.loads(paths["audit"].read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise RuntimeError(
                f"Formal Foldseek stage has an invalid execution audit: "
                f"{paths['audit']}"
            ) from error
        required = {
            "schema_version",
            "started_at_utc",
            "completed_at_utc",
            "allow_score_fallback",
            "score_field",
            "createdb_command",
            "easy_search_command",
            "raw_output",
            "matrix_output",
        }
        if audit.get("schema_version") != 1 or not required.issubset(audit):
            raise RuntimeError(
                "Formal Foldseek execution audit is incomplete or has an "
                "unsupported schema"
            )
        if (
            audit["allow_score_fallback"] is not False
            or audit["score_field"] != "alntmscore"
        ):
            raise RuntimeError(
                "Formal Foldseek execution changed the fixed alntmscore "
                "semantics"
            )
        createdb = audit["createdb_command"]
        command = audit["easy_search_command"]
        if (
            not isinstance(createdb, list)
            or createdb[:2] != ["foldseek", "createdb"]
            or not isinstance(command, list)
            or command[:2] != ["foldseek", "easy-search"]
            or "--format-output" not in command
        ):
            raise RuntimeError(
                "Formal Foldseek execution audit has invalid commands"
            )
        format_index = command.index("--format-output")
        if (
            format_index + 1 >= len(command)
            or command[format_index + 1]
            != "query,target,alntmscore,prob,qcov,tcov"
            or "qtmscore" in command
        ):
            raise RuntimeError(
                "Formal Foldseek execution audit does not prove the fixed "
                "alntmscore output field"
            )
        try:
            started = datetime.fromisoformat(audit["started_at_utc"])
            completed = datetime.fromisoformat(audit["completed_at_utc"])
        except (TypeError, ValueError) as error:
            raise RuntimeError(
                "Formal Foldseek execution audit has invalid UTC timestamps"
            ) from error
        if (
            started.tzinfo is None
            or completed.tzinfo is None
            or completed < started
        ):
            raise RuntimeError(
                "Formal Foldseek execution audit has inconsistent timestamps"
            )
        for logical_name in ("raw_output", "matrix_output"):
            observed = self._file_record(paths[logical_name])
            if not self._same_file_identity(audit[logical_name], observed):
                raise RuntimeError(
                    "Formal Foldseek execution audit does not match the fresh "
                    f"{logical_name}"
                )

    def _validate_clean_audit(self, run: PipelineRun) -> None:
        working = self.repo_root / "working" / run.name
        clean_root = self.repo_root / "tools" / "CLEAN" / "app"
        audit_path = working / "clean_execution.json"
        try:
            audit = json.loads(audit_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise RuntimeError(
                f"Formal CLEAN stage has an invalid execution audit: {audit_path}"
            ) from error
        required = {
            "schema_version",
            "command",
            "command_cwd",
            "started_at_utc",
            "completed_at_utc",
            "removed_preexisting_result",
            "isolated_runtime",
            "runtime_directory",
            "static_links",
            "nested_python",
            "esm_weights",
            "esm_weight_storage",
            "esm_extraction",
            "torch_home",
            "network_policy",
            "pythonpath",
            "python_no_user_site",
            "clean_entrypoint",
            "input_fasta",
            "working_input",
            "clean_result",
            "working_result",
        }
        if audit.get("schema_version") != 3 or not required.issubset(audit):
            raise RuntimeError(
                "Formal CLEAN execution audit is incomplete or has an "
                "unsupported schema"
            )
        runtime = Path(audit["runtime_directory"])
        self._assert_repo_local_path(
            runtime,
            label="Formal CLEAN isolated runtime",
        )
        try:
            runtime_relative = runtime.relative_to(working)
        except ValueError as error:
            raise RuntimeError(
                "Formal CLEAN runtime is outside its working directory: "
                f"{runtime}"
            ) from error
        if (
            runtime.is_symlink()
            or not runtime.is_dir()
            or len(runtime_relative.parts) != 1
            or not runtime_relative.name.startswith("clean-runtime-")
            or audit.get("isolated_runtime") is not True
        ):
            raise RuntimeError(
                f"Formal CLEAN isolated runtime is invalid: {runtime}"
            )
        paths = {
            "entrypoint": clean_root / "CLEAN_infer_fasta.py",
            "clean_input": (
                runtime
                / "data"
                / "inputs"
                / f"{run.name}_homoleft.fasta"
            ),
            "working_input": working / f"{run.name}_homoleft.fasta",
            "clean_result": (
                runtime
                / "results"
                / "inputs"
                / f"{run.name}_homoleft_maxsep.csv"
            ),
            "working_result": working / f"{run.name}_homoleft_maxsep.csv",
        }
        expected_command = [
            sys.executable,
            str(paths["entrypoint"]),
            "--fasta_data",
            f"{run.name}_homoleft",
        ]
        if (
            audit["command"] != expected_command
            or audit["command_cwd"] != str(runtime)
        ):
            raise RuntimeError(
                "Formal CLEAN execution audit does not record the exact "
                "expected command and isolated working directory"
            )
        try:
            started = datetime.fromisoformat(audit["started_at_utc"])
            completed = datetime.fromisoformat(audit["completed_at_utc"])
        except (TypeError, ValueError) as error:
            raise RuntimeError(
                "Formal CLEAN execution audit has invalid UTC timestamps"
            ) from error
        if (
            started.tzinfo is None
            or completed.tzinfo is None
            or completed < started
        ):
            raise RuntimeError(
                "Formal CLEAN execution audit has inconsistent timestamps"
            )
        if audit["removed_preexisting_result"] is not False:
            raise RuntimeError(
                "An isolated CLEAN runtime cannot remove a pre-existing result"
            )
        static_links = audit["static_links"]
        if not isinstance(static_links, list) or not static_links:
            raise RuntimeError(
                "Formal CLEAN execution audit lacks its static-link contract"
            )
        for record in static_links:
            if not isinstance(record, dict) or set(record) != {
                "path",
                "target",
            }:
                raise RuntimeError(
                    "Formal CLEAN execution audit has an invalid static link"
                )
            link = Path(record["path"])
            target = Path(record["target"])
            try:
                link.relative_to(runtime)
                target.relative_to(clean_root)
            except ValueError as error:
                raise RuntimeError(
                    "Formal CLEAN static link escapes its allowed roots: "
                    f"{link} -> {target}"
                ) from error
            if (
                not link.is_symlink()
                or link.resolve(strict=True) != target
                or target != target.resolve(strict=True)
            ):
                raise RuntimeError(
                    f"Formal CLEAN static link changed: {link}"
                )
        nested_python = audit["nested_python"]
        active_python = Path(sys.executable).resolve()
        active_invocation = Path(sys.executable).absolute()
        invocation = (
            Path(nested_python.get("invocation_path", ""))
            if isinstance(nested_python, dict)
            else Path()
        )
        if (
            not isinstance(nested_python, dict)
            or invocation.parent != active_invocation.parent
            or invocation.resolve(strict=True) != active_python
            or nested_python.get("resolved_path") != str(active_python)
            or nested_python.get("sys_prefix") != sys.prefix
            or nested_python
            != {
                "invocation_path": str(invocation),
                "resolved_path": str(active_python),
                "sys_prefix": sys.prefix,
                **self._runtime_file_record(active_python),
            }
        ):
            raise RuntimeError(
                "Formal CLEAN nested Python is not bound to the active "
                "interpreter"
            )
        weight_records = audit["esm_weights"]
        if not isinstance(weight_records, dict) or set(weight_records) != {
            "source_model",
            "source_regression",
        }:
            raise RuntimeError(
                "Formal CLEAN execution audit has an incomplete ESM weight "
                "contract"
            )
        source_model_path = Path(
            weight_records["source_model"].get("observed_path", "")
        )
        if (
            source_model_path.name != "esm1b_t33_650M_UR50S.pt"
            or source_model_path.parent.name != "checkpoints"
            or source_model_path.parent.parent.name != "hub"
        ):
            raise RuntimeError(
                "Formal CLEAN source ESM weight is outside the pinned "
                "TORCH_HOME/hub/checkpoints layout"
            )
        expected_torch_home = source_model_path.parent.parent.parent
        if (
            audit["torch_home"] != str(expected_torch_home)
            or audit["network_policy"] != "loopback_refuse_proxy"
            or audit["pythonpath"]
            != os.pathsep.join(
                (
                    str(clean_root / "src"),
                    str(clean_root / "esm"),
                )
            )
            or audit["python_no_user_site"] != "1"
        ):
            raise RuntimeError(
                "Formal CLEAN execution audit lacks its isolated offline "
                "model-loading contract"
            )
        esm_extraction = audit["esm_extraction"]
        canonical_extract = clean_root / "esm" / "scripts" / "extract.py"
        source_wrapper = (
            self.repo_root / "src" / "clean_esm_extract_wrapper.py"
        )
        runtime_wrapper = runtime / "esm" / "scripts" / "extract.py"
        if (
            not isinstance(esm_extraction, dict)
            or set(esm_extraction) != {
                "token_budget",
                "canonical_extract",
                "source_wrapper",
                "runtime_wrapper",
            }
            or esm_extraction["token_budget"] != run.esm_token_budget
            or esm_extraction["canonical_extract"]
            != {
                "observed_path": str(canonical_extract),
                **self._runtime_file_record(canonical_extract),
            }
            or esm_extraction["source_wrapper"]
            != {
                "observed_path": str(source_wrapper),
                **self._runtime_file_record(source_wrapper),
            }
            or esm_extraction["runtime_wrapper"]
            != {
                "observed_path": str(runtime_wrapper),
                **self._runtime_file_record(runtime_wrapper),
            }
            or not self._same_file_identity(
                esm_extraction["runtime_wrapper"],
                esm_extraction["source_wrapper"],
            )
        ):
            raise RuntimeError(
                "Formal CLEAN execution audit has an invalid ESM token-budget "
                "controller"
            )
        if audit["esm_weight_storage"] != "direct_pinned_source":
            raise RuntimeError(
                "Formal CLEAN must use the pinned source ESM weights directly"
            )
        preflight_resources = self.preflight_by_stage.get("clean", {})
        for logical_name, audit_name in (
            ("clean_esm_model_weights", "source_model"),
            ("clean_esm_regression_weights", "source_regression"),
        ):
            record = weight_records[audit_name]
            path = Path(record["observed_path"])
            if (
                not path.is_file()
                or path.is_symlink()
                or not self._same_file_identity(
                    preflight_resources.get(logical_name),
                    self._file_record(path),
                )
                or record
                != {
                    "observed_path": str(path),
                    **self._runtime_file_record(path),
                }
            ):
                raise RuntimeError(
                    "Formal CLEAN source ESM weight does not match "
                    f"preflight: {logical_name}"
                )
        for logical_name in (
            "entrypoint",
            "clean_input",
            "working_input",
            "clean_result",
            "working_result",
        ):
            path = paths[logical_name]
            self._assert_repo_local_path(
                path,
                label=f"Formal CLEAN {logical_name}",
                require_file=True,
            )
            if (
                not path.is_file()
                or path.is_symlink()
                or path.stat().st_size == 0
            ):
                raise RuntimeError(
                    "Formal CLEAN execution evidence is missing, empty, or a "
                    f"symlink: {logical_name}={path}"
                )
        observed = {}
        for logical_name in (
            "entrypoint",
            "clean_input",
            "working_input",
            "clean_result",
            "working_result",
        ):
            if logical_name in {"entrypoint", "clean_input", "clean_result"}:
                observed[logical_name] = self._runtime_file_record(
                    paths[logical_name]
                )
            else:
                observed[logical_name] = self._file_record(paths[logical_name])
        recorded_entrypoint = (
            self.preflight_by_stage.get("clean", {}).get("clean_entrypoint")
        )
        if (
            not self._same_file_identity(
                audit["clean_entrypoint"],
                observed["entrypoint"],
            )
            or audit["clean_entrypoint"] != observed["entrypoint"]
            or not self._same_file_identity(
                recorded_entrypoint,
                observed["entrypoint"],
            )
            or not self._same_file_identity(
                audit["input_fasta"],
                observed["clean_input"],
            )
            or audit["input_fasta"] != observed["clean_input"]
            or not self._same_file_identity(
                audit["working_input"],
                observed["working_input"],
            )
            or audit["working_input"] != observed["working_input"]
            or not self._same_file_identity(
                observed["clean_input"],
                observed["working_input"],
            )
            or not self._same_file_identity(
                audit["clean_result"],
                observed["clean_result"],
            )
            or audit["clean_result"] != observed["clean_result"]
            or not self._same_file_identity(
                audit["working_result"],
                observed["working_result"],
            )
            or not self._same_file_identity(
                observed["clean_result"],
                observed["working_result"],
            )
        ):
            raise RuntimeError(
                "Formal CLEAN execution evidence does not bind the preflight "
                "entrypoint, current non-homolog input, and fresh result"
            )

    def _validate_stage_audit(self, stage: str, run: PipelineRun) -> None:
        working = self.repo_root / "working" / run.name
        audit_paths = {
            "draft_model": working / "draft_model_audit.json",
            "merge_annotations": working / "annotation_merge_audit.json",
            "reaction_pool": working / "reaction_pool_audit.json",
        }
        path = audit_paths[stage]
        try:
            audit = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise RuntimeError(
                f"Stage {stage!r} has an invalid audit artifact: {path}"
            ) from error
        if stage == "draft_model":
            required = {
                "homology_reference_gene_count",
                "mapped_homology_reference_genes",
                "mapped_homology_reference_gene_count",
                "missing_homology_reference_genes",
                "missing_homology_reference_gene_count",
                "model1_reaction_count",
                "complex_candidate_reaction_count",
                "complex_retained_reaction_count",
                "missing_complex_subunits",
                "missing_complex_subunit_count",
                "model2_reaction_count",
            }
            if not required.issubset(audit):
                raise RuntimeError(
                    "Draft-model audit is missing required fields: "
                    f"{sorted(required - set(audit))}"
                )
            if (
                audit["homology_reference_gene_count"]
                != audit["mapped_homology_reference_gene_count"]
                + audit["missing_homology_reference_gene_count"]
                or audit["mapped_homology_reference_gene_count"]
                != len(audit["mapped_homology_reference_genes"])
                or audit["missing_homology_reference_gene_count"]
                != len(audit["missing_homology_reference_genes"])
                or set(audit["mapped_homology_reference_genes"])
                & set(audit["missing_homology_reference_genes"])
                or audit["missing_complex_subunit_count"]
                != len(audit["missing_complex_subunits"])
                or audit["complex_retained_reaction_count"]
                > audit["complex_candidate_reaction_count"]
                or audit["model1_reaction_count"] <= 0
                or audit["model2_reaction_count"] <= 0
            ):
                raise RuntimeError(
                    "Draft-model audit violates its conservation contract"
                )
            return
        if stage == "merge_annotations":
            required = {
                "eggnog_candidate_pair_count",
                "eggnog_rhea_mapped_pairs",
                "eggnog_unmapped_pairs",
                "eggnog_unmapped_pair_count",
                "merged_annotation_pair_count",
            }
            if not required.issubset(audit):
                raise RuntimeError(
                    "Annotation-merge audit is missing required fields: "
                    f"{sorted(required - set(audit))}"
                )
            if (
                audit["eggnog_candidate_pair_count"]
                != audit["eggnog_rhea_mapped_pairs"]
                + audit["eggnog_unmapped_pair_count"]
                or audit["eggnog_unmapped_pair_count"]
                != len(audit["eggnog_unmapped_pairs"])
                or audit["merged_annotation_pair_count"] < 0
            ):
                raise RuntimeError(
                    "Annotation-merge audit violates candidate conservation"
                )
            return
        allowed_statuses = {
            "added",
            "already_present",
            "already_present_conflict",
            "already_present_gpr_merged",
            "already_present_no_merge",
            "explicitly_excluded_released_code",
            "incomplete_stoichiometry",
            "missing_from_reaction_pool",
            "rejected_mass_or_charge_balance",
            "rejected_postmerge_metabolite_identity_conflict",
        }
        candidates = audit.get("candidates")
        status_counts = audit.get("status_counts")
        if not isinstance(candidates, list) or not isinstance(
            status_counts,
            dict,
        ):
            raise RuntimeError(
                "Reaction-pool audit lacks candidate/status records"
            )
        observed_counts = {}
        for record in candidates:
            if not isinstance(record, dict):
                raise RuntimeError(
                    "Reaction-pool audit has a non-object candidate record"
                )
            status = record.get("status")
            if status not in allowed_statuses:
                raise RuntimeError(
                    f"Reaction-pool audit has invalid status: {status!r}"
                )
            observed_counts[status] = observed_counts.get(status, 0) + 1
        if (
            type(audit.get("input_candidate_count")) is not int
            or any(type(value) is not int or value < 0 for value in status_counts.values())
            or
            audit.get("input_candidate_count") != len(candidates)
            or status_counts != observed_counts
            or sum(status_counts.values()) != len(candidates)
        ):
            raise RuntimeError(
                "Reaction-pool audit violates candidate conservation"
            )
        if self._is_exact_paper_profile(run):
            disallowed = {
                status: observed_counts.get(status, 0)
                for status in ("explicitly_excluded_released_code",)
                if observed_counts.get(status, 0)
            }
            if disallowed:
                raise RuntimeError(
                    "Formal paper reaction-pool audit contains candidates "
                    "that were not evaluated by the published QC contract: "
                    f"{disallowed}"
                )
            self._validate_formal_reaction_pool_audit(run, path, audit)

    def _validate_formal_reaction_pool_audit(
        self,
        run: PipelineRun,
        audit_path: Path,
        audit: dict,
    ) -> None:
        """Bind skipped and chemically rejected candidates to source files."""
        working = self.repo_root / "working" / run.name
        pool_path = self.repo_root / "models" / "reaction_pool_final.json"
        preflight = self.preflight_by_stage.get("reaction_pool", {}).get(
            "reaction_pool"
        )
        if (
            pool_path.is_symlink()
            or not pool_path.is_file()
            or not isinstance(preflight, dict)
            or preflight.get("observed_path") != str(pool_path)
        ):
            raise RuntimeError(
                "Formal reaction-pool audit lacks the frozen pool preflight"
            )
        pool_bytes = pool_path.read_bytes()
        pool_record = {
            "bytes": len(pool_bytes),
            "sha256": hashlib.sha256(pool_bytes).hexdigest(),
        }
        if (
            pool_record["sha256"] != FROZEN_REACTION_POOL_SHA256
            or not self._same_file_identity(preflight, pool_record)
        ):
            raise RuntimeError(
                "Formal reaction-pool audit pool differs from its frozen "
                "snapshot or preflight identity"
            )
        try:
            pool = json.loads(pool_bytes)
            pool_reactions = pool["reactions"]
        except (ValueError, KeyError, TypeError) as error:
            raise RuntimeError(
                "Formal reaction-pool audit cannot parse the frozen pool"
            ) from error
        if not isinstance(pool_reactions, list):
            raise RuntimeError("Frozen reaction pool lacks reaction records")
        pool_ids = set()
        pool_reactions_by_id = {}
        for reaction in pool_reactions:
            reaction_id = reaction.get("id") if isinstance(reaction, dict) else None
            if not isinstance(reaction_id, str) or not reaction_id:
                raise RuntimeError("Frozen reaction pool has an invalid reaction ID")
            if reaction_id in pool_ids:
                raise RuntimeError("Frozen reaction pool has duplicate reaction IDs")
            pool_ids.add(reaction_id)
            pool_reactions_by_id[reaction_id] = reaction
        pool_metabolites_by_id = {
            metabolite["id"]: metabolite
            for metabolite in pool.get("metabolites", [])
            if isinstance(metabolite, dict) and isinstance(metabolite.get("id"), str)
        }
        if len(pool_metabolites_by_id) != len(pool.get("metabolites", [])):
            raise RuntimeError("Frozen reaction pool has invalid metabolite IDs")
        if self._file_record(pool_path) != pool_record:
            raise RuntimeError("Frozen reaction pool changed during audit")

        # The author iterates one spreadsheet row per candidate. Internal
        # status counts alone cannot prove that no input row was omitted.
        score_path = working / f"{run.name}gpr_score.xlsx"
        score_rows = self._module("pandas").read_excel(score_path)
        candidates = audit["candidates"]
        if len(score_rows) != len(candidates):
            raise RuntimeError(
                "Formal reaction-pool audit does not conserve input rows"
            )
        missing = []
        identity_rejections = []
        needs_model = []
        draft_path = working / f"model2{run.name}.yml"
        cobra = self._module("cobra")
        model_path = working / f"{run.name}-GEM_withgaps.yml"
        model = cobra.io.load_yaml_model(str(model_path))
        reference_path = self.repo_root / "models" / REFERENCE_MODELS[run.refname]
        reference_preflight = self.preflight_by_stage.get("reaction_pool", {}).get(
            "reference_model"
        )
        reference_record = self._file_record(reference_path)
        if (
            reference_path.is_symlink()
            or not isinstance(reference_preflight, dict)
            or reference_preflight.get("observed_path") != str(reference_path)
            or not self._same_file_identity(reference_preflight, reference_record)
        ):
            raise RuntimeError(
                "Formal reaction-pool audit reference differs from preflight"
            )
        reference_model = None

        def load_reference_model():
            nonlocal reference_model
            if reference_model is None:
                reference_model = cobra.io.read_sbml_model(str(reference_path))
            return reference_model

        already_present_policy = run.reaction_pool_already_present_policy
        conflict_policy = run.reaction_pool_conflict_policy
        if already_present_policy not in REACTION_POOL_ALREADY_PRESENT_POLICIES:
            raise RuntimeError(
                "Unknown already-present policy: "
                f"{already_present_policy!r}"
            )
        if conflict_policy not in REACTION_POOL_CONFLICT_POLICIES:
            raise RuntimeError(f"Unknown conflict policy: {conflict_policy!r}")
        for record in candidates:
            status = record["status"]
            if (
                already_present_policy == "author_skip"
                and status
                in {
                    "already_present_conflict",
                    "already_present_gpr_merged",
                    "already_present_no_merge",
                }
            ):
                raise RuntimeError(
                    "Reaction-pool audit records already-present merge states "
                    "under the author_skip policy"
                )
            if (
                already_present_policy
                not in {
                    "merge_verified_identical_or_multisource",
                    "merge_threeway_or_multisource",
                }
                and status == "already_present_no_merge"
            ):
                raise RuntimeError(
                    "Reaction-pool audit records no-merge states without a "
                    "multisource policy"
                )
            if (
                status == "already_present_conflict"
                and record.get("conflict_kind") == "template_pool_divergence"
                and already_present_policy != "merge_threeway_or_multisource"
            ):
                raise RuntimeError(
                    "Reaction-pool audit records a template/pool divergence "
                    "without the three-way policy"
                )
            if (
                already_present_policy != "author_skip"
                and status == "already_present"
            ):
                raise RuntimeError(
                    "Reaction-pool audit records an unmerged already-present "
                    f"skip under policy {already_present_policy!r}"
                )
            if conflict_policy != "suffix_new_metabolite" and (
                "suffix_metabolites" in record
                or record.get("addition_kind") == "suffix_metabolite"
            ):
                raise RuntimeError(
                    "Reaction-pool audit records suffix-metabolite additions "
                    f"under conflict policy {conflict_policy!r}"
                )
            if conflict_policy != "remap_by_model_identity" and (
                "remapped_metabolites" in record
                or record.get("addition_kind") == "remapped_metabolite"
            ):
                raise RuntimeError(
                    "Reaction-pool audit records remapped-metabolite additions "
                    f"under conflict policy {conflict_policy!r}"
                )
            if record.get("addition_kind") and status != "added":
                raise RuntimeError(
                    "Reaction-pool audit records an addition kind for a "
                    f"non-added candidate: {record.get('reaction_id')}"
                )

        def normalized_names(name):
            lowered = (name or "").lower()
            return {lowered, lowered.replace("(", "").replace(")", "")}

        def candidate_source_summary(reaction_id, source):
            """Coefficients, identities, names and bounds of a record's source."""
            if source == "frozen_reaction_pool":
                source_reaction = pool_reactions_by_id.get(reaction_id)
                if source_reaction is None:
                    raise RuntimeError(
                        "Candidate source is absent from the frozen pool: "
                        f"{reaction_id}"
                    )
                source_metabolites = source_reaction.get("metabolites")
                if not isinstance(source_metabolites, dict):
                    raise RuntimeError(
                        "Frozen reaction-pool candidate has no stoichiometry"
                    )
                for metabolite_id in source_metabolites:
                    if metabolite_id not in pool_metabolites_by_id:
                        raise RuntimeError(
                            "Frozen reaction pool lacks metabolite "
                            f"{metabolite_id}"
                        )
                bounds = (
                    source_reaction.get("lower_bound"),
                    source_reaction.get("upper_bound"),
                )
                if any(
                    not isinstance(bound, (int, float)) for bound in bounds
                ):
                    raise RuntimeError(
                        "Frozen reaction-pool candidate lacks flux bounds: "
                        f"{reaction_id}"
                    )
                return {
                    "coefficients": {
                        metabolite_id: float(coefficient)
                        for metabolite_id, coefficient in source_metabolites.items()
                    },
                    "identities": {
                        metabolite_id: {
                            field: pool_metabolites_by_id[metabolite_id].get(field)
                            for field in ("formula", "charge", "compartment")
                        }
                        for metabolite_id in source_metabolites
                    },
                    "names": {
                        metabolite_id: pool_metabolites_by_id[metabolite_id].get(
                            "name"
                        )
                        for metabolite_id in source_metabolites
                    },
                    "lower_bound": float(bounds[0]),
                    "upper_bound": float(bounds[1]),
                }
            if source == "yeast9_reference":
                reference = load_reference_model()
                try:
                    reference_reaction = reference.reactions.get_by_id(reaction_id)
                except KeyError as error:
                    raise RuntimeError(
                        "Candidate source is absent from reference model"
                    ) from error
                return {
                    "coefficients": {
                        metabolite.id: float(coefficient)
                        for metabolite, coefficient in reference_reaction.metabolites.items()
                    },
                    "identities": {
                        metabolite.id: {
                            "formula": metabolite.formula,
                            "charge": metabolite.charge,
                            "compartment": metabolite.compartment,
                        }
                        for metabolite in reference_reaction.metabolites
                    },
                    "names": {
                        metabolite.id: metabolite.name
                        for metabolite in reference_reaction.metabolites
                    },
                    "lower_bound": float(reference_reaction.lower_bound),
                    "upper_bound": float(reference_reaction.upper_bound),
                }
            raise RuntimeError(f"Candidate lacks a known source: {source!r}")

        def saved_reaction_summary(reaction):
            return {
                "coefficients": {
                    metabolite.id: float(coefficient)
                    for metabolite, coefficient in reaction.metabolites.items()
                },
                "identities": {
                    metabolite.id: {
                        "formula": metabolite.formula,
                        "charge": metabolite.charge,
                        "compartment": metabolite.compartment,
                    }
                    for metabolite in reaction.metabolites
                },
                "names": {
                    metabolite.id: metabolite.name
                    for metabolite in reaction.metabolites
                },
                "lower_bound": float(reaction.lower_bound),
                "upper_bound": float(reaction.upper_bound),
            }

        def classify_existing_vs_source(existing, source):
            """Mirror the builder's already-present verification precedence."""
            existing_coefficients = existing["coefficients"]
            source_coefficients = source["coefficients"]
            if set(existing_coefficients) != set(source_coefficients):
                return "stoichiometry_or_identity_differ"
            for metabolite_id, coefficient in existing_coefficients.items():
                if not math.isclose(
                    coefficient,
                    source_coefficients[metabolite_id],
                    rel_tol=0.0,
                    abs_tol=1e-9,
                ):
                    return "stoichiometry_or_identity_differ"
            for metabolite_id in existing_coefficients:
                identities = (
                    existing["identities"][metabolite_id],
                    source["identities"][metabolite_id],
                )
                if any(
                    value is None or value == ""
                    for identity in identities
                    for value in identity.values()
                ):
                    return "unverifiable_identity"
            for metabolite_id in existing_coefficients:
                if (
                    existing["identities"][metabolite_id]
                    != source["identities"][metabolite_id]
                ):
                    return "stoichiometry_or_identity_differ"
            if not (
                math.isclose(
                    existing["lower_bound"],
                    source["lower_bound"],
                    rel_tol=0.0,
                    abs_tol=1e-9,
                )
                and math.isclose(
                    existing["upper_bound"],
                    source["upper_bound"],
                    rel_tol=0.0,
                    abs_tol=1e-9,
                )
            ):
                return "direction_or_bounds_differ"
            return None

        def rule_genes_and_complexity(rule):
            # cobra normalizes malformed rules to an empty GPR (body None)
            # with a warning instead of raising; mirror that here.
            if not isinstance(rule, str) or not rule.strip():
                return set(), False
            gpr = cobra.core.GPR.from_string(rule)
            if gpr.body is None:
                return set(), False
            has_complex = any(
                isinstance(node, ast.BoolOp) and isinstance(node.op, ast.And)
                for node in ast.walk(gpr.body)
            )
            return set(gpr.genes), has_complex

        def summary_nadp_family(name):
            """Mirror of the builder's heuristic cofactor-family label."""
            lowered = (name or "").lower()
            if "nadp" in lowered:
                return "nadp"
            if "nad" in lowered:
                return "nad"
            return None

        def summary_has_cofactor_pair(names_left, names_right):
            families_left = {
                summary_nadp_family(name) for name in names_left
            } - {None}
            families_right = {
                summary_nadp_family(name) for name in names_right
            } - {None}
            return (
                bool(families_left)
                and bool(families_right)
                and families_left != families_right
            )

        def classify_summary_divergence(left, right):
            """Mirror the builder's extended divergence classification."""
            left_coefficients = left["coefficients"]
            right_coefficients = right["coefficients"]
            if set(left_coefficients) != set(right_coefficients):
                only_left = sorted(set(left_coefficients) - set(right_coefficients))
                only_right = sorted(set(right_coefficients) - set(left_coefficients))
                detail = {
                    "left_only_metabolites": only_left,
                    "right_only_metabolites": only_right,
                }
                if summary_has_cofactor_pair(
                    [left["names"].get(m) for m in only_left],
                    [right["names"].get(m) for m in only_right],
                ):
                    return {"classification": "cofactor_pair", "detail": detail}
                return {
                    "classification": "stoichiometry_or_identity_differ",
                    "detail": detail,
                }
            coefficient_differences = sorted(
                metabolite_id
                for metabolite_id, coefficient in left_coefficients.items()
                if not math.isclose(
                    coefficient,
                    right_coefficients[metabolite_id],
                    rel_tol=0.0,
                    abs_tol=1e-9,
                )
            )
            if coefficient_differences:
                reversed_all = all(
                    math.isclose(
                        left_coefficients[m],
                        -right_coefficients[m],
                        rel_tol=0.0,
                        abs_tol=1e-9,
                    )
                    for m in left_coefficients
                ) and any(
                    not math.isclose(
                        left_coefficients[m], 0.0, rel_tol=0.0, abs_tol=1e-9
                    )
                    for m in left_coefficients
                )
                detail = {"coefficient_differences": coefficient_differences}
                if reversed_all:
                    return {
                        "classification": "reversed_stoichiometry",
                        "detail": detail,
                    }
                return {
                    "classification": "stoichiometry_or_identity_differ",
                    "detail": detail,
                }
            incomplete = []
            identity_differences = []
            for metabolite_id in sorted(left_coefficients):
                left_identity = left["identities"][metabolite_id]
                right_identity = right["identities"][metabolite_id]
                if any(
                    value is None or value == ""
                    for identity in (left_identity, right_identity)
                    for value in identity.values()
                ):
                    incomplete.append(metabolite_id)
                elif left_identity != right_identity:
                    identity_differences.append({
                        "metabolite_id": metabolite_id,
                        "left_identity": left_identity,
                        "right_identity": right_identity,
                    })
            if incomplete:
                return {
                    "classification": "unverifiable_identity",
                    "detail": {"incomplete_metabolites": incomplete},
                }
            if identity_differences:
                if summary_has_cofactor_pair(
                    [
                        left["names"].get(row["metabolite_id"])
                        for row in identity_differences
                    ],
                    [
                        right["names"].get(row["metabolite_id"])
                        for row in identity_differences
                    ],
                ):
                    return {
                        "classification": "cofactor_pair",
                        "detail": {"identity_differences": identity_differences},
                    }
                return {
                    "classification": "stoichiometry_or_identity_differ",
                    "detail": {"identity_differences": identity_differences},
                }
            if not (
                math.isclose(
                    left["lower_bound"], right["lower_bound"],
                    rel_tol=0.0, abs_tol=1e-9,
                )
                and math.isclose(
                    left["upper_bound"], right["upper_bound"],
                    rel_tol=0.0, abs_tol=1e-9,
                )
            ):
                return {
                    "classification": "direction_or_bounds_differ",
                    "detail": {
                        "left_bounds": [
                            left["lower_bound"], left["upper_bound"]
                        ],
                        "right_bounds": [
                            right["lower_bound"], right["upper_bound"]
                        ],
                    },
                }
            return None

        def threeway_gate_recompute(reaction_id):
            """Recompute the builder's three-way gate from frozen sources."""
            saved_summary = saved_reaction_summary(
                model.reactions.get_by_id(reaction_id)
            )
            template_summary = candidate_source_summary(
                reaction_id, "yeast9_reference"
            )
            pool_summary = candidate_source_summary(
                reaction_id, "frozen_reaction_pool"
            )
            pairs = (
                ("template_pool", template_summary, pool_summary),
                ("model_template", saved_summary, template_summary),
                ("model_pool", saved_summary, pool_summary),
            )
            comparisons = {}
            divergent_pairs = []
            for label, left, right in pairs:
                divergence = classify_summary_divergence(left, right)
                comparisons[label] = divergence
                if divergence is not None:
                    divergent_pairs.append(label)
            return comparisons, divergent_pairs

        def pool_coverage_check(record, reaction_id):
            """Bind the three-way pool_coverage marker to pool membership."""
            in_refmodel = (
                expected_candidate_source(reaction_id) == "yeast9_reference"
            )
            in_pool = reaction_id in pool_ids
            if already_present_policy != "merge_threeway_or_multisource":
                if "pool_coverage" in record:
                    raise RuntimeError(
                        "pool_coverage is recorded only under the three-way "
                        f"policy: {reaction_id}"
                    )
                return
            if in_refmodel and not in_pool:
                if record.get("pool_coverage") != "not_in_pool":
                    raise RuntimeError(
                        "Three-way candidate absent from the pool must record "
                        f"pool_coverage=not_in_pool: {reaction_id}"
                    )
            elif "pool_coverage" in record:
                raise RuntimeError(
                    "pool_coverage must be recorded only for pool-absent "
                    f"reference candidates: {reaction_id}"
                )

        def threeway_merge_gate_check(record, reaction_id):
            """Re-run the three-way gate for a merged/no-merge record."""
            if already_present_policy != "merge_threeway_or_multisource":
                return
            in_refmodel = record["candidate_source"] == "yeast9_reference"
            in_pool = reaction_id in pool_ids
            if in_refmodel and in_pool:
                _, divergent_pairs = threeway_gate_recompute(reaction_id)
                if divergent_pairs:
                    raise RuntimeError(
                        "Candidate merged despite a template/pool divergence: "
                        f"{reaction_id}: {divergent_pairs}"
                    )
                verification = record.get("verification") or {}
                if verification.get("three_way_check") != (
                    "model_template_pool_identical"
                ):
                    raise RuntimeError(
                        "Three-way merge lacks its gate marker: "
                        f"{reaction_id}"
                    )
            pool_coverage_check(record, reaction_id)

        multisource_votes = None
        multisource_loaded = False

        def load_multisource_votes():
            """(Gene, rhea)-pair counts of distinct non-empty anno sources.

            Mirrors the builder: the multiple_score table is keyed per
            (gene, rhea) pair and only rows pointing at the candidate's own
            rhea number count. With ensemble_score 0.4 above the largest
            single-tool weight (anno_B 0.3), passing candidates already carry
            >= 2 pair-level sources, so this recheck is a redundant safeguard
            on the current profile; a missing evidence file blocks all merges.
            """
            nonlocal multisource_votes, multisource_loaded
            if multisource_loaded:
                return multisource_votes
            multisource_loaded = True
            pandas = self._module("pandas")
            evidence_path = working / f"{run.name}multiple_score.xlsx"
            if not evidence_path.is_file():
                return None
            table = pandas.read_excel(evidence_path)
            if not {"reaction", "rhea"} <= set(table.columns):
                multisource_votes = {}
                return multisource_votes
            present_columns = [
                column
                for column in REACTION_POOL_MULTISOURCE_COLUMNS
                if column in table.columns
            ]
            sources = {}
            for _, evidence_row in table.iterrows():
                gene = evidence_row["reaction"]
                rhea = _rhea_number_key(evidence_row["rhea"])
                if not isinstance(gene, str) or not gene.strip() or rhea is None:
                    continue
                pair_sources = sources.setdefault((gene, rhea), set())
                for column in present_columns:
                    value = evidence_row[column]
                    if pandas.notna(value) and str(value).strip():
                        pair_sources.add(column)
            multisource_votes = {
                pair: len(columns) for pair, columns in sources.items()
            }
            return multisource_votes

        draft_model = None

        def load_draft_model():
            nonlocal draft_model
            if draft_model is None:
                draft_model = cobra.io.load_yaml_model(str(draft_path))
            return draft_model

        def expected_candidate_source(reaction_id):
            reference_ids = {
                reaction.id for reaction in load_reference_model().reactions
            }
            if reaction_id in reference_ids:
                return "yeast9_reference"
            if reaction_id in pool_ids:
                return "frozen_reaction_pool"
            return None

        gpr_states = {}
        gpr_tracked_reactions = []

        def pre_merge_gpr_state(reaction_id):
            """GPR rule/genes of a reaction just before the current row.

            Draft-model reactions give the stage-start state; reactions added
            by an earlier row carry that row's candidate rule (tracked in the
            candidate loop below).
            """
            if reaction_id in gpr_states:
                return gpr_states[reaction_id]
            draft = load_draft_model()
            try:
                draft_reaction = draft.reactions.get_by_id(reaction_id)
            except KeyError as error:
                raise RuntimeError(
                    "Already-present candidate was neither in the draft model "
                    f"nor added by an earlier row: {reaction_id}"
                ) from error
            try:
                genes, _ = rule_genes_and_complexity(
                    draft_reaction.gene_reaction_rule
                )
            except (SyntaxError, TypeError, ValueError) as error:
                raise RuntimeError(
                    f"Draft GPR cannot be recomputed: {reaction_id}"
                ) from error
            state = {"rule": draft_reaction.gene_reaction_rule, "genes": genes}
            gpr_states[reaction_id] = state
            return state

        def verification_summary_check(record, summary, source, reaction_id):
            """Recompute the recorded verification summary, names included."""
            # summary is the candidate SOURCE summary; the builder records
            # name mismatches between the saved model reaction's metabolite
            # names and the source reaction's same-ID metabolite names, both
            # compared under the author's metdict name equivalence
            # (lowercase, parentheses stripped) so add-path name swaps like
            # 'H(+)' vs 'H+' are not false mismatches.
            source_names = summary["names"]
            saved_names = saved_reaction_summary(
                model.reactions.get_by_id(reaction_id)
            )["names"]
            name_mismatch_ids = sorted(
                metabolite_id
                for metabolite_id, saved_name in saved_names.items()
                if saved_name
                and source_names.get(metabolite_id)
                and _author_name_key(source_names[metabolite_id])
                != _author_name_key(saved_name)
            )
            verification = record.get("verification")
            if (
                not isinstance(verification, dict)
                or verification.get("candidate_source") != source
                or bool(verification.get("name_mismatch"))
                != bool(name_mismatch_ids)
                or verification.get("name_mismatch_metabolites")
                != name_mismatch_ids
            ):
                raise RuntimeError(
                    "Candidate verification summary does not recompute from "
                    f"source and saved model: {reaction_id}"
                )

        def multisource_skipped_check(record, candidate_genes, row, reaction_id):
            """Recompute the pair-level skipped-gene list for a record."""
            skipped = record.get("skipped_insufficient_evidence")
            if not isinstance(skipped, list):
                raise RuntimeError(
                    "Multisource candidate lacks its skipped-gene list: "
                    f"{reaction_id}"
                )
            votes = load_multisource_votes()
            rhea = _rhea_number_key(row["rhea"])
            expected_skipped = sorted(
                gene
                for gene in candidate_genes
                if votes is None
                or votes.get((gene, rhea), 0)
                < REACTION_POOL_MULTISOURCE_MIN_SOURCES
            )
            if sorted(skipped) != expected_skipped:
                raise RuntimeError(
                    "Multisource skipped-gene list does not match the "
                    f"annotation evidence: {reaction_id}"
                )
            return set(skipped)

        def merged_record_check(record, row, reaction_id):
            """Recompute an already_present_gpr_merged record from sources."""
            source = record.get("candidate_source")
            if source != expected_candidate_source(reaction_id) or source is None:
                raise RuntimeError(
                    "Merged candidate source does not match the reference/"
                    f"pool resolution: {reaction_id}"
                )
            summary = candidate_source_summary(reaction_id, source)
            try:
                saved_reaction = model.reactions.get_by_id(reaction_id)
            except KeyError as error:
                raise RuntimeError(
                    "Merged candidate is absent from the saved model: "
                    f"{reaction_id}"
                ) from error
            if (
                classify_existing_vs_source(
                    saved_reaction_summary(saved_reaction), summary
                )
                is not None
            ):
                raise RuntimeError(
                    "Merged candidate no longer matches its source chemistry: "
                    f"{reaction_id}"
                )
            candidate_value = row.get("reaction")
            candidate_rule = (
                candidate_value if isinstance(candidate_value, str) else ""
            )
            try:
                candidate_genes, candidate_complex = rule_genes_and_complexity(
                    candidate_rule
                )
            except (SyntaxError, TypeError, ValueError) as error:
                raise RuntimeError(
                    "Merged candidate GPR cannot be recomputed: "
                    f"{reaction_id}"
                ) from error
            state = pre_merge_gpr_state(reaction_id)
            _, state_complex = rule_genes_and_complexity(state["rule"])
            if candidate_complex or state_complex:
                raise RuntimeError(
                    f"Merged candidate involves a complex GPR: {reaction_id}"
                )
            merged_genes = record.get("merged_genes")
            if (
                not isinstance(merged_genes, list)
                or any(not isinstance(gene, str) for gene in merged_genes)
                or merged_genes != sorted(merged_genes)
            ):
                raise RuntimeError(
                    "Merged candidate has an invalid merged_genes list: "
                    f"{reaction_id}"
                )
            if already_present_policy in {
                "merge_verified_identical_or_multisource",
                "merge_threeway_or_multisource",
            }:
                skipped_set = multisource_skipped_check(
                    record, candidate_genes, row, reaction_id
                )
            else:
                if "skipped_insufficient_evidence" in record:
                    raise RuntimeError(
                        "Non-multisource merged candidate records skipped "
                        f"genes: {reaction_id}"
                    )
                skipped_set = set()
            # The merge is bound exactly to the pre-merge GPR: nothing may be
            # dropped, added, or under/over-reported relative to it.
            eligible = candidate_genes - skipped_set
            state_genes = state["genes"]
            if merged_genes != sorted(eligible - state_genes):
                raise RuntimeError(
                    "Merged genes differ from eligible minus pre-existing "
                    f"genes: {reaction_id}"
                )
            if record.get("merge_kind") != (
                "into_empty_gpr" if not state_genes else "into_existing_gpr"
            ):
                raise RuntimeError(
                    f"Merge kind does not match the pre-merge GPR: {reaction_id}"
                )
            verification_summary_check(record, summary, source, reaction_id)
            threeway_merge_gate_check(record, reaction_id)
            gpr_states[reaction_id] = {
                "rule": " or ".join(sorted(state_genes | set(merged_genes))),
                "genes": state_genes | set(merged_genes),
            }
            gpr_tracked_reactions.append(reaction_id)
            return {
                "input_row": record["input_row"],
                "reaction_id": reaction_id,
                "merge_kind": record["merge_kind"],
                "merged_genes": merged_genes,
            }

        def no_merge_record_check(record, row, reaction_id):
            """Recompute an already_present_no_merge record from sources."""
            if "merged_genes" in record or "merge_kind" in record:
                raise RuntimeError(
                    f"No-merge candidate carries merge fields: {reaction_id}"
                )
            if (
                record.get("no_merge_reason")
                != "all_candidate_genes_below_min_sources"
            ):
                raise RuntimeError(
                    f"No-merge candidate has an unknown reason: {reaction_id}"
                )
            source = record.get("candidate_source")
            if source != expected_candidate_source(reaction_id) or source is None:
                raise RuntimeError(
                    "No-merge candidate source does not match the reference/"
                    f"pool resolution: {reaction_id}"
                )
            summary = candidate_source_summary(reaction_id, source)
            try:
                saved_reaction = model.reactions.get_by_id(reaction_id)
            except KeyError as error:
                raise RuntimeError(
                    "No-merge candidate is absent from the saved model: "
                    f"{reaction_id}"
                ) from error
            if (
                classify_existing_vs_source(
                    saved_reaction_summary(saved_reaction), summary
                )
                is not None
            ):
                raise RuntimeError(
                    "No-merge candidate no longer matches its source "
                    f"chemistry: {reaction_id}"
                )
            candidate_value = row.get("reaction")
            candidate_rule = (
                candidate_value if isinstance(candidate_value, str) else ""
            )
            try:
                candidate_genes, candidate_complex = rule_genes_and_complexity(
                    candidate_rule
                )
            except (SyntaxError, TypeError, ValueError) as error:
                raise RuntimeError(
                    "No-merge candidate GPR cannot be recomputed: "
                    f"{reaction_id}"
                ) from error
            state = pre_merge_gpr_state(reaction_id)
            _, state_complex = rule_genes_and_complexity(state["rule"])
            if candidate_complex or state_complex:
                raise RuntimeError(
                    f"No-merge candidate involves a complex GPR: {reaction_id}"
                )
            skipped_set = multisource_skipped_check(
                record, candidate_genes, row, reaction_id
            )
            if candidate_genes - skipped_set:
                raise RuntimeError(
                    "No-merge candidate actually had eligible genes: "
                    f"{reaction_id}"
                )
            verification_summary_check(record, summary, source, reaction_id)
            threeway_merge_gate_check(record, reaction_id)
            gpr_tracked_reactions.append(reaction_id)
            return {
                "input_row": record["input_row"],
                "reaction_id": reaction_id,
                "no_merge_reason": record["no_merge_reason"],
            }

        def conflict_record_check(record, row, reaction_id):
            """Recompute an already_present_conflict record from sources."""
            source = expected_candidate_source(reaction_id)
            try:
                saved_reaction = model.reactions.get_by_id(reaction_id)
            except KeyError as error:
                raise RuntimeError(
                    "Conflicted candidate is absent from the saved model: "
                    f"{reaction_id}"
                ) from error
            conflict_kind = record.get("conflict_kind")
            if conflict_kind == "template_pool_divergence":
                if (
                    source != "yeast9_reference"
                    or reaction_id not in pool_ids
                    or record.get("candidate_source") != source
                ):
                    raise RuntimeError(
                        "Template/pool divergence requires the candidate in "
                        f"both reference and pool: {reaction_id}"
                    )
                comparisons, divergent_pairs = threeway_gate_recompute(
                    reaction_id
                )
                if not divergent_pairs:
                    raise RuntimeError(
                        "Divergence record but all three versions are "
                        f"identical: {reaction_id}"
                    )
                expected_detail = {
                    "classification": comparisons[divergent_pairs[0]][
                        "classification"
                    ],
                    "divergent_pairs": divergent_pairs,
                    "pairs": comparisons,
                }
                if record.get("conflict_detail") != expected_detail:
                    raise RuntimeError(
                        "Template/pool divergence detail does not recompute "
                        f"from the frozen sources: {reaction_id}"
                    )
            elif source is None:
                if (
                    record.get("candidate_source") is not None
                    or conflict_kind != "unverifiable_identity"
                    or record.get("conflict_detail")
                    != {"reason": "candidate_source_absent"}
                ):
                    raise RuntimeError(
                        "Source-absent conflict record does not match the "
                        f"reference/pool resolution: {reaction_id}"
                    )
            else:
                if record.get("candidate_source") != source:
                    raise RuntimeError(
                        "Conflicted candidate source does not match the "
                        f"reference/pool resolution: {reaction_id}"
                    )
                summary = candidate_source_summary(reaction_id, source)
                expected_kind = classify_existing_vs_source(
                    saved_reaction_summary(saved_reaction), summary
                )
                if expected_kind is None:
                    candidate_value = row.get("reaction")
                    candidate_rule = (
                        candidate_value
                        if isinstance(candidate_value, str)
                        else ""
                    )
                    try:
                        _, candidate_complex = rule_genes_and_complexity(
                            candidate_rule
                        )
                        _, saved_complex = rule_genes_and_complexity(
                            saved_reaction.gene_reaction_rule
                        )
                    except (SyntaxError, TypeError, ValueError):
                        expected_kind = "unverifiable_identity"
                    else:
                        if candidate_complex or saved_complex:
                            expected_kind = "complex_gpr_not_merged"
                if expected_kind is None:
                    raise RuntimeError(
                        "Conflicted candidate is actually identical to its "
                        f"source: {reaction_id}"
                    )
                if conflict_kind != expected_kind:
                    raise RuntimeError(
                        "Conflict kind does not recompute from source and "
                        f"saved model: {reaction_id}: "
                        f"{conflict_kind!r} != {expected_kind!r}"
                    )
            if not isinstance(record.get("conflict_detail"), dict):
                raise RuntimeError(
                    f"Conflicted candidate lacks conflict detail: {reaction_id}"
                )
            pool_coverage_check(record, reaction_id)
            conflict_row = {
                "input_row": record["input_row"],
                "reaction_id": reaction_id,
                "conflict_kind": conflict_kind,
            }
            if conflict_kind == "template_pool_divergence":
                conflict_row["classification"] = record["conflict_detail"][
                    "classification"
                ]
            return conflict_row

        def addition_record_check(record, reaction_id):
            """Recompute a suffix/remap identity-conflict resolution."""
            source = record.get("candidate_source")
            if (
                source not in {"frozen_reaction_pool", "yeast9_reference"}
                or source != expected_candidate_source(reaction_id)
            ):
                raise RuntimeError(
                    "Resolved-conflict candidate source does not match the "
                    f"reference/pool resolution: {reaction_id}"
                )
            summary = candidate_source_summary(reaction_id, source)
            try:
                saved_reaction = model.reactions.get_by_id(reaction_id)
            except KeyError as error:
                raise RuntimeError(
                    "Resolved-conflict candidate is absent from the saved "
                    f"model: {reaction_id}"
                ) from error
            saved_coefficients = {
                metabolite.id: float(coefficient)
                for metabolite, coefficient in saved_reaction.metabolites.items()
            }

            def model_identity(metabolite_id):
                try:
                    metabolite = model.metabolites.get_by_id(metabolite_id)
                except KeyError as error:
                    raise RuntimeError(
                        f"Saved model lacks metabolite {metabolite_id} for "
                        f"resolved candidate {reaction_id}"
                    ) from error
                return {
                    "formula": metabolite.formula,
                    "charge": metabolite.charge,
                    "compartment": metabolite.compartment,
                }

            def check_common(entry):
                old_id = entry.get("metabolite_id")
                if not isinstance(old_id, str) or old_id not in summary[
                    "coefficients"
                ]:
                    raise RuntimeError(
                        "Resolved candidate references a metabolite absent "
                        f"from its source: {reaction_id}: {old_id!r}"
                    )
                source_identity = summary["identities"][old_id]
                if any(
                    value is None or value == ""
                    for value in source_identity.values()
                ):
                    raise RuntimeError(
                        "Resolved candidate has an incomplete source "
                        f"identity: {reaction_id}: {old_id}"
                    )
                existing_identity = model_identity(old_id)
                if not any(
                    value is None or value == ""
                    for value in existing_identity.values()
                ) and existing_identity == source_identity:
                    raise RuntimeError(
                        "Recorded conflict resolution had no identity "
                        f"conflict: {reaction_id}: {old_id}"
                    )
                if {
                    field: entry.get(field)
                    for field in ("formula", "charge", "compartment")
                } != source_identity:
                    raise RuntimeError(
                        "Resolution record attributes differ from the "
                        f"candidate source: {reaction_id}: {old_id}"
                    )
                if old_id in saved_coefficients:
                    raise RuntimeError(
                        "Conflicted metabolite is still referenced by the "
                        f"saved reaction: {reaction_id}: {old_id}"
                    )
                return old_id, source_identity

            addition_kind = record.get("addition_kind")
            if addition_kind == "suffix_metabolite":
                entries = record.get("suffix_metabolites")
                if not isinstance(entries, list) or not entries:
                    raise RuntimeError(
                        "Suffix addition lacks suffix_metabolites: "
                        f"{reaction_id}"
                    )
                for entry in entries:
                    old_id, source_identity = check_common(entry)
                    new_id = entry.get("new_metabolite_id")
                    if new_id != f"{old_id}__pool":
                        raise RuntimeError(
                            "Suffix metabolite ID does not follow the "
                            f"__pool convention: {reaction_id}: {new_id!r}"
                        )
                    if model_identity(new_id) != source_identity:
                        raise RuntimeError(
                            "Suffix metabolite identity differs from the "
                            f"candidate source: {reaction_id}: {new_id}"
                        )
                    observed = saved_coefficients.get(new_id)
                    if observed is None or not math.isclose(
                        observed,
                        summary["coefficients"][old_id],
                        rel_tol=0.0,
                        abs_tol=1e-9,
                    ):
                        raise RuntimeError(
                            "Saved reaction does not reference the suffix "
                            f"metabolite with the source coefficient: "
                            f"{reaction_id}: {new_id}"
                        )
            elif addition_kind == "remapped_metabolite":
                entries = record.get("remapped_metabolites")
                if not isinstance(entries, list) or not entries:
                    raise RuntimeError(
                        "Remap addition lacks remapped_metabolites: "
                        f"{reaction_id}"
                    )
                draft = load_draft_model()
                for entry in entries:
                    old_id, source_identity = check_common(entry)
                    target_id = entry.get("remapped_to")
                    if not isinstance(target_id, str) or target_id == old_id:
                        raise RuntimeError(
                            "Invalid remap target: "
                            f"{reaction_id}: {target_id!r}"
                        )
                    try:
                        target = model.metabolites.get_by_id(target_id)
                    except KeyError as error:
                        raise RuntimeError(
                            f"Remap target absent from the saved model: "
                            f"{reaction_id}: {target_id}"
                        ) from error
                    target_identity = {
                        "formula": target.formula,
                        "charge": target.charge,
                        "compartment": target.compartment,
                    }
                    if target_identity != source_identity:
                        raise RuntimeError(
                            "Remap target identity differs from the "
                            f"candidate source: {reaction_id}: {target_id}"
                        )
                    source_name = summary["names"].get(old_id)
                    if not normalized_names(
                        source_name
                    ) & normalized_names(target.name):
                        raise RuntimeError(
                            "Remap target name does not match the candidate "
                            f"source name: {reaction_id}: {old_id} -> "
                            f"{target_id}"
                        )
                    observed = saved_coefficients.get(target_id)
                    if observed is None or not math.isclose(
                        observed,
                        summary["coefficients"][old_id],
                        rel_tol=0.0,
                        abs_tol=1e-9,
                    ):
                        raise RuntimeError(
                            "Saved reaction does not reference the remap "
                            "target with the source coefficient: "
                            f"{reaction_id}: {target_id}"
                        )
                    # The runtime matched against the stage model; the draft
                    # model is the fixed pre-stage baseline used to rederive
                    # match uniqueness here.
                    draft_matches = [
                        metabolite.id
                        for metabolite in draft.metabolites
                        if metabolite.id != old_id
                        and normalized_names(metabolite.name)
                        & normalized_names(source_name)
                        and {
                            "formula": metabolite.formula,
                            "charge": metabolite.charge,
                            "compartment": metabolite.compartment,
                        }
                        == source_identity
                    ]
                    if draft_matches != [target_id]:
                        raise RuntimeError(
                            "Remap target is not the unique draft-model "
                            f"identity match: {reaction_id}: {draft_matches}"
                        )
            else:
                raise RuntimeError(
                    f"Unknown addition kind: {addition_kind!r}"
                )
            # The recorded blocked flag must recompute from the saved model:
            # the candidate's genuinely new metabolites (created __pool
            # suffixes plus previously unmatched ones) each appearing in
            # exactly one reaction.
            saved_use_count = {}
            for saved_stage_reaction in model.reactions:
                for metabolite in saved_stage_reaction.metabolites:
                    saved_use_count[metabolite.id] = (
                        saved_use_count.get(metabolite.id, 0) + 1
                    )
            new_ids = [
                entry["new_metabolite_id"]
                for entry in record.get("suffix_metabolites", [])
                if not entry.get("reused_existing_pool_metabolite")
            ]
            new_ids += record.get("unmatched_metabolites", [])
            expected_blocked = bool(new_ids) and all(
                saved_use_count.get(metabolite_id, 0) == 1
                for metabolite_id in new_ids
            )
            recorded_blocked = record.get("blocked_by_construction")
            if (
                not isinstance(recorded_blocked, bool)
                or recorded_blocked != expected_blocked
            ):
                raise RuntimeError(
                    "Recorded blocked_by_construction does not recompute from "
                    f"the saved model: {reaction_id}"
                )
            return {
                "input_row": record["input_row"],
                "reaction_id": reaction_id,
                "addition_kind": addition_kind,
                "blocked_by_construction": expected_blocked,
            }

        merged_rows = []
        no_merge_rows = []
        conflict_rows = []
        addition_rows = []
        for (index, row), record in zip(score_rows.iterrows(), candidates):
            expected_row = int(index) if isinstance(index, int) else str(index)
            reaction_id = "rhea_" + str(row["rhea"]) + "_c"
            if (
                record.get("input_row") != expected_row
                or record.get("reaction_id") != reaction_id
            ):
                raise RuntimeError(
                    "Formal reaction-pool audit candidate differs from "
                    f"input row {expected_row}"
                )
            status = record["status"]
            if status == "missing_from_reaction_pool":
                if reaction_id in pool_ids:
                    raise RuntimeError(
                        "Reaction-pool audit falsely marks a present "
                        f"reaction as missing: {reaction_id}"
                    )
                missing.append({
                    "input_row": expected_row,
                    "reaction_id": reaction_id,
                })
            elif (
                status
                not in {
                    "already_present",
                    "already_present_conflict",
                    "already_present_gpr_merged",
                    "already_present_no_merge",
                }
                and reaction_id not in pool_ids
            ):
                raise RuntimeError(
                    "Reaction-pool audit assigns a pool-dependent status "
                    f"to absent reaction {reaction_id}: {status}"
                )
            if status == "rejected_postmerge_metabolite_identity_conflict":
                source = record.get("candidate_source")
                if source == "frozen_reaction_pool":
                    source_reaction = pool_reactions_by_id[reaction_id]
                    source_metabolites = source_reaction.get("metabolites")
                    if not isinstance(source_metabolites, dict):
                        raise RuntimeError(
                            "Frozen reaction-pool candidate has no stoichiometry"
                        )
                    identities = {
                        metabolite_id: pool_metabolites_by_id[metabolite_id]
                        for metabolite_id in source_metabolites
                    }
                elif source == "yeast9_reference":
                    if reference_model is None:
                        reference_model = cobra.io.read_sbml_model(
                            str(reference_path)
                        )
                    try:
                        reference_reaction = reference_model.reactions.get_by_id(
                            reaction_id
                        )
                    except KeyError as error:
                        raise RuntimeError(
                            "Conflict candidate is absent from reference model"
                        ) from error
                    identities = {
                        metabolite.id: {
                            "formula": metabolite.formula,
                            "charge": metabolite.charge,
                            "compartment": metabolite.compartment,
                        }
                        for metabolite in reference_reaction.metabolites
                    }
                else:
                    raise RuntimeError(
                        "Identity-conflict candidate lacks a known source"
                    )
                expected_conflicts = []
                for metabolite_id, source_metabolite in identities.items():
                    try:
                        model_metabolite = model.metabolites.get_by_id(
                            metabolite_id
                        )
                    except KeyError:
                        continue
                    candidate_identity = {
                        field: source_metabolite.get(field)
                        for field in ("formula", "charge", "compartment")
                    }
                    model_identity = {
                        "formula": model_metabolite.formula,
                        "charge": model_metabolite.charge,
                        "compartment": model_metabolite.compartment,
                    }
                    incomplete = any(
                        value is None or value == ""
                        for identity in (candidate_identity, model_identity)
                        for value in identity.values()
                    )
                    if incomplete or candidate_identity != model_identity:
                        expected_conflicts.append({
                            "metabolite_id": metabolite_id,
                            "conflict_kind": (
                                "unverifiable_identity"
                                if incomplete else "identity_mismatch"
                            ),
                            "candidate_identity": candidate_identity,
                            "model_identity": model_identity,
                        })
                if (
                    record.get("rejection_reason")
                    != "postmerge_metabolite_identity_conflict"
                    or not expected_conflicts
                    or sorted(record.get("identity_conflicts", []), key=lambda x: x["metabolite_id"])
                    != sorted(expected_conflicts, key=lambda x: x["metabolite_id"])
                ):
                    raise RuntimeError(
                        "Reaction-pool identity-conflict record does not "
                        f"match source and saved-model chemistry: {reaction_id}"
                    )
                identity_rejections.append({
                    "input_row": expected_row,
                    "reaction_id": reaction_id,
                    "candidate_source": source,
                    "rejection_reason": record["rejection_reason"],
                    "identity_conflicts": expected_conflicts,
                })
            elif "identity_conflicts" in record or "rejection_reason" in record:
                raise RuntimeError(
                    "Reaction-pool non-conflict candidate has conflict metadata"
                )
            if status == "already_present_gpr_merged":
                merged_rows.append(merged_record_check(record, row, reaction_id))
            elif status == "already_present_no_merge":
                no_merge_rows.append(
                    no_merge_record_check(record, row, reaction_id)
                )
            elif status == "already_present_conflict":
                conflict_rows.append(
                    conflict_record_check(record, row, reaction_id)
                )
            if status == "added" and "addition_kind" in record:
                addition_rows.append(addition_record_check(record, reaction_id))
            if status == "added":
                # The builder sets an added reaction's GPR from its row; track
                # it so a later duplicate row's merge binds to this state.
                added_value = row.get("reaction")
                added_rule = (
                    added_value if isinstance(added_value, str) else ""
                )
                try:
                    added_genes, _ = rule_genes_and_complexity(added_rule)
                except (SyntaxError, TypeError, ValueError):
                    pass
                else:
                    gpr_states[reaction_id] = {
                        "rule": added_rule,
                        "genes": added_genes,
                    }
            if status in {
                "added", "already_present", "already_present_conflict",
                "already_present_gpr_merged", "already_present_no_merge",
                "missing_from_reaction_pool",
                "rejected_postmerge_metabolite_identity_conflict",
            }:
                needs_model.append((reaction_id, status))

        # A pool absence never licenses an imbalanced addition. Re-read the
        # saved stage model and check that each claimed addition is balanced.
        model_ids = {reaction.id for reaction in model.reactions}
        for reaction_id in gpr_tracked_reactions:
            # Every recorded merge/no-merge must reproduce the saved GPR
            # exactly from the pre-merge state: no dropped original genes, no
            # foreign genes, no unrecorded additions.
            try:
                tracked = model.reactions.get_by_id(reaction_id)
            except KeyError as error:
                raise RuntimeError(
                    "Merged/no-merge candidate is absent from the saved "
                    f"model: {reaction_id}"
                ) from error
            try:
                saved_genes, _ = rule_genes_and_complexity(
                    tracked.gene_reaction_rule
                )
            except (SyntaxError, TypeError, ValueError) as error:
                raise RuntimeError(
                    f"Saved GPR cannot be recomputed: {reaction_id}"
                ) from error
            if saved_genes != gpr_states[reaction_id]["genes"]:
                raise RuntimeError(
                    "Saved GPR does not equal the pre-merge state plus the "
                    f"recorded merged genes: {reaction_id}"
                )
        for reaction_id, status in needs_model:
            if status in {
                "missing_from_reaction_pool",
                "rejected_postmerge_metabolite_identity_conflict",
            }:
                if reaction_id in model_ids:
                    raise RuntimeError(
                        "Reaction marked absent or identity-conflicted was added to "
                        f"the model: {reaction_id}"
                    )
                continue
            if reaction_id not in model_ids:
                raise RuntimeError(
                    "Reaction-pool audit claims a model reaction absent "
                    f"from the saved model: {reaction_id}"
                )
            if status == "added":
                try:
                    imbalance = model.reactions.get_by_id(
                        reaction_id
                    ).check_mass_balance()
                except (TypeError, ValueError) as error:
                    raise RuntimeError(
                        "Added reaction mass/charge balance could not be "
                        f"checked: {reaction_id}"
                    ) from error
                if imbalance:
                    raise RuntimeError(
                        "Reaction-pool audit added an imbalanced reaction: "
                        f"{reaction_id}: {imbalance}"
                    )

        receipt = {
            "policy_version": REACTION_POOL_AUDIT_POLICY_VERSION,
            "author_behavior": "missing_pool_reaction_skipped",
            "reaction_pool_policies": {
                "already_present_policy": already_present_policy,
                "conflict_policy": conflict_policy,
            },
            "reaction_pool": pool_record,
            "candidate_audit": self._file_record(audit_path),
            "candidate_coverage": self._file_record(
                working / "reaction_pool_coverage.json"
            ),
            "candidate_spreadsheet": self._file_record(score_path),
            "draft_model": self._file_record(draft_path),
            "reference_model": reference_record,
            "saved_model": self._file_record(model_path),
            "input_candidate_count": len(candidates),
            "status_counts": audit["status_counts"],
            "missing_from_reaction_pool": missing,
            "postmerge_identity_conflicts": identity_rejections,
            "postmerge_identity_conflict_count": len(identity_rejections),
            "already_present_gpr_merged": merged_rows,
            "already_present_no_merge": no_merge_rows,
            "already_present_conflicts": conflict_rows,
            "identity_conflict_additions": addition_rows,
            "excluded_released_code_count": 0,
        }
        (working / "reaction_pool_validation_v3.json").write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @staticmethod
    def _file_record(path: Path):
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return {
            "bytes": path.stat().st_size,
            "sha256": digest.hexdigest(),
        }

    @classmethod
    def _runtime_file_record(cls, path: Path):
        path = Path(path)
        metadata = os.stat(path, follow_symlinks=False)
        return {
            **cls._file_record(path),
            "device": metadata.st_dev,
            "inode": metadata.st_ino,
            "mtime_ns": metadata.st_mtime_ns,
            "ctime_ns": metadata.st_ctime_ns,
        }

    @staticmethod
    def _same_file_identity(left, right):
        if left is None or right is None:
            return left is right
        return (
            left.get("bytes") == right.get("bytes")
            and left.get("sha256") == right.get("sha256")
        )

    @classmethod
    def _same_gapfill_input_identity(cls, left, right):
        if left is None or right is None:
            return left is right
        if set(left) != set(right):
            return False
        return all(
            cls._same_file_identity(left[key], right[key])
            for key in left
        )

    def _write_manifest(
        self,
        run: PipelineRun,
        completed: List[str],
        stage_records,
    ) -> None:
        working = self.repo_root / "working" / run.name
        working.mkdir(parents=True, exist_ok=True)
        source_fasta = self._resolve_input(run.fasta)
        self._assert_formal_run_identity(run)
        source_code = (
            self._run_source_identity
            if self.verify_outputs and self._is_exact_paper_profile(run)
            else self._source_code_identity()
        )
        runtime_environment = (
            self._run_runtime_environment
            if self.verify_outputs and self._is_exact_paper_profile(run)
            else self._runtime_environment()
        )
        payload = {
            "schema_version": self.MANIFEST_SCHEMA_VERSION,
            "updated_at": datetime.now(timezone.utc).isoformat(),
            "verification_mode": self._verification_mode(run, completed),
            "gapfill_method_boundary": (
                self._manifest_gapfill_method_boundary(run, completed)
            ),
            "source_code": source_code,
            "runtime_environment": runtime_environment,
            "paper_method_evidence": self._paper_method_evidence(run),
            "name": run.name,
            "refname": run.refname,
            "mode": run.mode,
            "profile": asdict(run.profile),
            "solver": run.solver,
            "resolved_solver": self.solver_record,
            "target_input_contract": self._target_input_contract_record(run),
            "gapfill_coordinate_contract": (
                self._gapfill_coordinate_contract_record(run)
            ),
            "blast_threads": run.blast_threads,
            "esm_token_budget": run.esm_token_budget,
            "clean": run.clean,
            "deepec": run.deepec,
            "deepec_python": (
                str(self._resolve_input(run.deepec_python))
                if run.deepec_python is not None
                else None
            ),
            "plm_annotation": run.plm_annotation,
            "growth_medium": run.growth_medium,
            "essential_reactions": sorted(run.essential_reactions or []),
            "reaction_pool_policies": self._reaction_pool_policy_record(run),
            "gapfill_tie_break": run.gapfill_tie_break,
            "gapfill_unmapped_hit_policy": run.gapfill_unmapped_hit_policy,
            "completed_stages": completed,
            "source_fasta": {
                "observed_path": str(source_fasta),
                **self._file_record(source_fasta),
            },
            "gapfill_sequence_input": self._gapfill_sequence_record(run),
            "stages": stage_records,
        }
        if self._checkpoint_import_record is not None:
            payload["imported_from_checkpoint"] = self._checkpoint_import_record
        manifest = working / "alphagem_run_manifest.json"
        temporary = manifest.with_name(manifest.name + ".tmp")
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, manifest)

    def _verification_mode(
        self,
        run: PipelineRun,
        completed: Optional[List[str]] = None,
    ) -> str:
        if not self.verify_outputs:
            return "unverified"
        if self._is_exact_paper_profile(run):
            if "gapfill" not in (completed or []):
                return (
                    "formal_paper_declared_method_contract_in_progress"
                )
            boundary = self._manifest_gapfill_method_boundary(
                run,
                completed or [],
            )
            if boundary["status"] == "declared_operationalization":
                return (
                    "formal_paper_declared_method_contract_with_"
                    "gapfill_operationalization"
                )
            return "formal_paper_declared_method_contract"
        return "verified_control"

    @staticmethod
    def _paper_method_evidence(run: PipelineRun):
        if not PipelineRunner._is_exact_paper_profile(run):
            return None
        return PAPER_METHOD_EVIDENCE_CONTRACT

    @staticmethod
    def _reaction_pool_policy_record(run: PipelineRun):
        return {
            "already_present_policy": run.reaction_pool_already_present_policy,
            "conflict_policy": run.reaction_pool_conflict_policy,
        }

    def _manifest_gapfill_method_boundary(
        self,
        run: PipelineRun,
        completed: List[str],
    ):
        if not self.verify_outputs or not self._is_exact_paper_profile(run):
            return None
        if "gapfill" not in completed:
            return {"status": "not_reached"}
        audit_path = (
            self.repo_root
            / "working"
            / run.name
            / "gapfill_audit.json"
        )
        try:
            audit = json.loads(audit_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise RuntimeError(
                "A completed formal gapfill stage requires a valid audit: "
                f"{audit_path}"
            ) from error
        return self._gapfill_method_boundary_record(audit)

    @staticmethod
    def _gapfill_method_boundary_record(audit):
        def format_validation_proves_growth(minimum):
            validation = audit.get("format_validation")
            if (
                not isinstance(validation, dict)
                or validation.get("formats_equivalent") is not True
            ):
                return False
            for output_format in ("xml", "json"):
                record = validation.get(output_format)
                objective = (
                    record.get("objective_value")
                    if isinstance(record, dict)
                    else None
                )
                if (
                    not isinstance(record, dict)
                    or record.get("status") != "optimal"
                    or not isinstance(objective, (int, float))
                    or isinstance(objective, bool)
                    or not math.isfinite(float(objective))
                    or float(objective) <= minimum
                ):
                    return False
            return True

        rule = audit.get("growth_retention_rule")
        if rule is None:
            forbidden_decision_fields = {
                "removed_unsupported_reactions",
                "retained_algorithmic_reactions",
                "minimum_retained_reaction_count",
                "removal_trials",
            }
            objective = audit.get("objective_value")
            valid_early_exit = (
                audit.get("gapfill_exit_reason")
                == "initial_model_already_grows"
                and audit.get("solver_status") == "optimal"
                and isinstance(objective, (int, float))
                and not isinstance(objective, bool)
                and math.isfinite(float(objective))
                and objective > 0.001
                and audit.get("candidate_reactions") == []
                and audit.get("net_added_reactions") == []
                and audit.get("sequence_search_performed") is False
                and format_validation_proves_growth(0.001)
                and not (forbidden_decision_fields & set(audit))
            )
            if not valid_early_exit:
                raise RuntimeError(
                    "Formal gapfill audit without a growth-retention rule "
                    "does not prove the initial-growth early exit"
                )
            return {
                "status": "not_invoked",
                "exit_reason": "initial_model_already_grows",
                "initial_objective_value": objective,
            }
        if audit.get("gapfill_exit_reason") != "gapfill_completed":
            raise RuntimeError(
                "Formal gapfill audit with a growth-retention rule lacks the "
                "completed gapfill exit reason"
            )
        required = {
            "source",
            "objective_fraction",
            "objective_cap",
            "absolute_threshold",
            "selection",
        }
        if not isinstance(rule, dict) or not required.issubset(rule):
            raise RuntimeError(
                "Formal gapfill audit has an incomplete growth-retention rule"
            )
        removed = audit.get("removed_unsupported_reactions")
        retained = audit.get("retained_algorithmic_reactions")
        trials = audit.get("removal_trials")
        minimum_count = audit.get("minimum_retained_reaction_count")
        if (
            not isinstance(removed, list)
            or not isinstance(retained, list)
            or not isinstance(trials, dict)
            or len(set(removed)) != len(removed)
            or len(set(retained)) != len(retained)
            or set(removed) & set(retained)
            or set(removed) | set(retained) != set(trials)
            or minimum_count != len(retained)
        ):
            raise RuntimeError(
                "Formal gapfill audit violates unsupported-reaction "
                "decision conservation"
            )
        objective = audit.get("objective_value")
        if (
            audit.get("solver_status") != "optimal"
            or not isinstance(objective, (int, float))
            or isinstance(objective, bool)
            or not math.isfinite(float(objective))
            or float(objective) <= 1e-5
        ):
            raise RuntimeError(
                "Formal completed gapfill audit lacks an optimal finite final "
                "growth result"
            )
        if not format_validation_proves_growth(1e-5):
            raise RuntimeError(
                "Formal completed gapfill audit does not prove equivalent "
                "growing XML and JSON outputs"
            )
        pruning_objective = audit.get("final_growth_objective")
        if (
            audit.get("final_growth_status") != "optimal"
            or not isinstance(pruning_objective, (int, float))
            or isinstance(pruning_objective, bool)
            or not math.isfinite(float(pruning_objective))
            or float(pruning_objective) <= float(rule["absolute_threshold"])
        ):
            raise RuntimeError(
                "Formal completed gapfill audit lacks an optimal finite "
                "post-pruning growth result"
            )
        essential_reactions = audit.get("essential_reactions")
        essential_records = audit.get("essential_reaction_enforcement")
        if (
            not isinstance(essential_reactions, list)
            or not isinstance(essential_records, list)
            or any(
                not isinstance(reaction_id, str)
                for reaction_id in essential_reactions
            )
            or len(set(essential_reactions)) != len(essential_reactions)
            or any(
                not isinstance(record, dict)
                or not isinstance(record.get("reaction_id"), str)
                for record in essential_records
            )
            or sorted(set(essential_reactions))
            != sorted(
                record.get("reaction_id")
                for record in essential_records
            )
            or any(
                not isinstance(record, dict)
                or record.get("status") not in {"already_present", "added"}
                or record.get("final_present") is not True
                for record in essential_records
            )
        ):
            raise RuntimeError(
                "Formal completed gapfill audit does not prove essential "
                "reaction restoration"
            )
        candidates = sorted(set(removed) | set(retained))
        status = (
            "declared_operationalization"
            if candidates
            else "invoked_without_unsupported_candidates"
        )
        return {
            "status": status,
            "unsupported_candidate_count": len(candidates),
            "removed_count": len(removed),
            "retained_count": len(retained),
            "source": rule["source"],
            "objective_fraction": rule["objective_fraction"],
            "objective_cap": rule["objective_cap"],
            "absolute_threshold": rule["absolute_threshold"],
            "selection": rule["selection"],
        }

    def _initialize_formal_run_identity(self, run: PipelineRun) -> None:
        if not (
            self.verify_outputs and self._is_exact_paper_profile(run)
        ):
            return
        self._run_source_identity = self._source_code_identity()
        self._run_runtime_environment = self._runtime_environment()

    def _assert_formal_run_identity(self, run: PipelineRun) -> None:
        if not (
            self.verify_outputs and self._is_exact_paper_profile(run)
        ):
            return
        if (
            self._run_source_identity is None
            or self._run_runtime_environment is None
        ):
            raise RuntimeError(
                "Formal paper run identity was not initialized by validation"
            )
        current_source = self._source_code_identity()
        if current_source != self._run_source_identity:
            raise RuntimeError(
                "Formal paper source identity changed during the run: "
                + json.dumps(
                    {
                        "expected": self._run_source_identity,
                        "observed": current_source,
                    },
                    sort_keys=True,
                )
            )
        current_runtime = self._runtime_environment()
        if current_runtime != self._run_runtime_environment:
            raise RuntimeError(
                "Formal paper runtime environment changed during the run: "
                + json.dumps(
                    {
                        "expected": self._run_runtime_environment,
                        "observed": current_runtime,
                    },
                    sort_keys=True,
                )
            )

    @staticmethod
    def _is_exact_paper_profile(run: PipelineRun) -> bool:
        if not run.profile.is_paper_method:
            return False
        expected = resolve_profile_for_reference(
            get_profile("paper"),
            run.refname,
        )
        return run.profile == expected

    def _source_code_identity(self):
        code_paths = []
        for relative in (
            "AlphaGEM.py",
            "README.md",
            "AlphaGEM.yml",
            "deepecenvironment.yml",
            "environment.yml",
            "requirements.txt",
        ):
            path = self.repo_root / relative
            if path.is_file():
                code_paths.append(path)
        source_root = self.repo_root / "src"
        if source_root.is_dir():
            code_paths.extend(
                path
                for path in source_root.rglob("*.py")
                if "__pycache__" not in path.parts
            )
        file_records = {
            str(path.relative_to(self.repo_root)): self._file_record(path)
            for path in sorted(set(code_paths))
        }
        aggregate = hashlib.sha256(
            json.dumps(
                file_records,
                sort_keys=True,
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest()
        record = {
            "git_repository": False,
            "commit": None,
            "dirty": None,
            "status_porcelain": [],
            "tracked_diff_sha256": None,
            "code_files_sha256": aggregate,
            "files": file_records,
        }
        if not (self.repo_root / ".git").exists():
            return record
        try:
            commit = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=self.repo_root,
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
            status_lines = subprocess.run(
                ["git", "status", "--porcelain=v1", "--untracked-files=all"],
                cwd=self.repo_root,
                check=True,
                capture_output=True,
                text=True,
            ).stdout.splitlines()
            tracked_diff = subprocess.run(
                [
                    "git",
                    "diff",
                    "--binary",
                    "HEAD",
                    "--",
                    "AlphaGEM.py",
                    "README.md",
                    "AlphaGEM.yml",
                    "deepecenvironment.yml",
                    "environment.yml",
                    "requirements.txt",
                    "src",
                ],
                cwd=self.repo_root,
                check=True,
                capture_output=True,
            ).stdout
        except (OSError, subprocess.CalledProcessError) as error:
            raise RuntimeError(
                f"Cannot establish source-code identity: {error}"
            ) from error
        record.update(
            {
                "git_repository": True,
                "commit": commit,
                "dirty": bool(status_lines),
                "status_porcelain": status_lines,
                "tracked_diff_sha256": hashlib.sha256(
                    tracked_diff
                ).hexdigest(),
            }
        )
        return record

    def _runtime_environment(self):
        packages = {}
        for logical_name, distribution in self.RUNTIME_DISTRIBUTIONS.items():
            try:
                packages[logical_name] = importlib.metadata.version(
                    distribution
                )
            except importlib.metadata.PackageNotFoundError:
                packages[logical_name] = None
        return {
            "python": {
                "implementation": sys.implementation.name,
                "version": sys.version,
                "executable": sys.executable,
            },
            "packages": packages,
        }

    def _gapfill_sequence_record(self, run: PipelineRun):
        records = {}
        for logical_name, supplied, input_kind in (
            ("raw_genome", run.gapfill_genome, "raw_genome"),
            (
                "coordinate_map",
                run.gapfill_coordinate_map,
                "protein_to_genome_coordinate_map",
            ),
            (
                "coordinate_segments",
                run.gapfill_coordinate_segments,
                "protein_to_genome_cds_segments",
            ),
            (
                "coordinate_contract",
                run.gapfill_coordinate_contract,
                "coordinate_recovery_contract",
            ),
            ("nucleotide_cds_proxy", run.gapfill_cds, "nucleotide_cds"),
        ):
            if supplied is None:
                continue
            path = self._resolve_input(supplied)
            records[logical_name] = {
                "observed_path": str(path),
                "input_kind": input_kind,
                **self._file_record(path),
            }
        return records or None

    def _target_input_contract_record(self, run: PipelineRun):
        if run.target_input_contract is None:
            return None
        contract_path = self._resolve_input(run.target_input_contract)
        validated = validate_target_input_contract(
            contract_path,
            self._resolve_input(run.fasta),
        )
        return {
            "observed_path": str(contract_path),
            "bytes": contract_path.stat().st_size,
            **validated,
        }

    def _gapfill_coordinate_contract_record(self, run: PipelineRun):
        if run.gapfill_coordinate_contract is None:
            return None
        if run.target_input_contract is None:
            raise ValueError(
                "A gap-fill coordinate contract requires a target input "
                "contract"
            )
        target_contract = validate_target_input_contract(
            self._resolve_input(run.target_input_contract),
            self._resolve_input(run.fasta),
        )
        contract_path = self._resolve_input(
            run.gapfill_coordinate_contract
        )
        validated = validate_gapfill_coordinate_contract(
            contract_path,
            fasta_path=self._resolve_input(run.fasta),
            genome_path=self._resolve_input(run.gapfill_genome),
            coordinate_map_path=self._resolve_input(
                run.gapfill_coordinate_map
            ),
            coordinate_segments_path=self._resolve_input(
                run.gapfill_coordinate_segments
            ),
            expected_case_id=target_contract["case_id"],
        )
        return {
            "observed_path": str(contract_path),
            "bytes": contract_path.stat().st_size,
            **validated,
        }

    @contextmanager
    def _working_directory(self):
        original = Path.cwd()
        os.chdir(self.repo_root)
        try:
            yield
        finally:
            os.chdir(original)
