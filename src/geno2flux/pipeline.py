"""Single-command genome/proteome-to-ranking and model-screening pipeline.

``geno2flux pipeline --config <json> --output <dir>`` runs Step 0..9b end to end.
Every intermediate artifact can be supplied pre-computed in the config; SHA256
verification then lets the pipeline reuse it instead of rerunning an external
tool.  External tools (gene predictor, DeepLoc, DIAMOND, BLAST+) are invoked
only through configured command interfaces -- nothing external is vendored.

Step map
--------
0   genome -> proteome               (reuse ``proteome_fasta`` or ``predictor_cmd``)
1-3 AlphaGEM homology/pool           (reuse ``pre_gapfill_dir``/``withgaps_model``)
4   localization + compartment fix   (reuse ``localization_tsv`` or ``deeploc_cmd``)
5   HiGHS gapfill on selected medium (AlphaGEM ``gapfilling.gapfill``)
6   mount target reaction / module   (``target_module.mount_target_module``)
7   HiGHS FSEOF 10-step scan         (``published_fseof.run_published_fseof``)
8   paralog family merge + ranking   (DIAMOND self-alignment or cached TSV)
9   knockout screening               (reaction, single-gene, and paralog-family KO)
9b  growth-essentiality filtering    (precomputed TSV or built-in biomass KO scan)

Localization orders:
``pre-gapfill`` (path A, default) refines compartments on the withgaps draft
before gap-filling; ``post-gapfill`` (path B) refines on the already-gapfilled
model; ``both`` runs A then B; ``skip`` runs Steps 6..9b directly on an
already-prepared ``gapfilled_model`` without needing AlphaGEM runtime files.
"""
from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import cobra

from .common import save_json, sha256, table
from .localization import (
    load_localization_table,
    load_models_for_refinement,
    refine_compartments,
    resolve_deeploc_command,
)
from .paralog_grouping import (
    build_families,
    load_self_hits,
    rank_paralog_families,
    run_diamond_self_alignment,
    write_family_tables,
    TIERS,
)
from .knockout import run_ko_suite, run_single_gene_ko, write_ko_tsv

from .target_module import mount_target_module
from .published_fseof import run_published_fseof, FSEOFConfig
from .essential_filter import (
    FSEOF_GENE_FIELDS,
    FAMILY_FIELDS,
    positional_filtered_rank,
    KO_FAMILY_FIELDS,
    KO_GENE_FIELDS,
    apply_essential_filter,
    essential_family_ids,
    load_essential_gene_ids,
    read_tsv_rows,
    validate_threshold,
    write_filter_manifest_json,
    write_filtered_tsv,
)


class PipelineError(RuntimeError):
    """A pipeline stage failed; the manifest records where."""


def _bundle_root() -> Path:
    return Path(__file__).resolve().parents[2] / "third_party" / "alphagem"


def _resolve(config_dir: Path, value: str | None) -> Path | None:
    if value in (None, ""):
        return None
    p = Path(value)
    return p if p.is_absolute() else (config_dir / p).resolve()


def _check_sha(path: Path, expected: str | None) -> None:
    if expected and sha256(path) != expected:
        raise PipelineError(f"SHA256 mismatch for {path}")


def _reject_output_inside_runtime(output: Path, external_runtime: Path) -> None:
    """Output must never be (or sit inside) the read-only runtime tree."""
    out_r = output.resolve()
    rt_r = external_runtime.resolve()
    if out_r == rt_r or rt_r in out_r.parents:
        raise PipelineError(
            f"output directory {out_r} is inside the read-only external "
            f"runtime {rt_r} - point --output outside the runtime tree"
        )


_REPO_ROOT = Path(__file__).resolve().parents[2]

_ALPHAGEM_MODES = {"plmsearch", "structure_alignment"}


def _looks_like_path(value: str) -> bool:
    return "/" in value or "\\" in value or value.startswith((".", "~"))


def resolve_alphagem_python(
    configured: str | None,
    *,
    base_dir: Path,
    env: Mapping[str, str] | None = None,
) -> tuple[str, str]:
    """Interpreter for fresh AlphaGEM steps 1-3: config -> G2F_ALPHAGEM_PYTHON -> current."""
    environ = os.environ if env is None else env

    def _resolve_one(raw: str, label: str) -> str:
        value = raw.strip()
        if _looks_like_path(value):
            candidate = Path(value).expanduser()
            if not candidate.is_absolute():
                candidate = (base_dir / candidate).resolve()
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return str(candidate)
            raise PipelineError(
                f"{label} points to a missing or non-executable Python "
                f"interpreter: {raw!r}"
            )
        found = shutil.which(value, path=environ.get("PATH"))
        if found:
            return found
        raise PipelineError(
            f"{label} interpreter not found on PATH: {raw!r} - pass an "
            "absolute path instead"
        )

    if configured is not None and str(configured).strip():
        return _resolve_one(str(configured), "alphagem_python"), "configured"
    env_val = (environ.get("G2F_ALPHAGEM_PYTHON") or "").strip()
    if env_val:
        return _resolve_one(env_val, "G2F_ALPHAGEM_PYTHON"), "env_G2F_ALPHAGEM_PYTHON"
    return sys.executable, "current_python"


def _tool_bin_dirs(
    *,
    extra_bin_dirs: Sequence[Path] = (),
    env: Mapping[str, str] | None = None,
    repo_root: Path | None = None,
) -> list[Path]:
    """Drop-in binary dirs in lookup order: config, G2F_TOOLS_BIN, repo tools/bin."""
    environ = os.environ if env is None else env
    root = (repo_root or _REPO_ROOT).resolve()
    dirs = [Path(d).expanduser().resolve() for d in extra_bin_dirs]
    env_dir = (environ.get("G2F_TOOLS_BIN") or "").strip()
    if env_dir:
        dirs.append(Path(env_dir).expanduser().resolve())
    dirs.append(root / "tools" / "bin")
    seen: set[Path] = set()
    return [d for d in dirs if not (d in seen or seen.add(d))]


def resolve_tool_executable(
    configured: str,
    *,
    extra_bin_dirs: Sequence[Path] = (),
    env: Mapping[str, str] | None = None,
    repo_root: Path | None = None,
) -> str:
    """Resolve a bare command name against drop-in bin dirs, else return it unchanged."""
    if _looks_like_path(configured):
        return configured
    for directory in _tool_bin_dirs(
        extra_bin_dirs=extra_bin_dirs, env=env, repo_root=repo_root
    ):
        candidate = directory / configured
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
    return configured


def _augmented_tool_env(
    *,
    extra_bin_dirs: Sequence[Path] = (),
    repo_root: Path | None = None,
) -> dict[str, str] | None:
    """``os.environ`` copy with existing tool bin dirs on PATH, or None when no-op."""
    original = os.environ.get("PATH", "")
    front = original.split(os.pathsep)[0] if original else ""
    dirs = [
        str(d)
        for d in _tool_bin_dirs(extra_bin_dirs=extra_bin_dirs, repo_root=repo_root)
        if d.is_dir() and str(d) != front
    ]
    if not dirs:
        return None
    augmented = dict(os.environ)
    augmented["PATH"] = os.pathsep.join(dirs + ([original] if original else []))
    return augmented


@contextmanager
def _tool_path_scope(
    *, extra_bin_dirs: Sequence[Path] = (), repo_root: Path | None = None
):
    """Temporarily prepend existing tool bin dirs to PATH for in-process tool calls."""
    original = os.environ.get("PATH", "")
    dirs = [
        str(d)
        for d in _tool_bin_dirs(extra_bin_dirs=extra_bin_dirs, repo_root=repo_root)
        if d.is_dir()
    ]
    if dirs:
        os.environ["PATH"] = os.pathsep.join(dirs + ([original] if original else []))
    try:
        yield
    finally:
        if dirs:
            os.environ["PATH"] = original


_HEX64 = __import__("re").compile(r"^[0-9a-f]{64}$")


def _validated_expected(value: object, label: str, target: str) -> str:
    """Expected hashes must be non-empty 64-hex strings."""
    if not isinstance(value, str) or not _HEX64.match(value):
        raise PipelineError(
            f"malformed expected hash for {label} ({target!r}): expected a "
            "64-character lowercase hex SHA256 string"
        )
    return value


def _gate_file(
    path: Path, expected: str | None, label: str, hash_key: str,
    *, allow_unverified: bool
) -> dict[str, Any]:
    """Identity gate for reused single-file inputs."""
    record: dict[str, Any] = {"path": str(path), "expected_sha256": expected}
    if expected is not None:
        expected = _validated_expected(expected, label, repr(expected))
        record["expected_sha256"] = expected
        _check_sha(path, expected)
        record["identity"] = "verified"
        return record
    if not allow_unverified:
        raise PipelineError(
            f"{label} is reused without an expected hash: set {hash_key} "
            "in the config, or explicitly set allow_unverified_reuse=true "
            "(the manifest will then record identity: unverified_opted_in)"
        )
    print(
        f"[pipeline] WARNING: {label} reused WITHOUT expected hash "
        f"(allow_unverified_reuse=true): identity UNVERIFIED",
        file=sys.stderr,
    )
    record["identity"] = "unverified_opted_in"
    return record


def _gate_dir(
    directory: Path, expected_map: Mapping[str, Any] | None, label: str,
    *, allow_unverified: bool
) -> dict[str, Any]:
    """Identity gate for reused directory inputs, per actually-consumed file."""
    expected_map = expected_map or {}
    inventory: dict[str, Any] = {}
    for f in sorted(directory.iterdir()):
        if not f.is_file():
            continue
        actual = sha256(f)
        exp = expected_map.get(f.name)
        if exp is not None:
            exp = _validated_expected(exp, f"{label}/{f.name}", repr(exp))
            _check_sha(f, exp)
            status = "verified"
        elif allow_unverified:
            print(
                f"[pipeline] WARNING: {label}/{f.name} reused WITHOUT "
                "expected hash (allow_unverified_reuse=true): UNVERIFIED",
                file=sys.stderr,
            )
            status = "unverified_opted_in"
        else:
            raise PipelineError(
                f"{label} file {f.name!r} has no expected hash: add it to "
                f"{label}_sha256_map, or explicitly set "
                "allow_unverified_reuse=true"
            )
        inventory[f.name] = {"sha256": actual, "identity": status}
    declared_missing = sorted(set(expected_map) - set(inventory))
    if declared_missing:
        raise PipelineError(
            f"{label}_sha256_map declares files absent from the directory: "
            f"{declared_missing}"
        )
    return {
        "path": str(directory),
        "identity": "per-file (see files)",
        "files": inventory,
    }


def _record_file(path: Path) -> dict[str, Any]:
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def _read_any_model(path: Path) -> cobra.Model:
    suffix = path.suffix.lower()
    if suffix in {".yml", ".yaml"}:
        return cobra.io.load_yaml_model(str(path))
    if suffix == ".json":
        return cobra.io.load_json_model(str(path))
    return cobra.io.read_sbml_model(str(path))


class Manifest:
    """Append-only per-step provenance ledger."""

    def __init__(self, output: Path) -> None:
        self.path = output / "pipeline_manifest.json"
        self.steps: list[dict[str, Any]] = []
        self.started = time.time()

    def add(
        self,
        step: str,
        *,
        inputs: dict[str, Any],
        outputs: dict[str, Any],
        reused: bool,
        elapsed: float,
        detail: dict[str, Any] | None = None,
    ) -> None:
        self.steps.append(
            {
                "step": step,
                "inputs": inputs,
                "outputs": outputs,
                "reused": reused,
                "elapsed_seconds": round(elapsed, 3),
                **({"detail": detail} if detail else {}),
            }
        )
        self.flush()

    def flush(self) -> None:
        save_json(
            self.path,
            {
                "steps": self.steps,
                "total_elapsed_seconds": round(time.time() - self.started, 3),
            },
        )


def _register_highs() -> None:
    src = _bundle_root() / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
    from optlang_highs_interface import register_highs_interface

    register_highs_interface()


def _use_default_solver(name: str) -> None:
    """Register ``name`` and make it the default for every new cobra.Model()."""
    if name == "highs":
        _register_highs()
    cobra.Configuration().solver = name


def _prepare_alphagem_runtime(runtime_root: Path, manifest: Manifest) -> Path:
    """Install the vendored AlphaGEM snapshot into ``runtime_root``."""
    script = (
        Path(__file__).resolve().parents[2]
        / "scripts"
        / "prepare_alphagem_runtime.py"
    )
    result = subprocess.run(
        [sys.executable, str(script), "--prepare", "--runtime-root", str(runtime_root)],
        capture_output=True,
        text=True,
    )
    if result.returncode not in (0,):
        raise PipelineError(
            f"prepare_alphagem_runtime failed:\n{result.stdout}\n{result.stderr}"
        )
    return runtime_root


def _ensure_data_available(
    runtime_root: Path, refname: str, reference_fasta: Path, reference_xlsx: Path
) -> None:
    data_dir = runtime_root / "data_available"
    data_dir.mkdir(parents=True, exist_ok=True)
    fasta = data_dir / f"{refname}.fasta"
    xlsx = data_dir / f"{refname}.xlsx"
    if not fasta.exists():
        shutil.copy2(reference_fasta, fasta)
    if not xlsx.exists():
        shutil.copy2(reference_xlsx, xlsx)


@contextmanager
def _chdir(path: Path):
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def _src_module_origin(module: object) -> str | None:
    file = getattr(module, "__file__", None)
    if file:
        return str(Path(file).resolve())
    paths = list(getattr(module, "__path__", []) or [])
    return str(Path(paths[0]).resolve()) if paths else None


def _purge_foreign_src_modules(roots: tuple[Path, ...]) -> None:
    """Drop cached src/src.* modules that resolve OUTSIDE the given roots."""
    resolved = [r.resolve() for r in roots]
    for name in [m for m in list(sys.modules) if m == "src" or m.startswith("src.")]:
        origin = _src_module_origin(sys.modules[name])
        if origin is None:
            del sys.modules[name]
            continue
        op = Path(origin)
        if not any(op == r or r in op.parents for r in resolved):
            del sys.modules[name]


def _repo_imports(runtime_root: Path, *, prefer_runtime: bool = False) -> None:
    """Make src.* and the bundled top-level modules importable."""
    if prefer_runtime:
        bundle = str(_bundle_root())
        sys.path[:] = [e for e in sys.path if e != bundle]
        _purge_foreign_src_modules((runtime_root,))
        s0 = str(runtime_root)
        if s0 not in sys.path:
            sys.path.insert(0, s0)
        return
    _purge_foreign_src_modules((runtime_root, _bundle_root()))
    for entry in (runtime_root, _bundle_root()):
        s = str(entry)
        if s not in sys.path:
            sys.path.insert(0, s)


def _write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {k: (";".join(v) if isinstance(v, list) else v) for k, v in row.items()}
            )


def _gene_report(
    proteins: list[dict[str, Any]],
    reaction_rows: list[dict[str, Any]],
    model_genes: set[str],
) -> list[dict[str, Any]]:
    """One row per input protein; only candidate reactions contribute slopes."""
    indexed: dict[str, list[dict[str, Any]]] = {}
    for reaction in reaction_rows:
        for gid in reaction["gene_ids"]:
            indexed.setdefault(gid, []).append(reaction)
    rows: list[dict[str, Any]] = []
    for protein in proteins:
        gid = protein["target_protein_id"]
        linked = indexed.get(gid, [])
        raw = [r for r in linked if r["raw_candidate"]]
        robust = [r for r in linked if r["robust_candidate"]]
        rows.append(
            {
                "gene_id": gid,
                "in_model": gid in model_genes,
                "raw_score": (sum(r["slope"] for r in raw) / len(raw)) if raw else None,
                "robust_score": (
                    sum(r["slope"] for r in robust) / len(robust)
                )
                if robust
                else None,
                "raw_candidate_reactions": ";".join(r["reaction_id"] for r in raw),
                "robust_candidate_reactions": ";".join(
                    r["reaction_id"] for r in robust
                ),
                "all_associated_reactions": ";".join(
                    r["reaction_id"] for r in linked
                ),
            }
        )
    for name in ("raw", "robust"):
        scores = sorted(
            {r[f"{name}_score"] for r in rows if r[f"{name}_score"] is not None},
            reverse=True,
        )
        ranks = {
            score: 1
            + sum(
                r[f"{name}_score"] > score
                for r in rows
                if r[f"{name}_score"] is not None
            )
            for score in scores
        }
        for row in rows:
            row[f"{name}_rank"] = ranks.get(row[f"{name}_score"])
    rows.sort(key=lambda r: (r["raw_rank"] is None, r["raw_rank"] or 0, r["gene_id"]))
    return rows


def _load_protein_map(id_map: Path | None, proteome_fasta: Path) -> list[dict[str, Any]]:
    """Return the protein roster used for the all-gene FSEOF report."""
    if id_map is not None and id_map.is_file():
        with id_map.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    from Bio import SeqIO

    from .paralog_grouping import protein_key

    return [
        {"target_protein_id": protein_key(record.id)}
        for record in SeqIO.parse(str(proteome_fasta), "fasta")
    ]


def _resolve_target_for_model(
    model: cobra.Model, cfg: Mapping[str, Any], cdir: Path
) -> tuple[str, dict[str, Any]]:
    """Mount or verify the configured target reaction/module on ``model``."""
    target_module_spec = cfg.get("target_module")
    target_rxn_cfg = cfg.get("target_reaction")
    receipt = mount_target_module(
        model,
        target_module_spec,
        target_reaction=target_rxn_cfg,
        base_dir=cdir,
    )
    return str(receipt["product_reaction"]), receipt


def run_pipeline(
    config_path: Path,
    output: Path,
    *,
    solver_override: str | None = None,
    localization_order_override: str | None = None,
) -> dict[str, Any]:
    """Run with the configured solver layer and restore the caller's state."""
    from . import highspy_backend

    previous_solver = cobra.Configuration().solver
    with highspy_backend.scoped_highspy_backend():
        try:
            return _run_pipeline_inner(
                config_path,
                output,
                solver_override=solver_override,
                localization_order_override=localization_order_override,
            )
        finally:
            cobra.Configuration().solver = previous_solver


def _run_pipeline_inner(
    config_path: Path,
    output: Path,
    *,
    solver_override: str | None = None,
    localization_order_override: str | None = None,
) -> dict[str, Any]:
    config_path = config_path.resolve()
    cfg = json.loads(config_path.read_text())
    cdir = config_path.parent
    output = output.resolve()
    ext_early = cfg.get("alphagem_runtime_root")
    if ext_early is not None:
        _reject_output_inside_runtime(output, cdir / ext_early)
    output.mkdir(parents=True, exist_ok=True)
    manifest = Manifest(output)

    solver = solver_override or cfg.get("solver", "highs")
    order = localization_order_override or cfg.get(
        "localization_order", "pre-gapfill"
    )
    if order not in {"pre-gapfill", "post-gapfill", "both", "skip"}:
        raise PipelineError(
            f"Invalid localization_order: {order!r} "
            "(expected 'pre-gapfill', 'post-gapfill', 'both', or 'skip')"
        )
    prune = bool(cfg.get("prune_relocated_genes", False))
    name = cfg.get("name", "sample")
    refname = cfg.get("refname", "yeast")
    biomass = cfg.get("biomass", "r_4041")
    threads = int(cfg.get("threads", 8))
    scan_cfg = cfg.get("scan", {})
    _use_default_solver(solver)

    audit_backend = cfg.get("audit_solver_backend", "highspy")
    if audit_backend == "highspy":
        from .highspy_backend import install_highspy_backend

        install_highspy_backend(
            threads=int(cfg.get("audit_solver_threads", 1)),
            mip_feasibility_tolerance=float(
                cfg.get("audit_solver_mip_feasibility_tolerance", 1e-9)
            ),
        )
    elif audit_backend in ("", "scipy"):
        pass
    else:
        raise PipelineError(
            f"Invalid audit_solver_backend: {audit_backend!r} "
            "(expected 'highspy' [default] or 'scipy' to opt out)"
        )

    def res(key: str) -> Path | None:
        return _resolve(cdir, cfg.get(key))

    proteome_fasta = res("proteome_fasta")
    genome_fasta = res("genome_fasta")
    localization_tsv = res("localization_tsv") or res("deeploc_csv")
    withgaps_model = res("withgaps_model")
    pre_gapfill_dir = res("pre_gapfill_dir")
    gapfilled_model_input = res("gapfilled_model")
    template_override = res("template_model")
    diamond_hits = res("diamond_self_hits")
    id_map = res("id_map")
    essential_tsv = res("essential_genes_tsv")
    structure_maplist = res("structure_maplist")
    structure_dir = res("structure_dir")
    prottrans_model_path = res("prottrans_model_path")
    compute_growth_essentiality = bool(cfg.get("compute_growth_essentiality", False))

    for key, path in (
        ("proteome_fasta", proteome_fasta),
        ("localization_tsv", localization_tsv),
        ("gapfilled_model", gapfilled_model_input),
        ("template_model", template_override),
        ("diamond_self_hits", diamond_hits),
        ("id_map", id_map),
        ("essential_genes_tsv", essential_tsv),
    ):
        if path is not None and not path.is_file():
            raise PipelineError(
                f"{key} is configured but not a readable file: {path}"
            )
    if essential_tsv is not None or compute_growth_essentiality:
        try:
            threshold = validate_threshold(cfg.get("essential_ko_loss_threshold", 0.5))
            if essential_tsv is not None:
                load_essential_gene_ids(essential_tsv, threshold)
        except ValueError as exc:
            raise PipelineError(str(exc)) from exc

    if withgaps_model is not None and not withgaps_model.is_file():
        raise PipelineError(
            f"withgaps_model is configured but not a readable file: "
            f"{withgaps_model} - fix the path or remove the key (fresh "
            "steps 1-3); refusing to guess"
        )
    if pre_gapfill_dir is not None and not pre_gapfill_dir.is_dir():
        raise PipelineError(
            f"pre_gapfill_dir is configured but not a directory: "
            f"{pre_gapfill_dir} - fix the path or remove the key (fresh "
            "steps 1-3); refusing to guess"
        )
    tools_bin_dir = res("tools_bin_dir")
    if tools_bin_dir is not None and not tools_bin_dir.is_dir():
        raise PipelineError(
            f"tools_bin_dir is configured but not a directory: "
            f"{tools_bin_dir} - fix the path or remove the key "
            "(<repo>/tools/bin and G2F_TOOLS_BIN are always searched); "
            "refusing to guess"
        )
    extra_bins: tuple[Path, ...] = (
        (tools_bin_dir,) if tools_bin_dir is not None else ()
    )
    alphagem_mode = cfg.get("alphagem_mode", "plmsearch")
    if alphagem_mode not in _ALPHAGEM_MODES:
        raise PipelineError(
            f"Invalid alphagem_mode: {alphagem_mode!r} "
            "(expected 'plmsearch' or 'structure_alignment')"
        )

    allow_unverified = cfg.get("allow_unverified_reuse", False)
    if not isinstance(allow_unverified, bool):
        raise PipelineError("allow_unverified_reuse must be a JSON boolean")
    reuse_identity: dict[str, Any] = {}
    if proteome_fasta is not None and proteome_fasta.is_file():
        reuse_identity["proteome_fasta"] = _gate_file(
            proteome_fasta, cfg.get("proteome_sha256"), "proteome_fasta",
            "proteome_sha256", allow_unverified=allow_unverified,
        )
    if withgaps_model is not None and withgaps_model.is_file():
        reuse_identity["withgaps_model"] = _gate_file(
            withgaps_model, cfg.get("withgaps_sha256"), "withgaps_model",
            "withgaps_sha256", allow_unverified=allow_unverified,
        )
    if pre_gapfill_dir is not None and pre_gapfill_dir.is_dir():
        reuse_identity["pre_gapfill_dir"] = _gate_dir(
            pre_gapfill_dir,
            cfg.get("pre_gapfill_dir_sha256_map"),
            "pre_gapfill_dir",
            allow_unverified=allow_unverified,
        )
    if localization_tsv is not None and localization_tsv.is_file():
        reuse_identity["localization_tsv"] = _gate_file(
            localization_tsv, cfg.get("localization_sha256"),
            "localization_tsv", "localization_sha256",
            allow_unverified=allow_unverified,
        )
    if gapfilled_model_input is not None and gapfilled_model_input.is_file():
        reuse_identity["gapfilled_model"] = _gate_file(
            gapfilled_model_input, cfg.get("gapfilled_model_sha256"),
            "gapfilled_model", "gapfilled_model_sha256",
            allow_unverified=allow_unverified,
        )
    if template_override is not None and template_override.is_file():
        reuse_identity["template_model"] = _gate_file(
            template_override, cfg.get("template_model_sha256"),
            "template_model", "template_model_sha256",
            allow_unverified=allow_unverified,
        )
    if diamond_hits is not None and diamond_hits.is_file():
        reuse_identity["diamond_self_hits"] = _gate_file(
            diamond_hits, cfg.get("diamond_self_hits_sha256"),
            "diamond_self_hits", "diamond_self_hits_sha256",
            allow_unverified=allow_unverified,
        )

    # Determine whether AlphaGEM runtime staging is needed:
    # - Screening-only mode (`order == "skip"`) with `gapfilled_model` needs no AlphaGEM runtime.
    # - `order == "post-gapfill"` with `gapfilled_model` AND explicit `template_model` also needs no AlphaGEM runtime.
    needs_alphagem_runtime = not (
        (order == "skip" and gapfilled_model_input is not None)
        or (order == "post-gapfill" and gapfilled_model_input is not None and template_override is not None)
    )
    fresh_steps_1_3 = (
        needs_alphagem_runtime
        and withgaps_model is None
        and pre_gapfill_dir is None
    )
    if fresh_steps_1_3 and alphagem_mode == "structure_alignment":
        missing_keys: list[str] = []
        if structure_maplist is None:
            missing_keys.append("structure_maplist")
        elif not structure_maplist.is_file():
            missing_keys.append(
                "structure_maplist (not a readable file: "
                f"{structure_maplist})"
            )
        if structure_dir is None:
            missing_keys.append("structure_dir")
        elif not structure_dir.is_dir():
            missing_keys.append(
                f"structure_dir (not a directory: {structure_dir})"
            )
        elif not any(structure_dir.glob("*.pdb")):
            missing_keys.append("structure_dir (contains no *.pdb files)")
        if prottrans_model_path is None:
            missing_keys.append("prottrans_model_path")
        elif not prottrans_model_path.is_dir():
            missing_keys.append(
                "prottrans_model_path (not a directory: "
                f"{prottrans_model_path})"
            )
        if missing_keys:
            raise PipelineError(
                "alphagem_mode='structure_alignment' with fresh steps 1-3 "
                "(AlphaGEM.py --profile paper) requires config keys: "
                + ", ".join(missing_keys)
            )

    if order == "skip" and gapfilled_model_input is None:
        raise PipelineError(
            "localization_order='skip' requires 'gapfilled_model' to be configured"
        )

    external_runtime = res("alphagem_runtime_root")
    runtime_root: Path | None = None
    working: Path | None = None
    if needs_alphagem_runtime:
        runtime_detail: dict[str, Any] = {}
        if external_runtime is not None:
            _reject_output_inside_runtime(output, external_runtime)
            required = [
                external_runtime / "AlphaGEM.py",
                external_runtime / "src",
                external_runtime / "models" / "yeast-GEM.xml",
                external_runtime / "models" / "reaction_pool_final.json",
            ]
            missing = [str(r) for r in required if not r.exists()]
            if missing:
                raise PipelineError(
                    "alphagem_runtime_root lacks required entries: "
                    f"{missing} — prepare it with scripts/"
                    "prepare_alphagem_runtime.py"
                )
            data_dir = external_runtime / "data_available"
            missing_refs = [
                str(p)
                for p in (
                    data_dir / f"{refname}.fasta",
                    data_dir / f"{refname}.xlsx",
                )
                if not p.is_file()
            ]
            if missing_refs:
                raise PipelineError(
                    "alphagem_runtime_root/data_available lacks prepared "
                    f"references (read-only, nothing is written): {missing_refs}"
                )
            view = output / "alphagem_runtime_view"
            if view.exists() or (output / "alphagem_runtime").exists():
                raise PipelineError(
                    "output already contains a previous run's runtime tree "
                    "(alphagem_runtime_view or alphagem_runtime) — one attempt "
                    "per output directory"
                )
            view.mkdir(parents=True, exist_ok=False)
            if fresh_steps_1_3:
                shutil.copytree(
                    external_runtime,
                    view,
                    dirs_exist_ok=True,
                    symlinks=False,
                    ignore=shutil.ignore_patterns("__pycache__", "working", "temp", "var"),
                )
            else:
                shutil.copytree(external_runtime / "src", view / "src")
                for entry in ("AlphaGEM.py", "models"):
                    (view / entry).symlink_to(
                        external_runtime / entry,
                        target_is_directory=entry != "AlphaGEM.py",
                    )
                shutil.copytree(
                    external_runtime / "data_available", view / "data_available"
                )
            runtime_root = view
            working = view / "working" / name
            working.mkdir(parents=True, exist_ok=True)
            (view / "temp").mkdir(exist_ok=True)
            runtime_detail.update(
                runtime_source="external_configured_runtime_read_only",
                fresh_steps_1_3=(
                    "view is a COMPLETE per-run writable copy of the prepared "
                    "runtime (entry resolves REPO_ROOT inside the view; tools/ "
                    "rhea/ plmsearchtools/ and source identity included); "
                    "upstream tool EXECUTION not validated in this environment"
                    if fresh_steps_1_3
                    else False
                ),
                external_runtime_root=str(external_runtime),
                per_run_view=str(view),
                view_contract=(
                    "fresh: full per-run copy (caches/scratch excluded); "
                    "reuse: src+data_available copied, AlphaGEM.py/models "
                    "symlinked; all writes land in the per-run view"
                ),
            )
            si = external_runtime / "SOURCE_IDENTITY.json"
            if si.is_file():
                runtime_detail["SOURCE_IDENTITY_sha256"] = sha256(si)
        else:
            runtime_root = output / "alphagem_runtime"
            if runtime_root.exists() or (output / "alphagem_runtime_view").exists():
                raise PipelineError(
                    "output already contains a previous run's runtime tree — "
                    "one attempt per output directory; use a fresh output"
                )
            _prepare_alphagem_runtime(runtime_root, manifest)
            working = runtime_root / "working" / name
            working.mkdir(parents=True, exist_ok=True)
            (runtime_root / "temp").mkdir(exist_ok=True)
            runtime_detail["runtime_source"] = "vendored_bundle_copy"
        runtime_detail["resource_verification"] = (
            "runtime resource lock check not wired into the unified pipeline; "
            "see requirements.txt 'AlphaGEM runtime resources' (lock-resources)"
        )
        t0 = time.time()
        reference_fasta = res("reference_fasta")
        reference_xlsx = res("reference_xlsx")
        if reference_fasta is None or reference_xlsx is None:
            raise PipelineError(
                "reference_fasta and reference_xlsx (data_available/<refname>.fasta/.xlsx) are required"
            )
        if external_runtime is None:
            _ensure_data_available(
                runtime_root, refname, reference_fasta, reference_xlsx
            )
        manifest.add(
            "runtime_prepare",
            inputs={},
            outputs={"runtime_root": str(runtime_root)},
            reused=False,
            elapsed=time.time() - t0,
            detail=runtime_detail,
        )
    else:
        if (output / "pipeline_result.json").exists():
            raise PipelineError(
                "output already contains a previous run's pipeline_result.json; use a fresh output"
            )

    # ----- Step 0: proteome -----
    t0 = time.time()
    if proteome_fasta is not None and proteome_fasta.is_file():
        manifest.add(
            "step0_proteome",
            inputs={
                "proteome_fasta": _record_file(proteome_fasta),
                "identity": reuse_identity["proteome_fasta"],
            },
            outputs={"proteome_fasta": str(proteome_fasta)},
            reused=True,
            elapsed=time.time() - t0,
        )
    elif genome_fasta is not None and cfg.get("predictor_cmd"):
        proteome_fasta = output / "proteome.fasta"
        cmd = list(cfg["predictor_cmd"]) + [str(genome_fasta), str(proteome_fasta)]
        subprocess.run(cmd, check=True)
        manifest.add(
            "step0_proteome",
            inputs={"genome_fasta": _record_file(genome_fasta)},
            outputs={"proteome_fasta": _record_file(proteome_fasta)},
            reused=False,
            elapsed=time.time() - t0,
        )
    else:
        raise PipelineError(
            "No proteome_fasta supplied and no predictor_cmd configured for Step 0"
        )
    assert proteome_fasta is not None

    # ----- Steps 1-3: AlphaGEM through reaction_pool -> withgaps -----
    withgaps_in_working: Path | None = None
    pristine_withgaps: Path | None = None
    if needs_alphagem_runtime:
        assert runtime_root is not None and working is not None
        t0 = time.time()
        withgaps_in_working = working / f"{name}-GEM_withgaps.yml"
        if withgaps_model is not None and withgaps_model.is_file():
            shutil.copy2(withgaps_model, withgaps_in_working)
            shutil.copy2(proteome_fasta, working / f"{name}.fasta")
            manifest.add(
                "step1_3_withgaps",
                inputs={
                    "withgaps_model": _record_file(withgaps_model),
                    "identity": reuse_identity["withgaps_model"],
                },
                outputs={"withgaps": _record_file(withgaps_in_working)},
                reused=True,
                elapsed=time.time() - t0,
            )
        elif pre_gapfill_dir is not None and pre_gapfill_dir.is_dir():
            src = pre_gapfill_dir / f"{name}-GEM_withgaps.yml"
            if not src.is_file():
                raise PipelineError(f"pre_gapfill_dir lacks {name}-GEM_withgaps.yml")
            for f in pre_gapfill_dir.iterdir():
                if f.is_file():
                    shutil.copy2(f, working / f.name)
            shutil.copy2(proteome_fasta, working / f"{name}.fasta")
            manifest.add(
                "step1_3_withgaps",
                inputs={"identity": reuse_identity["pre_gapfill_dir"]},
                outputs={"withgaps": _record_file(withgaps_in_working)},
                reused=True,
                elapsed=time.time() - t0,
            )
        else:
            alphagem_python, alphagem_python_source = resolve_alphagem_python(
                cfg.get("alphagem_python"), base_dir=cdir
            )
            cmd = [
                alphagem_python,
                str(runtime_root / "AlphaGEM.py"),
                "--profile", "paper",
                "--mode", alphagem_mode,
                "--refname", refname,
                "--name", name,
                "--fasta", str(proteome_fasta),
                "--solver", solver,
                "--stop-after", "reaction_pool",
            ]
            for flag, key in (
                ("--target-input-contract", "target_input_contract"),
                ("--esm-model-path", "esm_model_path"),
                ("--deepec-python", "deepec_python"),
            ):
                value = res(key)
                if value is not None:
                    cmd += [flag, str(value)]
            if cfg.get("esm_token_budget") is not None:
                raw_budget = cfg["esm_token_budget"]
                if isinstance(raw_budget, bool) or not isinstance(raw_budget, int) or raw_budget <= 0:
                    raise PipelineError(
                        f"esm_token_budget must be a positive integer, got {raw_budget!r}"
                    )
                cmd += ["--esm-token-budget", str(raw_budget)]
            blast_threads = cfg.get("blast_threads")
            if blast_threads is not None:
                if isinstance(blast_threads, bool) or not isinstance(blast_threads, int) or blast_threads <= 0:
                    raise PipelineError(
                        f"blast_threads must be a positive integer, got {blast_threads!r}"
                    )
                cmd += ["--blast-threads", str(blast_threads)]
            if alphagem_mode == "structure_alignment":
                for flag, key in (
                    ("--prottrans-model-path", "prottrans_model_path"),
                    ("--maplist", "structure_maplist"),
                    ("--structure", "structure_dir"),
                ):
                    value = res(key)
                    if value is not None:
                        cmd += [flag, str(value)]
            tool_env = _augmented_tool_env(extra_bin_dirs=extra_bins)
            run_kwargs: dict[str, Any] = {"cwd": runtime_root, "check": True}
            if tool_env is not None:
                run_kwargs["env"] = tool_env
            try:
                subprocess.run(cmd, **run_kwargs)
            except subprocess.CalledProcessError as exc:
                message = (
                    f"AlphaGEM steps 1-3 subprocess exited with code "
                    f"{exc.returncode}: {cmd!r}"
                )
                if alphagem_python_source == "current_python":
                    message += (
                        " - fresh AlphaGEM reconstruction ran under the "
                        "current Geno2Flux Python interpreter, which lacks "
                        "AlphaGEM's deep-learning stack (PyTorch, fair-esm, "
                        "CLEAN, ...). Set \"alphagem_python\" in the config or "
                        "export G2F_ALPHAGEM_PYTHON to a Python that has it "
                        "installed, then run 'bash run.sh check-env "
                        "--runtime-root <dir>' to verify."
                    )
                raise PipelineError(message) from exc
            manifest.add(
                "step1_3_withgaps",
                inputs={"proteome_fasta": _record_file(proteome_fasta)},
                outputs={"withgaps": _record_file(withgaps_in_working)},
                reused=False,
                elapsed=time.time() - t0,
                detail={
                    "alphagem_python": alphagem_python,
                    "alphagem_python_source": alphagem_python_source,
                    "alphagem_mode": alphagem_mode,
                },
            )
        pristine_withgaps = output / "_pristine" / withgaps_in_working.name
        pristine_withgaps.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(withgaps_in_working, pristine_withgaps)

    # ----- Step 4 preparation: localization table -----
    localization: dict[str, set[str]] = {}
    loc_reused = False
    loc_status = "skipped"
    loc_source = "none"
    require_deeploc = bool(cfg.get("require_deeploc", False))
    if order != "skip":
        if localization_tsv is not None and localization_tsv.is_file():
            localization = load_localization_table(localization_tsv)
            loc_reused = True
            loc_status = "reused_table"
            loc_source = "localization_tsv"
        else:
            try:
                resolved_deeploc, loc_source = resolve_deeploc_command(
                    cfg.get("deeploc_cmd")
                )
            except FileNotFoundError as exc:
                raise PipelineError(str(exc)) from exc
            if resolved_deeploc is not None:
                deeploc_out = output / "deeploc_out"
                deeploc_out.mkdir(exist_ok=True)
                cmd = list(resolved_deeploc) + [
                    "-f",
                    str(proteome_fasta),
                    "-o",
                    str(deeploc_out),
                ]
                subprocess.run(cmd, check=True)
                csvs = sorted(deeploc_out.glob("results_*.csv"))
                if not csvs:
                    raise PipelineError("deeploc2 produced no results_*.csv")
                localization = load_localization_table(csvs[-1])
                localization_tsv = csvs[-1]
                loc_reused = False
                loc_status = "predicted_with_deeploc"
            elif require_deeploc:
                raise PipelineError(
                    "DeepLoc was required (require_deeploc=true), but no "
                    "localization_tsv was supplied and no DeepLoc binary was found "
                    "(checked deeploc_cmd, G2F_DEEPLOC_BIN, tools/deeploc/, and PATH)"
                )
            else:
                loc_status = "fallback_default_compartments"
                loc_source = "not_found_fallback_default"
                print(
                    "[geno2flux] DeepLoc not detected and no localization_tsv "
                    "supplied; keeping default model compartments (install DeepLoc "
                    "into tools/deeploc/ or set G2F_DEEPLOC_BIN to enable refinement).",
                    file=sys.stderr,
                )

    template_path = (
        template_override
        if template_override is not None
        else (runtime_root / "models" / "yeast-GEM.xml" if runtime_root is not None else None)
    )

    results: dict[str, Any] = {"paths": {}, "manifest": str(manifest.path)}
    orders = (
        ["pre-gapfill", "post-gapfill"]
        if order == "both"
        else (["screening"] if order == "skip" else [order])
    )
    for od in orders:
        path_dir = output / od
        path_dir.mkdir(exist_ok=True)

        if od == "screening":
            assert gapfilled_model_input is not None
            t0 = time.time()
            final_xml = path_dir / f"{name}-GEM_input.xml"
            if gapfilled_model_input.suffix.lower() == ".xml":
                shutil.copy2(gapfilled_model_input, final_xml)
            else:
                loaded_in = _read_any_model(gapfilled_model_input)
                cobra.io.write_sbml_model(loaded_in, str(final_xml))
            refined_path = final_xml
            audit = {"pulled_reaction_count": 0, "skipped_localization": True}
            manifest.add(
                "step5_model_input[screening]",
                inputs={
                    "gapfilled_model": _record_file(gapfilled_model_input),
                    "identity": reuse_identity["gapfilled_model"],
                },
                outputs={"final_model": _record_file(final_xml)},
                reused=True,
                elapsed=time.time() - t0,
                detail={"mode": "screening_only"},
            )
        else:
            assert template_path is not None
            t0 = time.time()
            if od == "pre-gapfill":
                assert withgaps_in_working is not None
                base_model_path = withgaps_in_working
            else:
                if gapfilled_model_input is not None:
                    base_model_path = gapfilled_model_input
                else:
                    assert pristine_withgaps is not None and withgaps_in_working is not None
                    assert runtime_root is not None and working is not None
                    shutil.copy2(pristine_withgaps, withgaps_in_working)
                    gapfill_base_dir = output / "_post_gapfill_base"
                    gapfill_base_dir.mkdir(exist_ok=True)
                    gapfilled_base = _run_gapfill(
                        name=name,
                        refname=refname,
                        runtime_root=runtime_root,
                        prefer_runtime=external_runtime is not None,
                        working=working,
                        cfg=cfg,
                        cdir=cdir,
                        solver=solver,
                        extra_bin_dirs=extra_bins,
                        out_dir=gapfill_base_dir,
                        tag="post_base",
                    )
                    base_model_path = gapfilled_base
            if loc_status == "fallback_default_compartments":
                target = _read_any_model(base_model_path)
                audit = {
                    "pulled_reaction_count": 0,
                    "pulled_reactions": [],
                    "prune_relocated_genes": prune,
                    "final_reaction_count": len(target.reactions),
                    "deeploc_status": loc_status,
                    "deeploc_source": loc_source,
                    "fallback": "deeploc_unavailable_default_compartments",
                }
            else:
                target, template = load_models_for_refinement(base_model_path, template_path)
                audit = refine_compartments(
                    target, template, localization, prune_relocated_genes=prune
                )
                audit["deeploc_status"] = loc_status
                audit["deeploc_source"] = loc_source
            if od == "pre-gapfill":
                refined_yml = path_dir / f"{name}-GEM_withgaps_loc_refined.yml"
                cobra.io.save_yaml_model(target, str(refined_yml))
                refined_path = refined_yml
            else:
                refined_xml = path_dir / f"{name}-GEM_loc_post_gapfill.xml"
                cobra.io.write_sbml_model(target, str(refined_xml))
                refined_path = refined_xml
            save_json(path_dir / "localization_refinement_audit.json", audit)
            step4_inputs: dict[str, Any] = {
                "base_model": _record_file(base_model_path),
                "template": _record_file(template_path),
            }
            if localization_tsv is not None:
                step4_inputs["localization_tsv"] = _record_file(localization_tsv)
                step4_inputs["identity"] = reuse_identity.get(
                    "localization_tsv", {"identity": "fresh_prediction"}
                )
            else:
                step4_inputs["localization_tsv"] = None
                step4_inputs["identity"] = {"identity": "fallback_default_compartments"}
            manifest.add(
                f"step4_localization[{od}]",
                inputs=step4_inputs,
                outputs={
                    "refined_model": _record_file(refined_path),
                    "audit": _record_file(path_dir / "localization_refinement_audit.json"),
                },
                reused=loc_reused,
                elapsed=time.time() - t0,
                detail={
                    "pulled": audit["pulled_reaction_count"],
                    "order": od,
                    "deeploc_status": loc_status,
                    "deeploc_source": loc_source,
                },
            )

            # ----- Step 5: HiGHS gapfill -----
            t0 = time.time()
            gapfill_input_yml = refined_path
            if withgaps_in_working is not None:
                shutil.copy2(refined_path, withgaps_in_working)
            final_xml = path_dir / f"{name}-GEM_loc_{od.replace('-gapfill','_gapfill')}.xml"
            if od == "post-gapfill" and gapfilled_model_input is not None:
                final_model = target
                cobra.io.write_sbml_model(final_model, str(final_xml))
                cobra.io.save_json_model(
                    final_model, str(path_dir / f"{name}-GEM_loc_post_gapfill.json")
                )
                gapfill_audit = {"reused": "already_gapfilled_input"}
                save_json(path_dir / "gapfill_audit.json", gapfill_audit)
            elif od == "post-gapfill" and gapfilled_model_input is None:
                final_xml = refined_path.with_suffix(".xml")
                cobra.io.write_sbml_model(target, str(final_xml))
                cobra.io.save_json_model(
                    target, str(path_dir / f"{name}-GEM_loc_post_gapfill.json")
                )
                base_audit = output / "_post_gapfill_base" / "gapfill_audit.json"
                if base_audit.is_file():
                    shutil.copy2(base_audit, path_dir / "gapfill_audit.json")
                else:
                    save_json(path_dir / "gapfill_audit.json", {"reused": "post_gapfill_base"})
            else:
                assert runtime_root is not None and working is not None
                final_xml = _run_gapfill(
                    name=name,
                    refname=refname,
                    runtime_root=runtime_root,
                    prefer_runtime=external_runtime is not None,
                    working=working,
                    cfg=cfg,
                    cdir=cdir,
                    solver=solver,
                    extra_bin_dirs=extra_bins,
                    out_dir=path_dir,
                    tag="pre" if od == "pre-gapfill" else "post",
                )
            gapfill_reused = (
                od == "post-gapfill" and gapfilled_model_input is not None
            )
            manifest.add(
                f"step5_gapfill[{od}]",
                inputs={
                    "model": _record_file(gapfill_input_yml),
                    **(
                        {"identity": reuse_identity["gapfilled_model"]}
                        if gapfill_reused
                        else {}
                    ),
                },
                outputs={"final_model": _record_file(final_xml)},
                reused=gapfill_reused,
                elapsed=time.time() - t0,
                detail=(
                    {"reused_input": "already_gapfilled_model"}
                    if gapfill_reused
                    else {"rerun": "alphagem_gapfilling"}
                ),
            )

        # ----- Step 6: mount target reaction / product module -----
        t0 = time.time()
        model = cobra.io.read_sbml_model(str(final_xml))
        target_rxn_id, target_receipt = _resolve_target_for_model(model, cfg, cdir)
        _attach_solver(model, solver)
        target_xml = path_dir / f"target_model_{od.replace('-', '_')}.xml"
        cobra.io.write_sbml_model(model, str(target_xml))
        save_json(path_dir / "target_module_receipt.json", target_receipt)
        manifest.add(
            f"step6_target[{od}]",
            inputs={"model": _record_file(final_xml)},
            outputs={"target_model": _record_file(target_xml)},
            reused=False,
            elapsed=time.time() - t0,
            detail={
                "target_reaction": target_rxn_id,
                "reaction_count": len(model.reactions),
            },
        )

        # ----- Step 7: HiGHS FSEOF -----
        t0 = time.time()
        fseof_dir = path_dir / "fseof"
        fseof_dir.mkdir(exist_ok=True)
        fseof_result = run_published_fseof(
            model,
            biomass,
            target_rxn_id,
            config=scan_cfg or None,
            source_identity={"pipeline_step": "step7_fseof", "order": od},
        )
        save_json(fseof_dir / "fseof_full_witnesses.json", fseof_result)
        proteins = _load_protein_map(id_map, proteome_fasta)
        gene_rows = _gene_report(
            proteins, fseof_result["reaction_rows"], {g.id for g in model.genes}
        )
        _write_tsv(
            fseof_dir / "all_genes_fseof.tsv",
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
        _write_tsv(
            fseof_dir / "all_reactions_fseof.tsv",
            fseof_result["reaction_rows"],
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
        save_json(
            fseof_dir / "fseof_summary.json",
            {
                "order": od,
                "solver": solver,
                "biomass_reaction": biomass,
                "target_reaction": target_rxn_id,
                "maximum_product": fseof_result["maximum_product"],
                "raw_candidates": fseof_result["raw_candidate_count"],
                "robust_candidates": fseof_result["robust_candidate_count"],
                "scan_growth": [s["point_optimal_growth"] for s in fseof_result["scans"]],
                "cobra": cobra.__version__,
            },
        )
        manifest.add(
            f"step7_fseof[{od}]",
            inputs={"target_model": _record_file(target_xml)},
            outputs={
                "all_genes_fseof": _record_file(fseof_dir / "all_genes_fseof.tsv"),
                "all_reactions_fseof": _record_file(
                    fseof_dir / "all_reactions_fseof.tsv"
                ),
            },
            reused=False,
            elapsed=time.time() - t0,
            detail={"raw_candidates": fseof_result["raw_candidate_count"]},
        )

        # ----- Step 8: paralog families -----
        t0 = time.time()
        para_dir = path_dir / "paralog"
        para_dir.mkdir(exist_ok=True)
        if diamond_hits is not None and diamond_hits.is_file():
            hits_path = diamond_hits
            hits_reused = True
        else:
            hits_path = para_dir / "self_hits.tsv"
            run_diamond_self_alignment(
                proteome_fasta,
                hits_path,
                diamond_bin=resolve_tool_executable(
                    cfg.get("diamond_bin", "diamond"), extra_bin_dirs=extra_bins
                ),
                threads=threads,
                workdir=para_dir,
            )
            hits_reused = False
        hits = load_self_hits(hits_path)
        gene_scores = {
            row["gene_id"]: (
                float(row["raw_score"]) if row["raw_score"] is not None else None
            )
            for row in gene_rows
        }
        all_families = {}
        for tier_name, tier in TIERS.items():
            fams = build_families(hits, tier)
            all_families[tier_name] = rank_paralog_families(
                gene_scores, fams, all_protein_ids=gene_scores.keys()
            )
        write_family_tables(
            all_families["default"],
            para_dir / "paralog_family_rankings.json",
            para_dir / "paralog_family_rankings.tsv",
        )
        save_json(para_dir / "paralog_family_rankings_all_tiers.json", all_families)
        manifest.add(
            f"step8_paralog[{od}]",
            inputs={
                "self_hits": _record_file(hits_path),
                "proteome": _record_file(proteome_fasta),
                **(
                    {"identity": reuse_identity["diamond_self_hits"]}
                    if hits_reused
                    else {}
                ),
            },
            outputs={
                "family_rankings": _record_file(
                    para_dir / "paralog_family_rankings.json"
                )
            },
            reused=hits_reused,
            elapsed=time.time() - t0,
            detail={"families": len(all_families["default"]), "hits": len(hits)},
        )

        # ----- Step 9: knockout prediction -----
        t0 = time.time()
        ko_dir = path_dir / "knockout"
        ko_dir.mkdir(exist_ok=True)
        ko_families = {}
        for fam in all_families["default"]:
            ko_families[fam["family_id"]] = fam.get("members", [])
        for fam in all_families["default"]:
            fid = fam["family_id"]
            if fid not in ko_families or sorted(ko_families[fid]) != sorted(fam["members"]):
                raise PipelineError(
                    f"Family KO input mismatch for {fid}: "
                    f"ko_families={sorted(ko_families.get(fid,[]))} "
                    f"vs default_tier={sorted(fam['members'])}"
                )
        ko_result = run_ko_suite(
            model,
            target_rxn_id,
            families=ko_families,
            out_dir=ko_dir,
        )
        save_json(ko_dir / "ko_summary.json", ko_result)
        manifest.add(
            f"step9_knockout[{od}]",
            inputs={"target_model": _record_file(target_xml)},
            outputs={
                "single_reaction_ko": _record_file(ko_dir / "single_reaction_ko.tsv"),
                "single_gene_ko": _record_file(ko_dir / "single_gene_ko.tsv"),
                "family_ko": _record_file(ko_dir / "family_ko.tsv"),
            },
            reused=False,
            elapsed=time.time() - t0,
            detail={"top_reaction_loss": ko_result.get("top_reaction_loss")},
        )

        # ----- Optional built-in growth single-gene KO for Step 9b -----
        effective_essential_tsv = essential_tsv
        if effective_essential_tsv is None and compute_growth_essentiality:
            t0g = time.time()
            growth_opt, growth_nonopt = run_single_gene_ko(model, biomass, ko_dir)
            effective_essential_tsv = ko_dir / "single_gene_ko_growth.tsv"
            write_ko_tsv(effective_essential_tsv, growth_opt, "gene_id")
            write_ko_tsv(
                ko_dir / "single_gene_ko_growth_nonoptimal.tsv",
                growth_nonopt,
                "gene_id",
                nonoptimal=True,
            )
            manifest.add(
                f"step9a_growth_essentiality[{od}]",
                inputs={"target_model": _record_file(target_xml)},
                outputs={"single_gene_ko_growth": _record_file(effective_essential_tsv)},
                reused=False,
                elapsed=time.time() - t0g,
                detail={
                    "biomass_reaction": biomass,
                    "optimal_rows": len(growth_opt),
                    "nonoptimal_rows": len(growth_nonopt),
                },
            )

        # ----- Step 9b: essentiality-filtered rankings -----
        essential_summary = _write_essential_filtered_outputs(
            fseof_dir=fseof_dir,
            para_dir=para_dir,
            ko_dir=ko_dir,
            gene_rows=gene_rows,
            family_rows=all_families["default"],
            essential_tsv=effective_essential_tsv,
            cfg=cfg,
            manifest=manifest,
            od=od,
        )

        results["paths"][od] = {
            "refined_model": str(refined_path),
            "final_model": str(final_xml),
            "target_model": str(target_xml),
            "target_reaction": target_rxn_id,
            "fseof_dir": str(fseof_dir),
            "paralog_dir": str(para_dir),
            "ko_dir": str(ko_dir),
            "pulled": audit["pulled_reaction_count"],
            "final_reactions": len(model.reactions),
            "essential_filter": essential_summary,
        }
    manifest.flush()
    save_json(output / "pipeline_result.json", results)
    return results


def _write_essential_filtered_outputs(
    *,
    fseof_dir: Path,
    para_dir: Path,
    ko_dir: Path,
    gene_rows: list[dict[str, Any]],
    family_rows: list[dict[str, Any]],
    essential_tsv: Path | None,
    cfg: Mapping[str, Any],
    manifest: "Manifest",
    od: str,
) -> dict[str, Any]:
    """Step 9b: essentiality-filtered rankings (growth-derived rule)."""
    if essential_tsv is None:
        return {"essentiality_source": "not_configured"}
    if not essential_tsv.is_file():
        raise PipelineError(
            f"essential_genes_tsv is configured but not a readable "
            f"file: {essential_tsv}"
        )
    try:
        threshold = validate_threshold(
            cfg.get("essential_ko_loss_threshold", 0.5)
        )
    except ValueError as exc:
        raise PipelineError(str(exc)) from exc
    record = load_essential_gene_ids(essential_tsv, threshold)
    essential_genes = record["ids"]
    essential_fams = essential_family_ids(family_rows, essential_genes)
    t0e = time.time()
    gres = apply_essential_filter(
        gene_rows, essential_genes, "gene_id", "raw_rank",
        order="gene_v2",
    )
    write_filtered_tsv(
        fseof_dir / "all_genes_fseof_essential_filtered.tsv",
        gres["kept"],
        FSEOF_GENE_FIELDS,
    )
    fres = apply_essential_filter(
        family_rows, essential_fams, "family_id", "max_score_rank"
    )
    positional_filtered_rank(fres["kept"], "max_score_rank", "family_id",
                             "filtered_rank_max")
    positional_filtered_rank(fres["kept"], "mean_score_rank", "family_id",
                             "filtered_rank_mean")
    write_filtered_tsv(
        para_dir / "paralog_family_rankings_essential_filtered.tsv",
        fres["kept"],
        FAMILY_FIELDS + ["filtered_rank_max", "filtered_rank_mean"],
    )
    gene_ko_fields, gene_ko_rows = read_tsv_rows(ko_dir / "single_gene_ko.tsv")
    kg = apply_essential_filter(
        gene_ko_rows, essential_genes, "gene_id", "ko_rank"
    )
    write_filtered_tsv(
        ko_dir / "single_gene_ko_essential_filtered.tsv",
        kg["kept"],
        gene_ko_fields,
    )
    fam_ko_fields, fam_ko_rows = read_tsv_rows(ko_dir / "family_ko.tsv")
    kf = apply_essential_filter(
        fam_ko_rows, essential_fams, "family_id", "ko_rank"
    )
    write_filtered_tsv(
        ko_dir / "family_ko_essential_filtered.tsv",
        kf["kept"],
        fam_ko_fields,
    )

    growth_fields, growth_rows = read_tsv_rows(essential_tsv)
    gg = apply_essential_filter(
        [r for r in growth_rows if r.get("status") == "optimal"],
        essential_genes, "gene_id", "ko_rank",
    )
    write_filtered_tsv(
        ko_dir / "single_gene_ko_growth_essential_filtered.tsv",
        gg["kept"],
        growth_fields if record["mode"] != "explicit_id_list" else KO_GENE_FIELDS,
    )
    summary = {
        "essentiality_source": record["mode"],
        "essentiality_input": str(essential_tsv),
        "threshold": record["threshold"],
        "essential_genes": len(essential_genes),
        "essential_families": len(essential_fams),
        "fseof_genes_kept": len(gres["kept"]),
        "families_kept": len(fres["kept"]),
        "gene_ko_kept": len(kg["kept"]),
        "growth_ko_kept": len(gg["kept"]),
    }
    write_filter_manifest_json(
        ko_dir / "essential_filter_detail.json",
        {
            **summary,
            "rule": (
                "growth-objective single-gene KO loss >= threshold "
                "(optimal rows only, boundary inclusive); family "
                "essential = any member; product-objective KO tables "
                "are never an essentiality source"
            ),
        },
    )
    manifest.add(
        f"step9b_essential_filter[{od}]",
        inputs={
            "essential_genes_tsv": _record_file(essential_tsv),
            "fseof_genes": _record_file(fseof_dir / "all_genes_fseof.tsv"),
            "family_rankings": _record_file(
                para_dir / "paralog_family_rankings.tsv"
            ),
            "single_gene_ko": _record_file(ko_dir / "single_gene_ko.tsv"),
            "family_ko": _record_file(ko_dir / "family_ko.tsv"),
        },
        outputs={
            "fseof_gene_filtered": _record_file(
                fseof_dir / "all_genes_fseof_essential_filtered.tsv"
            ),
            "family_ranking_filtered": _record_file(
                para_dir / "paralog_family_rankings_essential_filtered.tsv"
            ),
            "single_gene_ko_filtered": _record_file(
                ko_dir / "single_gene_ko_essential_filtered.tsv"
            ),
            "family_ko_filtered": _record_file(
                ko_dir / "family_ko_essential_filtered.tsv"
            ),
            "detail": _record_file(ko_dir / "essential_filter_detail.json"),
            "growth_ko_filtered": _record_file(
                ko_dir / "single_gene_ko_growth_essential_filtered.tsv"
            ),
        },
        reused=False,
        elapsed=time.time() - t0e,
        detail=summary,
    )
    return summary


def _attach_solver(model: cobra.Model, name: str) -> None:
    if name == "highs":
        _register_highs()
    model.solver = name
    model.solver.configuration.tolerances.feasibility = 1e-9
    if name == "highs":
        model.solver.configuration.continuous_lp_backend = "linprog"


def _run_gapfill(
    *,
    name: str,
    refname: str,
    runtime_root: Path,
    prefer_runtime: bool = False,
    working: Path,
    cfg: Mapping[str, Any],
    cdir: Path,
    solver: str,
    extra_bin_dirs: Sequence[Path] = (),
    out_dir: Path,
    tag: str,
) -> Path:
    """Run AlphaGEM ``gapfilling.gapfill`` in-process under HiGHS."""
    _repo_imports(runtime_root, prefer_runtime=prefer_runtime)
    _use_default_solver(solver)

    def res(key):
        return _resolve(cdir, cfg.get(key))

    with _chdir(runtime_root):
        from src import gapfilling as gapfill_module

        if cfg.get("audit_checkpoint_dir"):
            _ckpt_dir = res("audit_checkpoint_dir")
            _ckpt_dir.mkdir(parents=True, exist_ok=True)
            _orig_prune = gapfill_module._prune_unsupported_reactions

            def _dumping_prune(model_obj, reaction_ids, **kw):
                _xml_path = _ckpt_dir / f"{name}_model2_pre_audit.xml"
                cobra.io.write_sbml_model(model_obj, str(_xml_path))
                _sha = sha256(_xml_path)
                _ur_path = _ckpt_dir / f"{name}_unsupported_reactions.json"
                save_json(
                    _ur_path,
                    {
                        "unsupported_reactions": list(reaction_ids),
                        "model2_sha256": _sha,
                        "model2_path": str(_xml_path),
                    },
                )
                print(f"[checkpoint] model2 -> {_xml_path} sha256={_sha}", flush=True)
                return _orig_prune(model_obj, reaction_ids, **kw)

            gapfill_module._prune_unsupported_reactions = _dumping_prune

        try:
            with _tool_path_scope(extra_bin_dirs=extra_bin_dirs):
                gapfill_module.gapfill(
                    name,
                    refname,
                    f"{refname}-GEM.xml",
                    grothmedium=cfg.get("grothmedium", "min"),
                    search_backend=cfg.get("gapfill_search_backend", "tblastn_genome"),
                    gapfill_genome=(
                        str(res("gapfill_genome")) if res("gapfill_genome") else None
                    ),
                    gapfill_coordinate_map=(
                        str(res("gapfill_coordinate_map"))
                        if res("gapfill_coordinate_map")
                        else None
                    ),
                    gapfill_coordinate_segments=(
                        str(res("gapfill_coordinate_segments"))
                        if res("gapfill_coordinate_segments")
                        else None
                    ),
                    gapfill_coordinate_contract=(
                        str(res("gapfill_coordinate_contract"))
                        if res("gapfill_coordinate_contract")
                        else None
                    ),
                    tie_break=cfg.get("gapfill_tie_break", "deterministic"),
                    unmapped_hit_policy=cfg.get("gapfill_unmapped_hit_policy", "skip_with_audit"),
                    bitscore_threshold=float(cfg.get("gapfill_bitscore_threshold", 200.0)),
                    bitscore_operator=cfg.get("gapfill_bitscore_operator", "gt"),
                    strategy=cfg.get("gapfill_strategy", "selected_medium"),
                    require_solver_tolerances=bool(cfg.get("require_solver_tolerances", True)),
                    require_optimal_solutions=bool(cfg.get("require_optimal_solutions", True)),
                    require_essential_reactions=bool(
                        cfg.get("require_essential_reactions", False)
                    ),
                )
        finally:
            if cfg.get("audit_checkpoint_dir"):
                gapfill_module._prune_unsupported_reactions = _orig_prune
    produced_xml = working / f"{name}-GEM.xml"
    produced_json = working / f"{name}-GEM.json"
    produced_audit = working / "gapfill_audit.json"
    for produced in (produced_xml, produced_json, produced_audit):
        if not produced.is_file():
            raise PipelineError(f"gapfill did not produce {produced}")
    final_xml = out_dir / f"{name}-GEM_loc_{tag}_gapfill.xml"
    shutil.copy2(produced_xml, final_xml)
    shutil.copy2(produced_json, out_dir / f"{name}-GEM_loc_{tag}_gapfill.json")
    shutil.copy2(produced_audit, out_dir / "gapfill_audit.json")
    return final_xml
