from __future__ import annotations

import argparse
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="geno2flux",
        description="Geno2Flux genome-to-GEM reconstruction and metabolic target prioritization",
    )
    s = p.add_subparsers(dest="command", required=True)

    sub = s.add_parser(
        "pipeline",
        help="Run the unified reconstruction and/or target-prioritization pipeline",
    )
    sub.add_argument("--config", type=Path, required=True)
    sub.add_argument("--output", type=Path, required=True)
    sub.add_argument("--solver", type=str)
    sub.add_argument(
        "--localization-order",
        type=str,
        choices=["pre-gapfill", "post-gapfill", "both", "skip"],
    )

    sub = s.add_parser(
        "lock-resources",
        help="Generate a SHA-256 manifest for an external resource directory",
    )
    sub.add_argument("--root", type=Path, required=True)
    sub.add_argument("--output", type=Path, required=True)

    s.add_parser(
        "predict-proteome",
        add_help=False,
        help="Predict or standardize genome -> proteome + AlphaGEM contracts (delegates to scripts/predict_proteome.py)",
    )

    if argv is None:
        import sys as _sys
        raw_argv = list(_sys.argv[1:])
    else:
        raw_argv = list(argv)

    if raw_argv and raw_argv[0] == "predict-proteome":
        import sys as _sys
        scripts_dir = Path(__file__).resolve().parents[2] / "scripts"
        if str(scripts_dir) not in _sys.path:
            _sys.path.insert(0, str(scripts_dir))
        import predict_proteome  # type: ignore[import-not-found]
        return predict_proteome.main(raw_argv[1:])

    a = p.parse_args(raw_argv)
    if a.command == "lock-resources":
        from .common import save_json, sha256

        if a.output.exists():
            raise FileExistsError(a.output)
        save_json(
            a.output,
            {
                str(f.relative_to(a.root)): sha256(f)
                for f in sorted(a.root.rglob("*"))
                if f.is_file()
                and ".git" not in f.parts
                and "__pycache__" not in f.parts
            },
        )
        return 0

    if a.command == "pipeline":
        from .pipeline import run_pipeline

        run_pipeline(
            a.config,
            a.output,
            solver_override=a.solver,
            localization_order_override=a.localization_order,
        )
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
