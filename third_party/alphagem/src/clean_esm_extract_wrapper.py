"""Run CLEAN's canonical ESM extractor with an explicit token budget."""

import os
import runpy
import sys


budget = os.environ.get("ALPHAGEM_CLEAN_ESM_TOKEN_BUDGET")
if budget is None or not budget.isdigit() or int(budget) <= 0:
    raise RuntimeError("Missing positive ALPHAGEM_CLEAN_ESM_TOKEN_BUDGET")
if "--toks_per_batch" in sys.argv:
    indices = [
        index
        for index, value in enumerate(sys.argv)
        if value == "--toks_per_batch"
    ]
    if (
        len(indices) != 1
        or indices[0] + 1 >= len(sys.argv)
        or sys.argv[indices[0] + 1] != budget
    ):
        raise RuntimeError(
            "CLEAN ESM token budget conflicts with the formal run contract"
        )
else:
    sys.argv.extend(["--toks_per_batch", budget])
extract_script = os.environ.get("ALPHAGEM_CLEAN_ESM_EXTRACT")
if not extract_script:
    raise RuntimeError("Missing ALPHAGEM_CLEAN_ESM_EXTRACT")
runpy.run_path(extract_script, run_name="__main__")
