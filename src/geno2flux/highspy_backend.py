"""Native-highspy audit backend for the unified pipeline.

The vendored ``optlang_highs_interface`` routes every MILP through
``scipy.optimize.milp``. On large fungal MILPs, that ``scipy.optimize.milp``
wrapper can return false ``infeasible`` verdicts on feasible problems, while
direct native ``highspy`` calls solve them reliably. This module packages the
solver-layer replacement as a pipeline configuration option:

    "audit_solver_backend": "highspy",   # DEFAULT
    "audit_solver_threads": 1,
    "audit_solver_mip_feasibility_tolerance": 1e-9

``"highspy"`` is the DEFAULT: the pipeline installs this backend unless
the config explicitly opts out with ``"scipy"`` (which keeps the vendored
scipy.milp call layer). When the highspy package is missing the install
fails loudly — there is no silent fallback.

Only the ``Problem.optimize`` entry point is replaced — the audit logic
(``gapfilling._prune_unsupported_reactions``) is not modified.
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import Any

_DEFAULT_THREADS = 1
_DEFAULT_TOLERANCE = 1e-9

_installed_receipt: dict[str, Any] | None = None
_original_optimize: Any = None


def _load_interface_module():
    """Import the vendored optlang_highs_interface on sys.path."""
    import sys
    from pathlib import Path

    bundle_src = (
        Path(__file__).resolve().parents[2] / "third_party" / "alphagem" / "src"
    )
    if str(bundle_src) not in sys.path:
        sys.path.insert(0, str(bundle_src))
    import optlang_highs_interface

    return optlang_highs_interface


def is_installed() -> bool:
    """True once the native backend has been installed in this process."""
    return _installed_receipt is not None


def install_receipt() -> dict[str, Any] | None:
    """Copy of the receipt describing the installed backend, if any."""
    return dict(_installed_receipt) if _installed_receipt else None


def install_highspy_backend(
    *,
    threads: int = _DEFAULT_THREADS,
    mip_feasibility_tolerance: float = _DEFAULT_TOLERANCE,
) -> dict[str, Any]:
    """Patch ``optlang_highs_interface.Problem.optimize`` to native highspy.

    This replacement is idempotent and raises RuntimeError when the highspy
    package is missing — there is deliberately NO silent fallback to the
    vendored scipy route.
    """
    global _installed_receipt

    try:
        import highspy
    except ImportError as exc:  # no silent fallback to the scipy route
        raise RuntimeError(
            "audit_solver_backend='highspy' requires the highspy package "
            "(pip install highspy==1.15.1); refusing to silently fall "
            "back to the vendored scipy.milp route"
        ) from exc
    import numpy as np
    import scipy.sparse as sp
    # Validate before publishing the process-global patch: a failed install
    # must leave both the callable and its receipt untouched.
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 0:
        raise ValueError("highspy threads must be a nonnegative integer")
    mip_feasibility_tolerance = float(mip_feasibility_tolerance)
    if not np.isfinite(mip_feasibility_tolerance) or mip_feasibility_tolerance < 1e-10:
        raise ValueError("highspy MIP feasibility tolerance must be finite and >=1e-10")
    from optlang import interface as optlang_interface

    module = _load_interface_module()
    inf = highspy.kHighsInf

    def highs_native_optimize(
        self,
        *,
        verbosity=0,
        timeout=None,
        presolve=True,
        mip_rel_gap=None,
        node_limit=None,
        continuous_lp_backend="milp",
        feasibility_tolerance=None,
    ):
        # Never expose a previous solve after a rejected option/model/run.
        self._status = None
        self._var_primals = None
        self._f = None
        self.last_result_receipt = None
        if continuous_lp_backend not in {"milp", "linprog"}:
            raise ValueError("Unknown continuous LP backend")
        var_names = list(self._variables)
        n_v = len(var_names)
        obj_vec = np.array(
            [self._objective.get(n, 0.0) for n in var_names], dtype=float
        )
        if self._direction == "max":
            obj_vec = -obj_vec
        lb = np.array(
            [
                self.bounds[n][0] if self.bounds[n][0] is not None else -inf
                for n in var_names
            ]
        )
        ub = np.array(
            [
                self.bounds[n][1] if self.bounds[n][1] is not None else inf
                for n in var_names
            ]
        )
        integ = np.array(
            [
                0 if self._variable_types[n] == "continuous" else 1
                for n in var_names
            ],
            dtype=np.int32,
        )
        starts, indices, values = [0], [], []
        rl, ru = [], []
        for row in self._constraints.values():
            for vn, coefficient in row["coefficients"].items():
                indices.append(self._variables[vn])
                values.append(float(coefficient))
            starts.append(len(indices))
            rl.append(row["lb"] if row["lb"] is not None else -inf)
            ru.append(row["ub"] if row["ub"] is not None else inf)
        # Preserve sparsity throughout assembly; a dense genome-scale
        # constraint matrix can require gigabytes on every knockout solve.
        A = sp.csr_matrix((values, indices, starts), shape=(len(rl), n_v))
        h = highspy.Highs()

        def _set_option(name: str, value: Any) -> None:
            status = h.setOptionValue(name, value)
            # HighsStatus.kError means the option name/value was rejected;
            # failing loudly beats silently solving with different settings
            # than the caller requested.
            status_text = str(status)
            if "error" in status_text.lower():
                raise RuntimeError(
                    f"highspy rejected option {name}={value!r} "
                    f"(status={status_text})"
                )

        _set_option("threads", threads)
        _set_option("mip_feasibility_tolerance", mip_feasibility_tolerance)
        _set_option("output_flag", bool(verbosity))
        # Forward the same solve options the vendored scipy route
        # forwards, so optlang's auto-retry (presolve toggle) and
        # user-configured limits keep their meaning. Native HiGHS names
        # the node limit `mip_max_nodes` (scipy.milp calls it node_limit).
        forwarded: dict[str, Any] = {}
        if timeout is not None:
            _set_option("time_limit", float(timeout))
            forwarded["time_limit"] = float(timeout)
        if presolve is not None:
            presolve_value = "on" if presolve else "off"
            _set_option("presolve", presolve_value)
            forwarded["presolve"] = presolve_value
        if mip_rel_gap is not None:
            _set_option("mip_rel_gap", float(mip_rel_gap))
            forwarded["mip_rel_gap"] = float(mip_rel_gap)
        if node_limit is not None:
            _set_option("mip_max_nodes", int(node_limit))
            forwarded["mip_max_nodes"] = int(node_limit)
        if feasibility_tolerance is not None:
            tolerance = float(feasibility_tolerance)
            if not np.isfinite(tolerance) or tolerance < 1e-10:
                raise RuntimeError(
                    "highspy feasibility tolerance must be finite and "
                    f">=1e-10, got {feasibility_tolerance!r}"
                )
            _set_option("primal_feasibility_tolerance", tolerance)
            _set_option("dual_feasibility_tolerance", tolerance)
            forwarded["feasibility_tolerance"] = tolerance
        lp = highspy.HighsLp()
        lp.num_col_ = n_v
        lp.num_row_ = len(rl)
        lp.col_cost_ = obj_vec
        lp.col_lower_ = lb
        lp.col_upper_ = ub
        lp.row_lower_ = np.asarray(rl, dtype=float)
        lp.row_upper_ = np.asarray(ru, dtype=float)
        lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
        lp.a_matrix_.start_ = A.indptr.astype(np.int32)
        lp.a_matrix_.index_ = A.indices.astype(np.int32)
        lp.a_matrix_.value_ = A.data
        lp.integrality_ = [
            highspy.HighsVarType.kInteger if integ[i]
            else highspy.HighsVarType.kContinuous
            for i in range(n_v)
        ]
        model_status = h.passModel(lp)
        if "error" in str(model_status).lower():
            raise RuntimeError(f"highspy passModel failed (status={model_status})")
        run_status = h.run()
        if "error" in str(run_status).lower():
            raise RuntimeError(f"highspy run failed (status={run_status})")
        status_raw = str(h.getModelStatus()).lower()
        if "unboundedorinfeasible" in status_raw:
            # HiGHS kUnboundedOrInfeasible is AMBIGUOUS after presolve;
            # claiming infeasible would repeat the false-verdict class
            # this backend exists to fix. Report undecided instead.
            self._status = optlang_interface.UNDEFINED
        elif "infeasible" in status_raw:
            self._status = optlang_interface.INFEASIBLE
        elif "optimal" in status_raw:
            self._status = optlang_interface.OPTIMAL
        elif "unbounded" in status_raw:
            self._status = optlang_interface.UNBOUNDED
        elif "timelimit" in status_raw or "time_limit" in status_raw:
            self._status = optlang_interface.TIME_LIMIT
        elif "iterationlimit" in status_raw or "iteration_limit" in status_raw:
            self._status = optlang_interface.ITERATION_LIMIT
        else:
            self._status = optlang_interface.UNDEFINED
        sol = h.getSolution()
        primals = np.asarray(sol.col_value, dtype=float)
        primals_valid = bool(
            sol.value_valid
            and primals.shape == (n_v,)
            and np.isfinite(primals).all()
        )
        objective = (
            float(h.getObjectiveValue())
            if self._status == optlang_interface.OPTIMAL else None
        )
        if self._status == optlang_interface.OPTIMAL and (
            not primals_valid or not np.isfinite(objective)
        ):
            self._status = optlang_interface.UNDEFINED
            raise RuntimeError("highspy returned an invalid optimal solution")
        self._var_primals = primals if primals_valid else np.zeros(n_v)
        self._f = objective
        self.last_result_receipt = {
            "backend": "highspy.Highs native",
            "threads": threads,
            "mip_feasibility_tolerance": mip_feasibility_tolerance,
            "forwarded_options": dict(forwarded),
            "raw_status": status_raw,
            "status": self._status,
            "objective_value": self._f,
            "primals_from_backend": primals_valid,
            "variable_count": n_v,
            "constraint_count": len(rl),
        }
        return self._status

    global _original_optimize
    if _original_optimize is None:
        _original_optimize = module.Problem.optimize
    module.Problem.optimize = highs_native_optimize
    _installed_receipt = {
        "backend": "highspy.Highs native",
        "patch_target": "optlang_highs_interface.Problem.optimize",
        "threads": int(threads),
        "mip_feasibility_tolerance": float(mip_feasibility_tolerance),
        "audit_function": "original gapfilling._prune_unsupported_reactions (unmodified)",
    }
    return dict(_installed_receipt)


def uninstall_highspy_backend() -> bool:
    """Restore the vendored ``Problem.optimize`` if the patch is active.

    Returns True when a patch was removed. The pipeline installs the
    backend for its whole run on purpose; this entry point exists for
    tests and scoped use so the process-global patch never leaks into
    code that expects the default scipy route.
    """
    global _installed_receipt, _original_optimize
    if _installed_receipt is None:
        return False
    module = _load_interface_module()
    if _original_optimize is not None:
        module.Problem.optimize = _original_optimize
    _installed_receipt = None
    _original_optimize = None
    return True


@contextmanager
def scoped_highspy_backend():
    """Start from scipy and restore the caller's exact backend on exit.

    Each scope may install its own native settings or keep scipy. Nested
    scopes and failures restore the enclosing callable and install state.
    """
    global _installed_receipt, _original_optimize
    module = _load_interface_module()
    previous_optimize = module.Problem.optimize
    previous_receipt = _installed_receipt
    previous_original = _original_optimize
    uninstall_highspy_backend()
    try:
        yield
    finally:
        module.Problem.optimize = previous_optimize
        _installed_receipt = previous_receipt
        _original_optimize = previous_original
