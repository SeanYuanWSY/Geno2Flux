"""Unit tests for the optional native-highspy audit backend.

Covers:
1. Module import does NOT patch anything by default.
2. ``install_highspy_backend`` replaces ``Problem.optimize`` and records
   a receipt with the requested threads/tolerance.
3. The patched ``optimize`` forwards deterministic options to highspy
   and maps statuses (optimal / infeasible / unbounded) to optlang
   constants.
4. Pipeline wiring: native HiGHS is the default and scipy is an explicit
   opt-out; nested scopes restore the caller's backend on success/failure.
"""

from __future__ import annotations

import json
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from geno2flux import highspy_backend, pipeline


# ---------------------------------------------------------------------------
# Helpers: fake optlang problem + fake highspy
# ---------------------------------------------------------------------------


class _FakeProblem:
    """Minimal stand-in exposing the internals the patch reads."""

    def __init__(self, status="HighsModelStatus.kOptimal", objective=7.0):
        self._variables = {"x": 0, "k": 1}
        self._objective = {"x": 1.0, "k": 1.0}
        self._direction = "min"
        self.bounds = {"x": (0.0, 5.0), "k": (0.0, 1.0)}
        self._variable_types = {"x": "continuous", "k": "binary"}
        self._constraints = {
            "c0": {"coefficients": {"x": 1.0}, "lb": 1.0, "ub": None},
        }
        self._status = None
        self._var_primals = None
        self._f = None
        self.last_result_receipt = None
        self._highs_status = status
        self._highs_objective = objective


def _make_fake_highspy(fake_problem):
    """Build a fake ``highspy`` module wired to the fake problem."""
    mod = MagicMock()
    mod.kHighsInf = float("inf")

    sol = MagicMock()
    sol.col_value = [2.0, 1.0]
    sol.value_valid = True
    h = MagicMock()
    h.getSolution.return_value = sol
    h.getModelStatus.return_value = fake_problem._highs_status
    h.getObjectiveValue.return_value = fake_problem._highs_objective
    mod.Highs.return_value = h
    mod.HighsLp = MagicMock()
    mod.MatrixFormat.kRowwise = MagicMock()
    mod.HighsVarType.kInteger = MagicMock()
    mod.HighsVarType.kContinuous = MagicMock()
    return mod, h

@pytest.fixture(autouse=True)
def _restore_backend_after_test():
    """Never leak the process-global patch into other suite tests."""
    yield
    highspy_backend.uninstall_highspy_backend()
    iface = highspy_backend._load_interface_module()
    assert highspy_backend.is_installed() is False


class _OptlangConstants:
    OPTIMAL = "optimal"
    INFEASIBLE = "infeasible"
    UNBOUNDED = "unbounded"
    UNDEFINED = "undefined"
    TIME_LIMIT = "time_limit"
    ITERATION_LIMIT = "iteration_limit"



def test_uninstall_restores_original_and_is_idempotent():
    iface = highspy_backend._load_interface_module()
    original = iface.Problem.optimize
    fake_hs = MagicMock()
    fake_hs.kHighsInf = float("inf")
    with patch.dict("sys.modules", {"highspy": fake_hs}):
        highspy_backend.install_highspy_backend()
    assert iface.Problem.optimize is not original
    assert highspy_backend.uninstall_highspy_backend() is True
    assert iface.Problem.optimize is original
    assert highspy_backend.is_installed() is False
    # second uninstall is a no-op returning False
    assert highspy_backend.uninstall_highspy_backend() is False
    assert iface.Problem.optimize is original


def test_patched_optimize_forwards_timeout_presolve_and_limits():
    problem = _FakeProblem(status="HighsModelStatus.kTimeLimit")
    fake_hs, highs = _make_fake_highspy(problem)
    fake_optlang = MagicMock()
    fake_optlang.interface = _OptlangConstants

    with patch.dict(
        "sys.modules",
        {"highspy": fake_hs, "optlang": fake_optlang},
    ):
        highspy_backend.install_highspy_backend()
        iface = highspy_backend._load_interface_module()
        result = iface.Problem.optimize(
            problem,
            timeout=12.5,
            presolve=False,
            mip_rel_gap=0.01,
            node_limit=7,
        )

    assert result == "time_limit"
    highs.setOptionValue.assert_any_call("time_limit", 12.5)
    highs.setOptionValue.assert_any_call("presolve", "off")
    highs.setOptionValue.assert_any_call("mip_rel_gap", 0.01)
    highs.setOptionValue.assert_any_call("mip_max_nodes", 7)
    receipt = problem.last_result_receipt
    assert receipt["forwarded_options"] == {
        "time_limit": 12.5,
        "presolve": "off",
        "mip_rel_gap": 0.01,
        "mip_max_nodes": 7,
    }

def test_patched_optimize_ambiguous_unbounded_or_infeasible_is_undefined():
    """kUnboundedOrInfeasible must stay undecided, never 'infeasible'."""
    problem = _FakeProblem(status="HighsModelStatus.kUnboundedOrInfeasible")
    fake_hs, _highs = _make_fake_highspy(problem)
    fake_optlang = MagicMock()
    fake_optlang.interface = _OptlangConstants
    with patch.dict(
        "sys.modules", {"highspy": fake_hs, "optlang": fake_optlang}
    ):
        highspy_backend.install_highspy_backend()
        iface = highspy_backend._load_interface_module()
        result = iface.Problem.optimize(problem)
    assert result == "undefined"
    assert problem._status == "undefined"
    assert problem._f is None


def test_patched_optimize_forwards_feasibility_tolerance():
    problem = _FakeProblem(status="HighsModelStatus.kOptimal")
    fake_hs, highs = _make_fake_highspy(problem)
    fake_optlang = MagicMock()
    fake_optlang.interface = _OptlangConstants
    with patch.dict(
        "sys.modules", {"highspy": fake_hs, "optlang": fake_optlang}
    ):
        highspy_backend.install_highspy_backend()
        iface = highspy_backend._load_interface_module()
        result = iface.Problem.optimize(problem, feasibility_tolerance=1e-9)
    assert result == "optimal"
    highs.setOptionValue.assert_any_call("primal_feasibility_tolerance", 1e-9)
    highs.setOptionValue.assert_any_call("dual_feasibility_tolerance", 1e-9)
    assert problem.last_result_receipt["forwarded_options"] == {
        "presolve": "on",
        "feasibility_tolerance": 1e-9,
    }


def test_patched_optimize_rejects_invalid_feasibility_tolerance():
    problem = _FakeProblem(status="HighsModelStatus.kOptimal")
    fake_hs, _highs = _make_fake_highspy(problem)
    fake_optlang = MagicMock()
    fake_optlang.interface = _OptlangConstants
    with patch.dict(
        "sys.modules", {"highspy": fake_hs, "optlang": fake_optlang}
    ):
        highspy_backend.install_highspy_backend()
        iface = highspy_backend._load_interface_module()
        with pytest.raises(RuntimeError, match="feasibility tolerance"):
            iface.Problem.optimize(problem, feasibility_tolerance=1e-12)



def test_rejected_option_raises_instead_of_silent_ignore():
    problem = _FakeProblem(status="HighsModelStatus.kOptimal")
    fake_hs, highs = _make_fake_highspy(problem)
    fake_optlang = MagicMock()
    fake_optlang.interface = _OptlangConstants
    highs.setOptionValue.return_value = "HighsStatus.kError"
    with patch.dict(
        "sys.modules", {"highspy": fake_hs, "optlang": fake_optlang}
    ):
        highspy_backend.install_highspy_backend()
        iface = highspy_backend._load_interface_module()
        with pytest.raises(RuntimeError, match="rejected option"):
            iface.Problem.optimize(problem)


def test_patched_optimize_presolve_on_forwarded():
    problem = _FakeProblem(status="HighsModelStatus.kOptimal")
    fake_hs, highs = _make_fake_highspy(problem)
    fake_optlang = MagicMock()
    fake_optlang.interface = _OptlangConstants
    with patch.dict(
        "sys.modules", {"highspy": fake_hs, "optlang": fake_optlang}
    ):
        highspy_backend.install_highspy_backend()
        iface = highspy_backend._load_interface_module()
        iface.Problem.optimize(problem, presolve=True)
    highs.setOptionValue.assert_any_call("presolve", "on")
    assert problem.last_result_receipt["forwarded_options"] == {
        "presolve": "on"
    }

# ---------------------------------------------------------------------------
# 1-2: install semantics
# ---------------------------------------------------------------------------


def test_import_alone_does_not_patch():
    iface = highspy_backend._load_interface_module()
    original = iface.Problem.optimize
    assert highspy_backend.is_installed() is False
    assert highspy_backend.install_receipt() is None
    assert iface.Problem.optimize is original


def test_install_patches_and_records_receipt():
    fake_hs = MagicMock()
    fake_hs.kHighsInf = float("inf")
    with patch.dict("sys.modules", {"highspy": fake_hs}):
        receipt = highspy_backend.install_highspy_backend(
            threads=1, mip_feasibility_tolerance=1e-9
        )
    iface = highspy_backend._load_interface_module()
    assert highspy_backend.is_installed() is True
    assert receipt["backend"] == "highspy.Highs native"
    assert receipt["threads"] == 1
    assert receipt["mip_feasibility_tolerance"] == 1e-9
    assert highspy_backend.install_receipt() == receipt
    assert callable(iface.Problem.optimize)


# ---------------------------------------------------------------------------
# 3: patched optimize behavior (fake highspy)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "status,expected",
    [
        ("HighsModelStatus.kOptimal", "optimal"),
        ("HighsModelStatus.kInfeasible", "infeasible"),
        ("HighsModelStatus.kUnbounded", "unbounded"),
    ],
)
def test_patched_optimize_maps_status_and_options(status, expected):
    problem = _FakeProblem(status=status)
    fake_hs, highs = _make_fake_highspy(problem)
    fake_optlang = MagicMock()
    fake_optlang.interface = _OptlangConstants

    with patch.dict(
        "sys.modules",
        {"highspy": fake_hs, "optlang": fake_optlang},
    ):
        highspy_backend.install_highspy_backend(
            threads=1, mip_feasibility_tolerance=1e-9
        )
        iface = highspy_backend._load_interface_module()
        result = iface.Problem.optimize(problem)

    assert result == expected
    assert problem._status == expected
    highs.setOptionValue.assert_any_call("threads", 1)
    highs.setOptionValue.assert_any_call("mip_feasibility_tolerance", 1e-9)
    highs.setOptionValue.assert_any_call("output_flag", False)
    receipt = problem.last_result_receipt
    assert receipt["backend"] == "highspy.Highs native"
    assert receipt["raw_status"] == status.lower()
    if expected == "optimal":
        assert problem._f == 7.0
        assert receipt["objective_value"] == 7.0
    else:
        assert problem._f is None


# ---------------------------------------------------------------------------
# 4: pipeline wiring
# ---------------------------------------------------------------------------


def _write_config(tmp_path: Path, extra: dict, stem: str = "config") -> Path:
    inputs = tmp_path / "inputs"
    inputs.mkdir(exist_ok=True)
    for name in (
        "proteome.fasta",
        "withgaps.yml",
        "loc.tsv",
        "ref.fasta",
        "ref.xlsx",
    ):
        (inputs / name).write_text("x")
    cfg = {
        "name": "test",
        "refname": "yeast",
        "solver": "highs",
        "allow_unverified_reuse": True,
        "proteome_fasta": str(inputs / "proteome.fasta"),
        "withgaps_model": str(inputs / "withgaps.yml"),
        "localization_tsv": str(inputs / "loc.tsv"),
        "reference_fasta": str(inputs / "ref.fasta"),
        "reference_xlsx": str(inputs / "ref.xlsx"),
    }
    cfg.update(extra)
    config_path = tmp_path / f"{stem}.json"
    config_path.write_text(json.dumps(cfg))
    return config_path


def _light_pipeline_patches(install_mock):
    """Patch everything heavy; the wiring under test runs before all of it."""
    return [
        patch.object(pipeline, "_use_default_solver", MagicMock()),
        patch.object(
            pipeline, "_prepare_alphagem_runtime", return_value=None
        ),
        patch.object(pipeline, "_ensure_data_available", return_value=None),
        patch(
            "geno2flux.highspy_backend.install_highspy_backend",
            install_mock,
        ),
    ]


def test_pipeline_config_installs_backend(tmp_path: Path):
    config_path = _write_config(
        tmp_path,
        {
            "audit_solver_backend": "highspy",
            "audit_solver_threads": 2,
            "audit_solver_mip_feasibility_tolerance": 5e-9,
        },
    )
    install_mock = MagicMock()
    with ExitStack() as stack:
        for p in _light_pipeline_patches(install_mock):
            stack.enter_context(p)
        try:
            pipeline.run_pipeline(
                config_path,
                tmp_path / "out",
                localization_order_override="pre-gapfill",
            )
        except Exception:
            # wiring happens before any model work; later failures with
            # placeholder inputs are irrelevant to this assertion.
            pass
    assert install_mock.called
    kwargs = install_mock.call_args.kwargs
    assert kwargs["threads"] == 2
    assert kwargs["mip_feasibility_tolerance"] == 5e-9


def test_pipeline_invalid_backend_rejected(tmp_path: Path):
    config_path = _write_config(
        tmp_path, {"audit_solver_backend": "gurobi-native"}
    )
    with ExitStack() as stack:
        for p in _light_pipeline_patches(MagicMock()):
            stack.enter_context(p)
        with pytest.raises(pipeline.PipelineError, match="audit_solver_backend"):
            pipeline.run_pipeline(
                config_path,
                tmp_path / "out",
                localization_order_override="pre-gapfill",
            )


def test_pipeline_default_installs_highspy(tmp_path: Path):
    """Zero-config (new default) routes through install_highspy_backend."""
    config_path = _write_config(tmp_path, {})
    install_mock = MagicMock()
    with ExitStack() as stack:
        for p in _light_pipeline_patches(install_mock):
            stack.enter_context(p)
        try:
            pipeline.run_pipeline(
                config_path,
                tmp_path / "out",
                localization_order_override="pre-gapfill",
            )
        except Exception:
            pass
    assert install_mock.called
    kwargs = install_mock.call_args.kwargs
    assert kwargs == {"threads": 1, "mip_feasibility_tolerance": 1e-9}


def test_pipeline_explicit_scipy_opts_out(tmp_path: Path):
    config_path = _write_config(
        tmp_path, {"audit_solver_backend": "scipy"}
    )
    install_mock = MagicMock()
    with ExitStack() as stack:
        for p in _light_pipeline_patches(install_mock):
            stack.enter_context(p)
        try:
            pipeline.run_pipeline(
                config_path,
                tmp_path / "out",
                localization_order_override="pre-gapfill",
            )
        except Exception:
            pass
    assert not install_mock.called


def test_install_missing_highspy_raises_no_silent_fallback():
    """highspy absent -> RuntimeError; never a silent scipy fallback."""
    with patch.dict("sys.modules", {"highspy": None}):
        with pytest.raises(RuntimeError, match="requires the highspy package"):
            highspy_backend.install_highspy_backend()
    assert highspy_backend.is_installed() is False

def test_pipeline_restores_solver_layer_across_runs(tmp_path: Path):
    """Backend installs per run and never leaks between runs (finally-restore)."""
    fake_hs = MagicMock()
    fake_hs.kHighsInf = float("inf")
    iface = highspy_backend._load_interface_module()
    original = iface.Problem.optimize
    high_cfg = _write_config(
        tmp_path, {"audit_solver_backend": "highspy"}, stem="config_high"
    )
    default_cfg = _write_config(tmp_path, {}, stem="config_default")
    install_calls = []
    real_install = highspy_backend.install_highspy_backend

    def _tracked_install(**kwargs):
        install_calls.append(kwargs)
        return real_install(**kwargs)

    with patch.dict(
        "sys.modules", {"highspy": fake_hs}
    ), patch.object(
        pipeline, "_use_default_solver", MagicMock()
    ), patch.object(
        pipeline, "_prepare_alphagem_runtime", return_value=None
    ), patch.object(
        pipeline, "_ensure_data_available", return_value=None
    ), patch.object(
        pipeline, "_check_sha", return_value=None
    ), patch.object(
        highspy_backend, "install_highspy_backend", _tracked_install
    ):
        try:
            pipeline.run_pipeline(
                high_cfg, tmp_path / "out1",
                localization_order_override="pre-gapfill",
            )
        except Exception:
            pass
        # install happened during the run ...
        assert install_calls
        # ... and the finally-uninstall already restored the layer.
        assert highspy_backend.is_installed() is False
        assert iface.Problem.optimize is original
        try:
            pipeline.run_pipeline(
                default_cfg, tmp_path / "out2",
                localization_order_override="pre-gapfill",
            )
        except Exception:
            pass
        assert iface.Problem.optimize is original
        assert highspy_backend.is_installed() is False


@pytest.mark.parametrize("fail", [False, True])
def test_nested_backend_scopes_restore_exact_caller_state(fail):
    iface = highspy_backend._load_interface_module()
    scipy_optimize = iface.Problem.optimize
    fake_hs = MagicMock()
    fake_hs.kHighsInf = float("inf")
    with patch.dict("sys.modules", {"highspy": fake_hs}):
        caller_receipt = highspy_backend.install_highspy_backend(threads=2)
        caller_optimize = iface.Problem.optimize
        with highspy_backend.scoped_highspy_backend():
            assert iface.Problem.optimize is scipy_optimize
            assert highspy_backend.install_receipt() is None
            outer_receipt = highspy_backend.install_highspy_backend(threads=1)
            outer_optimize = iface.Problem.optimize
            try:
                with highspy_backend.scoped_highspy_backend():
                    assert iface.Problem.optimize is scipy_optimize
                    assert highspy_backend.is_installed() is False
                    highspy_backend.install_highspy_backend(threads=3)
                    if fail:
                        raise RuntimeError("inner failure")
            except RuntimeError as exc:
                assert fail and str(exc) == "inner failure"
            assert iface.Problem.optimize is outer_optimize
            assert highspy_backend.install_receipt() == outer_receipt
        assert iface.Problem.optimize is caller_optimize
        assert highspy_backend.install_receipt() == caller_receipt
        assert highspy_backend.uninstall_highspy_backend() is True
        assert iface.Problem.optimize is scipy_optimize


def test_backend_scope_failure_restores_uninstalled_callable():
    iface = highspy_backend._load_interface_module()
    original = iface.Problem.optimize
    fake_hs = MagicMock()
    fake_hs.kHighsInf = float("inf")
    with patch.dict("sys.modules", {"highspy": fake_hs}):
        with pytest.raises(RuntimeError, match="scope failure"):
            with highspy_backend.scoped_highspy_backend():
                highspy_backend.install_highspy_backend()
                raise RuntimeError("scope failure")
    assert iface.Problem.optimize is original
    assert highspy_backend.install_receipt() is None
    assert highspy_backend.uninstall_highspy_backend() is False


@pytest.mark.parametrize(
    "options",
    [
        {"threads": -1},
        {"threads": 1.5},
        {"threads": True},
        {"mip_feasibility_tolerance": float("nan")},
        {"mip_feasibility_tolerance": float("inf")},
        {"mip_feasibility_tolerance": 1e-12},
        {"mip_feasibility_tolerance": "invalid"},
    ],
)
def test_invalid_install_preserves_existing_backend(options):
    iface = highspy_backend._load_interface_module()
    original = iface.Problem.optimize
    fake_hs = MagicMock()
    fake_hs.kHighsInf = float("inf")
    with patch.dict("sys.modules", {"highspy": fake_hs}):
        receipt = highspy_backend.install_highspy_backend()
        installed = iface.Problem.optimize
        with pytest.raises(ValueError):
            highspy_backend.install_highspy_backend(**options)
    assert iface.Problem.optimize is installed
    assert highspy_backend.install_receipt() == receipt
    highspy_backend.uninstall_highspy_backend()
    assert iface.Problem.optimize is original


@pytest.mark.parametrize("operation", ["passModel", "run"])
def test_native_error_clears_previous_solution_and_raises(operation):
    problem = _FakeProblem()
    problem._status = "optimal"
    problem._f = 123.0
    problem._var_primals = [123.0, 0.0]
    problem.last_result_receipt = {"status": "optimal"}
    fake_hs, highs = _make_fake_highspy(problem)
    getattr(highs, operation).return_value = "HighsStatus.kError"
    with patch.dict("sys.modules", {"highspy": fake_hs}):
        highspy_backend.install_highspy_backend()
        iface = highspy_backend._load_interface_module()
        with pytest.raises(RuntimeError, match=f"highspy {operation} failed"):
            iface.Problem.optimize(problem)
    assert problem._status is None
    assert problem._f is None
    assert problem._var_primals is None
    assert problem.last_result_receipt is None
    if operation == "passModel":
        highs.run.assert_not_called()


@pytest.mark.parametrize("invalid", ["flag", "length", "primal_nan", "objective_nan"])
def test_invalid_optimal_solution_is_not_published(invalid):
    problem = _FakeProblem()
    fake_hs, highs = _make_fake_highspy(problem)
    sol = highs.getSolution.return_value
    if invalid == "flag":
        sol.value_valid = False
    elif invalid == "length":
        sol.col_value = [2.0]
    elif invalid == "primal_nan":
        sol.col_value = [float("nan"), 1.0]
    else:
        highs.getObjectiveValue.return_value = float("nan")
    with patch.dict("sys.modules", {"highspy": fake_hs}):
        highspy_backend.install_highspy_backend()
        iface = highspy_backend._load_interface_module()
        with pytest.raises(RuntimeError, match="invalid optimal solution"):
            iface.Problem.optimize(problem)
    assert problem._status == "undefined"
    assert problem._f is None
    assert problem._var_primals is None


def test_missing_nonoptimal_primals_are_explicit_standins():
    import numpy as np

    problem = _FakeProblem(status="HighsModelStatus.kInfeasible")
    fake_hs, highs = _make_fake_highspy(problem)
    highs.getSolution.return_value.value_valid = False
    highs.getSolution.return_value.col_value = np.array([99.0, 99.0])
    with patch.dict("sys.modules", {"highspy": fake_hs}):
        highspy_backend.install_highspy_backend()
        iface = highspy_backend._load_interface_module()
        assert iface.Problem.optimize(problem) == "infeasible"
    assert list(problem._var_primals) == [0.0, 0.0]
    assert problem.last_result_receipt["primals_from_backend"] is False
    assert problem._f is None


def test_constraint_assembly_preserves_sparse_rows_and_bounds():
    import numpy as np
    import scipy.sparse as sp

    problem = _FakeProblem()
    problem._constraints["c1"] = {
        "coefficients": {"k": -2.0}, "lb": None, "ub": 3.0,
    }
    fake_hs, highs = _make_fake_highspy(problem)
    with patch.dict("sys.modules", {"highspy": fake_hs}), patch.object(
        sp, "csr_matrix", wraps=sp.csr_matrix
    ) as sparse_constructor:
        highspy_backend.install_highspy_backend()
        iface = highspy_backend._load_interface_module()
        iface.Problem.optimize(problem, verbosity=1)
    # CSR arrays go directly to scipy, rather than a rows-by-columns array.
    assert isinstance(sparse_constructor.call_args.args[0], tuple)
    lp = highs.passModel.call_args.args[0]
    np.testing.assert_array_equal(lp.a_matrix_.start_, [0, 1, 2])
    np.testing.assert_array_equal(lp.a_matrix_.index_, [0, 1])
    np.testing.assert_array_equal(lp.a_matrix_.value_, [1.0, -2.0])
    np.testing.assert_array_equal(lp.row_lower_, [1.0, -np.inf])
    np.testing.assert_array_equal(lp.row_upper_, [np.inf, 3.0])
    highs.setOptionValue.assert_any_call("output_flag", True)


@pytest.mark.parametrize("direction,expected", [("max", 14.0), ("min", 7.0)])
def test_native_integer_model_objective_direction_and_offset(direction, expected):
    pytest.importorskip("highspy")
    iface = highspy_backend._load_interface_module()
    problem = iface.Problem()
    problem.add_variable("x")
    problem.add_variable("k", variable_type="binary")
    problem.set_variable_bounds("x", 0.0, 5.0)
    problem.set_variable_bounds("k", 0.0, 1.0)
    problem.add_constraint("range", {"x": 1.0, "k": 1.0}, lb=1.0, ub=4.0)
    problem.objective = {"x": 2.0, "k": 3.0}
    problem.direction = direction
    problem.offset = 5.0
    highspy_backend.install_highspy_backend()
    assert problem.optimize(feasibility_tolerance=1e-9) == "optimal"
    assert problem.objective_value == pytest.approx(expected)
    assert problem.get_var_primal("k") == pytest.approx(1.0 if direction == "max" else 0.0)
    assert problem.last_result_receipt["primals_from_backend"]
    decoded_receipt = json.loads(json.dumps(problem.last_result_receipt, allow_nan=False))
    assert decoded_receipt["primals_from_backend"] is True


def test_native_binary_variable_is_not_relaxed():
    pytest.importorskip("highspy")
    iface = highspy_backend._load_interface_module()
    problem = iface.Problem()
    problem.add_variable("k", variable_type="binary")
    problem.set_variable_bounds("k", 0.0, 1.0)
    problem.add_constraint("fractional_cap", {"k": 1.0}, ub=0.5)
    problem.objective = {"k": 1.0}
    problem.direction = "max"
    highspy_backend.install_highspy_backend()
    assert problem.optimize() == "optimal"
    assert problem.objective_value == pytest.approx(0.0)
    assert problem.get_var_primal("k") == pytest.approx(0.0)


@pytest.mark.parametrize("verdict", ["infeasible", "unbounded"])
def test_native_nonoptimal_verdict_has_no_objective(verdict):
    pytest.importorskip("highspy")
    iface = highspy_backend._load_interface_module()
    problem = iface.Problem()
    problem.add_variable("x")
    problem.objective = {"x": 1.0}
    if verdict == "infeasible":
        problem.set_variable_bounds("x", 0.0, 0.0)
        problem.add_constraint("impossible", {"x": 1.0}, lb=1.0)
    highspy_backend.install_highspy_backend()
    assert problem.optimize(presolve=False) == verdict
    assert problem.objective_value is None
    assert problem.last_result_receipt["status"] == verdict
