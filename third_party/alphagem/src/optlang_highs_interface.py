"""optlang-compatible HiGHS interface with explicit MILP/continuous-LP APIs.

This module registers an optlang interface named ``highs`` into the
cobra/optlang solver lookup so ``--solver highs`` works end to end. It exists
because optlang 1.9.1 ships no HiGHS interface; the backend is SciPy's public
``milp`` wrapper of HiGHS. Developed against cobra 0.31.1 / scipy 1.15.3 and
checked line-by-line against cobra 0.29.1 (server): the GapFiller and solver
selection call sites are identical between those two cobra versions.

Honesty notes (kept deliberately explicit, per project contract):

- The public ``scipy.optimize.milp`` API exposes no thread parameter, so this
  interface never claims to set a thread count. HiGHS runs with its own
  internal defaults for anything scipy does not expose.
- Only ``time_limit`` (from ``Configuration.timeout``) and ``mip_rel_gap``
  (from ``tolerances.optimality``) are forwarded to the backend. The
  ``feasibility`` and ``integrality`` tolerances are stored and reported
  verbatim for audit read-back but are NOT forwarded on the milp route: scipy's milp options do
  not expose HiGHS's primal/integrality tolerance knobs, which therefore stay
  at HiGHS defaults. Every optimize() result records exactly which options
  were forwarded.
- Ordinary continuous LPs may explicitly select ``continuous_lp_backend="linprog"``.
  That route uses the public linprog(method="highs") API and forwards the
  requested primal/dual feasibility tolerance; integer problems retain milp.
  Historical/default calls remain on milp. This choice changes numerical
  solving, not the mathematical bounds, objective or constraints.
- Dual values (reduced costs, shadow prices) are not available through
  scipy.milp and raise NotImplementedError here. cobra skips duals for
  integer models and tolerates their absence for plain LPs.
- Non-decisive backend states (iteration/time limit, numeric trouble, "other")
  are mapped to the corresponding optlang statuses verbatim; nothing is ever
  reported as optimal unless HiGHS proved optimality (status 0).
"""

from collections import OrderedDict

import numpy as np
import optlang.interface as interface
import scipy
import scipy.optimize
import scipy.sparse
from optlang import symbolics
from optlang.expression_parsing import parse_optimization_expression
from optlang.util import inheritdocstring


# scipy.optimize.milp status codes (SciPy 1.14.1 and 1.15.x agree; the
# mip_node_count/mip_dual_bound/mip_gap result fields exist in both).
SCIPY_MILP_STATUS = {
    0: interface.OPTIMAL,
    2: interface.INFEASIBLE,
    3: interface.UNBOUNDED,
    4: interface.UNDEFINED,
}


class Problem(object):
    """Sparse problem state translated into explicit SciPy HiGHS calls."""

    def __init__(self):
        self._variables = OrderedDict()
        self._variable_types = {}
        self.bounds = OrderedDict()
        self._constraints = OrderedDict()
        self._objective = OrderedDict()
        self._offset = 0.0
        self._direction = "min"
        self._reset_solution()
        self.last_result_receipt = None

    def _reset_solution(self):
        self._status = None
        self._var_primals = None
        self._f = None

    def add_variable(self, name, variable_type="continuous"):
        if name in self._variables:
            raise ValueError(f"A variable named {name} already exists.")
        self._variables[name] = len(self._variables)
        self._variable_types[name] = variable_type
        self.bounds[name] = (None, None)
        self._reset_solution()

    def set_variable_type(self, name, variable_type):
        self._variable_types[name] = variable_type
        self._reset_solution()

    def set_variable_bounds(self, name, lower, upper):
        self.bounds[name] = (lower, upper)
        self._reset_solution()

    def add_constraint(self, name, coefficients, lb=None, ub=None):
        if name in self._constraints:
            raise ValueError(f"A constraint named {name} already exists.")
        self._constraints[name] = {
            "coefficients": dict(coefficients),
            "lb": lb,
            "ub": ub,
        }
        self._reset_solution()

    def set_constraint_bounds(self, name, lb, ub):
        self._constraints[name]["lb"] = lb
        self._constraints[name]["ub"] = ub
        self._reset_solution()

    def set_constraint_coefficients(self, name, coefficients):
        """Replace a stored constraint row (cobra updates rows in place)."""
        self._constraints[name]["coefficients"] = dict(coefficients)
        self._reset_solution()

    def change_var_name(self, name, value):
        if name == value:
            return
        if value in self._variables:
            raise ValueError("A variable with that name already exists")
        for mapping in (self._variables, self.bounds):
            mapping[value] = mapping.pop(name)
        self._variable_types[value] = self._variable_types.pop(name)
        self._objective = OrderedDict(
            (value if key == name else key, coef)
            for key, coef in self._objective.items()
        )
        for constraint in self._constraints.values():
            if name in constraint["coefficients"]:
                constraint["coefficients"][value] = constraint[
                    "coefficients"
                ].pop(name)
        # pop+reinsert moves the renamed variable to the end of the
        # OrderedDict, so re-enumerate stored indices in iteration order:
        # optimize() builds columns in iteration order and constraints index
        # by these stored values, and the two must stay the same mapping.
        self._variables = OrderedDict(
            (key, position) for position, key in enumerate(self._variables)
        )
        self._reset_solution()

    def change_constraint_name(self, name, value):
        if name != value and value in self._constraints:
            raise ValueError("A constraint with that name already exists")
        self._constraints[value] = self._constraints.pop(name)

    def remove_variable(self, name):
        del self._variables[name]
        del self._variable_types[name]
        del self.bounds[name]
        self._objective.pop(name, None)
        for constraint in list(self._constraints.values()):
            constraint["coefficients"].pop(name, None)
        self._variables = OrderedDict(
            (key, position) for position, key in enumerate(self._variables)
        )
        self._reset_solution()

    def remove_constraint(self, name):
        del self._constraints[name]
        self._reset_solution()

    def get_var_primal(self, name):
        if self._var_primals is None:
            return None
        return self._var_primals[self._variables[name]]

    def get_constraint_primal(self, name):
        if self._var_primals is None:
            return None
        coefficients = self._constraints[name]["coefficients"]
        return sum(
            coefficient * self._var_primals[self._variables[variable_name]]
            for variable_name, coefficient in coefficients.items()
        )

    @property
    def objective(self):
        return self._objective

    @objective.setter
    def objective(self, value):
        self._objective = value
        self._reset_solution()

    @property
    def offset(self):
        return self._offset

    @offset.setter
    def offset(self, value):
        self._offset = value
        self._reset_solution()

    @property
    def direction(self):
        return self._direction

    @direction.setter
    def direction(self, value):
        # cobra's Model.optimize() restores the objective direction after
        # solving; a no-change assignment must not invalidate the solution
        # that downstream code (GapFiller, pruning) still reads via primals.
        if value != self._direction:
            self._direction = value
            self._reset_solution()

    @property
    def status(self):
        if self._status is None:
            raise RuntimeError("Problem has not been optimized yet")
        return self._status

    @property
    def objective_value(self):
        # None when no optimal value exists (unsolved, infeasible, limited):
        # optlang/cobra read objective.value on non-optimal paths too.
        if self._f is None:
            return None
        value = self._f if self._direction == "min" else -self._f
        return value + self._offset

    def export_problem(self):
        """Audit copy of the full problem data in plain JSON-able types."""
        return {
            "variables": [
                {
                    "name": name,
                    "type": self._variable_types[name],
                    "lb": self.bounds[name][0],
                    "ub": self.bounds[name][1],
                }
                for name in self._variables
            ],
            "constraints": [
                {
                    "name": name,
                    "coefficients": dict(row["coefficients"]),
                    "lb": row["lb"],
                    "ub": row["ub"],
                }
                for name, row in self._constraints.items()
            ],
            "objective": {
                "coefficients": dict(self._objective),
                "offset": self._offset,
                "direction": self._direction,
            },
        }

    def optimize(
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
        """Solve through the selected public SciPy API and store its actual options."""
        variable_names = list(self._variables)
        objective_vector = np.array(
            [self._objective.get(name, 0.0) for name in variable_names],
            dtype=float,
        )
        if self._direction == "max":
            objective_vector = -objective_vector
        bounds = scipy.optimize.Bounds(
            [
                -np.inf if self.bounds[name][0] is None else self.bounds[name][0]
                for name in variable_names
            ],
            [
                np.inf if self.bounds[name][1] is None else self.bounds[name][1]
                for name in variable_names
            ],
        )
        integrality = np.array(
            [
                0
                if self._variable_types[name] == "continuous"
                else 1
                for name in variable_names
            ],
            dtype=int,
        )
        constraints = []
        constraint_names = []
        for constraint_name, row in self._constraints.items():
            coefficients = row["coefficients"]
            # Unknown variable names raise KeyError here on purpose: silently
            # dropping a coefficient would change the mathematical problem.
            row_names = list(coefficients)
            dense_row = scipy.sparse.coo_matrix(
                (
                    [coefficients[name] for name in row_names],
                    (
                        np.zeros(len(row_names), dtype=int),
                        [self._variables[name] for name in row_names],
                    ),
                ),
                shape=(1, len(variable_names)),
            )
            constraints.append(
                scipy.optimize.LinearConstraint(
                    dense_row,
                    lb=-np.inf if row["lb"] is None else row["lb"],
                    ub=np.inf if row["ub"] is None else row["ub"],
                )
            )
            constraint_names.append(constraint_name)
        options = {"disp": bool(verbosity)}
        forwarded = {"disp": bool(verbosity)}
        if presolve is not None:
            options["presolve"] = bool(presolve)
            forwarded["presolve"] = bool(presolve)
        if timeout is not None:
            options["time_limit"] = float(timeout)
            forwarded["time_limit"] = float(timeout)
        if mip_rel_gap is not None:
            options["mip_rel_gap"] = float(mip_rel_gap)
            forwarded["mip_rel_gap"] = float(mip_rel_gap)
        if node_limit is not None:
            options["node_limit"] = int(node_limit)
            forwarded["node_limit"] = int(node_limit)
        if continuous_lp_backend not in {"milp", "linprog"}:
            raise ValueError("Unknown continuous LP backend")
        lp_route = continuous_lp_backend == "linprog" and not integrality.any()
        if lp_route:
            eq_rows, eq_rhs, ub_rows, ub_rhs = [], [], [], []
            for constraint in constraints:
                lower, upper = float(np.asarray(constraint.lb).item()), float(np.asarray(constraint.ub).item())
                if np.isfinite(lower) and lower == upper:
                    eq_rows.append(constraint.A); eq_rhs.append(lower)
                else:
                    if np.isfinite(upper):
                        ub_rows.append(constraint.A); ub_rhs.append(upper)
                    if np.isfinite(lower):
                        ub_rows.append(-constraint.A); ub_rhs.append(-lower)
            # Gap and node limits are MILP-only and must not be advertised
            # as forwarded LP settings.
            options = {k: v for k, v in options.items() if k not in {"mip_rel_gap", "node_limit"}}
            if feasibility_tolerance is not None:
                tolerance = float(feasibility_tolerance)
                if not np.isfinite(tolerance) or tolerance < 1e-10:
                    raise ValueError("linprog HiGHS feasibility tolerance must be finite and >= 1e-10")
                options.update(primal_feasibility_tolerance=tolerance,
                               dual_feasibility_tolerance=tolerance)
            forwarded = dict(options)
            result = scipy.optimize.linprog(objective_vector,
                A_eq=scipy.sparse.vstack(eq_rows, format="csr") if eq_rows else None,
                b_eq=np.asarray(eq_rhs) if eq_rows else None,
                A_ub=scipy.sparse.vstack(ub_rows, format="csr") if ub_rows else None,
                b_ub=np.asarray(ub_rhs) if ub_rows else None,
                bounds=list(zip(bounds.lb, bounds.ub)), method="highs", options=options)
        else:
            result = scipy.optimize.milp(
                objective_vector, constraints=constraints, integrality=integrality,
                bounds=bounds, options=options,
            )
        raw_status = int(result.status)
        if raw_status == 1:
            status = (
                interface.TIME_LIMIT if timeout is not None
                else interface.ITERATION_LIMIT
            )
        else:
            status = SCIPY_MILP_STATUS.get(raw_status, interface.UNDEFINED)
        self._status = status
        # HiGHS proves optimality only at status 0; every other state keeps
        # its honest mapping above and never masquerades as optimal. When the
        # backend returns no point (x is None), cobra's Solution builder still
        # expects float primals (GLPK returns its last iterate); zeros are the
        # deterministic stand-in and the receipt marks them as not from the
        # backend. Only status == optimal carries HiGHS's actual solution.
        self._var_primals = (
            np.zeros(len(variable_names), dtype=float)
            if result.x is None
            else np.asarray(result.x, dtype=float)
        )
        self._f = None if result.fun is None else float(result.fun)
        self.last_result_receipt = {
            "backend": "scipy.optimize.linprog" if lp_route else "scipy.optimize.milp",
            "method": "highs" if lp_route else "HiGHS MILP",
            "scipy_version": scipy.__version__,
            "raw_status": raw_status,
            "status": status,
            "success": bool(result.success),
            "message": str(result.message),
            "objective_value": self._f,
            "primals_from_backend": result.x is not None,
            "mip_node_count": getattr(result, "mip_node_count", None),
            "mip_dual_bound": getattr(result, "mip_dual_bound", None),
            "mip_rel_gap": getattr(result, "mip_gap", None),
            "forwarded_options": forwarded,
            "integrality_tolerance_note": (
                "continuous LP; no integer-feasibility parameter applies" if lp_route else
                "HiGHS internal primal/integrality tolerances are not exposed "
                "through scipy.optimize.milp and remain at HiGHS defaults"
            ),
            "variable_count": len(variable_names),
            "constraint_count": len(constraint_names),
        }
        return status


def _symbolic_sum(terms, offset=0):
    """optlang.symbolics.add breaks on a single bare term; pad with Zero."""
    terms = list(terms)
    if not terms:
        return symbolics.Zero + offset
    return symbolics.add(symbolics.Zero, *terms) + offset


class Variable(interface.Variable, metaclass=inheritdocstring):
    def __init__(self, name, lb=None, ub=None, type="continuous", *args, **kwargs):
        if type not in ("continuous", "integer", "binary"):
            raise ValueError(f"Unsupported variable type: {type}")
        super(Variable, self).__init__(name, lb, ub, type, *args, **kwargs)

    @interface.Variable.lb.setter
    def lb(self, value):
        interface.Variable.lb.fset(self, value)
        if self.problem:
            self.problem.problem.set_variable_bounds(self.name, self.lb, self.ub)

    @interface.Variable.ub.setter
    def ub(self, value):
        interface.Variable.ub.fset(self, value)
        if self.problem:
            self.problem.problem.set_variable_bounds(self.name, self.lb, self.ub)

    def set_bounds(self, lb, ub):
        super(Variable, self).set_bounds(lb, ub)
        if self.problem:
            self.problem.problem.set_variable_bounds(self.name, lb, ub)

    @interface.Variable.type.setter
    def type(self, value):
        if value not in ("continuous", "integer", "binary"):
            raise ValueError(f"Unsupported variable type: {value}")
        super(Variable, Variable).type.fset(self, value)
        if self.problem:
            self.problem.problem.set_variable_type(self.name, value)

    def _get_primal(self):
        # GapFiller calls this low-level hook directly; the public `primal`
        # property in optlang's base class delegates here.
        if getattr(self, "problem", None) is not None:
            return self.problem.problem.get_var_primal(self.name)
        return None

    @property
    def dual(self):
        if getattr(self, "problem", None) is not None:
            raise NotImplementedError(
                "Duals are not available through scipy.optimize.milp."
            )
        return None

    @interface.Variable.name.setter
    def name(self, value):
        if getattr(self, "problem", None) is not None:
            self.problem.problem.change_var_name(self.name, value)
        super(Variable, Variable).name.fset(self, value)


class Constraint(interface.Constraint, metaclass=inheritdocstring):
    _INDICATOR_CONSTRAINT_SUPPORT = False

    def __init__(self, expression, sloppy=False, *args, **kwargs):
        super(Constraint, self).__init__(expression, sloppy=sloppy, *args, **kwargs)
        if not sloppy:
            if not self.is_Linear:
                raise ValueError(
                    f"HiGHS only supports linear constraints. {self} is not linear."
                )

    @property
    def primal(self):
        if getattr(self, "problem", None) is not None:
            return self.problem.problem.get_constraint_primal(self.name)
        return None

    @property
    def dual(self):
        if getattr(self, "problem", None) is not None:
            raise NotImplementedError(
                "Duals are not available through scipy.optimize.milp."
            )
        return None

    @interface.Constraint.lb.setter
    def lb(self, value):
        self._check_valid_lower_bound(value)
        if getattr(self, "problem", None) is not None:
            self.problem.problem.set_constraint_bounds(
                self.name, value, self.ub
            )
        self._lb = value

    @interface.Constraint.ub.setter
    def ub(self, value):
        self._check_valid_upper_bound(value)
        if getattr(self, "problem", None) is not None:
            self.problem.problem.set_constraint_bounds(
                self.name, self.lb, value
            )
        self._ub = value

    def coefficient_dict(self, names=True):
        if self.expression.is_Add:
            coefficient_dict = {
                variable: coefficient
                for variable, coefficient
                in self.expression.as_coefficients_dict().items()
                if variable.is_Symbol
            }
            coefficient_dict = {
                variable: float(coefficient)
                for variable, coefficient in coefficient_dict.items()
            }
        elif self.expression.is_Atom and self.expression.is_Symbol:
            coefficient_dict = {self.expression: 1}
        elif self.expression.is_Mul and len(self.expression.args) <= 2:
            args = self.expression.args
            coefficient_dict = {args[1]: float(args[0])}
        elif self.expression.is_Number:
            coefficient_dict = {}
        else:
            raise ValueError("Invalid expression: " + str(self.expression))
        if names is True:
            coefficient_dict = {
                variable.name: coefficient
                for variable, coefficient in coefficient_dict.items()
            }
        return coefficient_dict

    def set_linear_coefficients(self, coefficients):
        if self.problem is not None:
            self.problem.update()
            coefficients_dict = self.coefficient_dict(names=False)
            coefficients_dict.update(coefficients)
            self._expression = _symbolic_sum(
                variable * coefficient
                for variable, coefficient in coefficients_dict.items()
            )
            # cobra updates existing mass-balance rows through this method
            # (e.g. when a reaction is added); the stored backend row must
            # follow or the solved problem silently loses those terms.
            self.problem.problem.set_constraint_coefficients(
                self.name,
                {
                    variable.name: float(coefficient)
                    for variable, coefficient in coefficients_dict.items()
                },
            )
        else:
            raise Exception(
                "Can't change coefficients if constraint is not associated with a model."
            )

    def get_linear_coefficients(self, variables):
        if self.problem is not None:
            self.problem.update()
            coefficients_dict = self.coefficient_dict(names=False)
            return {v: coefficients_dict.get(v, 0) for v in variables}
        raise Exception(
            "Can't get coefficients if constraint is not associated with a model."
        )


class Objective(interface.Objective, metaclass=inheritdocstring):
    def __init__(self, expression, sloppy=False, **kwargs):
        super(Objective, self).__init__(expression, sloppy=sloppy, **kwargs)
        if not (sloppy or self.is_Linear):
            raise ValueError(
                f"HiGHS only supports linear objectives. {self} is not linear."
            )

    def _get_expression(self):
        if self.problem is not None:
            coefficients_dict = self.problem.problem.objective
            coefficients_dict = {
                self.problem._variables[name]: coefficient
                for name, coefficient in coefficients_dict.items()
                if name in self.problem._variables
            }
            self._expression = _symbolic_sum(
                (
                    variable * coefficient
                    for variable, coefficient in coefficients_dict.items()
                ),
                offset=self.problem.problem.offset,
            )
        return self._expression

    @property
    def value(self):
        if getattr(self, "problem", None) is not None:
            return self.problem.problem.objective_value
        return None

    @interface.Objective.direction.setter
    def direction(self, value):
        super(Objective, self.__class__).direction.fset(self, value)
        if getattr(self, "problem", None) is not None:
            self.problem.problem.direction = value

    def coefficient_dict(self):
        if self.expression.is_Add:
            coefficient_dict = {
                variable: coefficient
                for variable, coefficient
                in self.expression.as_coefficients_dict().items()
                if variable.is_Symbol
            }
        elif self.expression.is_Atom and self.expression.is_Symbol:
            coefficient_dict = {self.expression: 1}
        elif (
            self.expression.is_Mul
            and len(self.expression.args) <= 2
            and self.expression.args[1].is_Symbol
        ):
            args = self.expression.args
            coefficient_dict = {args[1]: float(args[0])}
        elif self.expression.is_Number:
            coefficient_dict = {}
        else:
            raise ValueError("Invalid expression: " + str(self.expression))
        return {
            variable.name: float(coefficient)
            for variable, coefficient in coefficient_dict.items()
        }

    def set_linear_coefficients(self, coefficients):
        if self.problem is not None:
            self.problem.update()
            self.problem.problem.objective.update(
                {variable.name: coefficient for variable, coefficient in coefficients.items()}
            )
        else:
            raise Exception(
                "Can't change coefficients if objective is not associated with a model."
            )

    def get_linear_coefficients(self, variables):
        if self.problem is not None:
            self.problem.update()
            return {
                variable: float(
                    self.problem.problem.objective.get(variable.name, 0)
                )
                for variable in variables
            }
        raise Exception(
            "Can't get coefficients if objective is not associated with a model."
        )


class Configuration(interface.MathematicalProgrammingConfiguration, metaclass=inheritdocstring):
    def __init__(
        self,
        problem=None,
        verbosity=0,
        timeout=None,
        presolve=True,
        tolerances=None,
        mip_node_limit=None,
        continuous_lp_backend="milp",
        *args,
        **kwargs,
    ):
        self._verbosity = verbosity
        self._timeout = timeout
        self._presolve = True
        self._presolve_auto_mapped = False
        self._mip_node_limit = mip_node_limit
        self.continuous_lp_backend = continuous_lp_backend
        self._stored_tolerances = {
            "feasibility": 1e-7,
            "integrality": 1e-6,
            "optimality": 0.0,
        }
        super(Configuration, self).__init__(problem=problem, *args, **kwargs)
        # Route through the validating setter: Configuration.clone() injects
        # other backends' defaults (e.g. GLPK's presolve="auto"), which must
        # not bypass validation when a model switches solver mid-life.
        self.presolve = presolve
        if tolerances:
            for name, value in tolerances.items():
                setattr(self.tolerances, name, value)

    def _tolerance_functions(self):
        stored = self._stored_tolerances
        return {
            name: (
                (lambda name=name: stored[name]),
                (lambda value, name=name: stored.__setitem__(name, float(value))),
            )
            for name in ("feasibility", "integrality", "optimality")
        }

    @classmethod
    def clone(cls, config, problem=None, **kwargs):
        # optlang's base clone copies only its standard properties. Preserve
        # this interface's explicit LP selection on model copies as well.
        cloned = super(Configuration, cls).clone(config, problem=problem, **kwargs)
        cloned.continuous_lp_backend = getattr(config, "continuous_lp_backend", "milp")
        cloned.mip_node_limit = getattr(config, "mip_node_limit", None)
        return cloned

    @property
    def continuous_lp_backend(self):
        return self._continuous_lp_backend

    @continuous_lp_backend.setter
    def continuous_lp_backend(self, value):
        if value not in {"milp", "linprog"}:
            raise ValueError("continuous_lp_backend must be milp or linprog")
        self._continuous_lp_backend = value

    @property
    def verbosity(self):
        return self._verbosity

    @verbosity.setter
    def verbosity(self, value):
        self._verbosity = value

    @property
    def timeout(self):
        return self._timeout

    @timeout.setter
    def timeout(self, value):
        self._timeout = value

    @property
    def presolve(self):
        return self._presolve

    @presolve.setter
    def presolve(self, value):
        # optlang's Model.optimize() auto-retry toggles presolve through this
        # setter with "auto"; scipy.optimize.milp exposes no such retry mode,
        # so "auto" maps to enabled and the mapping is recorded for audit.
        if value == "auto":
            value = True
            self._presolve_auto_mapped = True
        if value not in (True, False):
            raise ValueError(
                "The HiGHS interface accepts only True/False presolve "
                "('auto' is mapped to True and recorded)."
            )
        self._presolve = value

    @property
    def presolve_auto_mapped(self):
        """True when an 'auto' presolve value was mapped to True."""
        return self._presolve_auto_mapped

    @property
    def mip_node_limit(self):
        return self._mip_node_limit

    @mip_node_limit.setter
    def mip_node_limit(self, value):
        self._mip_node_limit = value


class Model(interface.Model, metaclass=inheritdocstring):
    def __getstate__(self):
        state = super(Model, self).__getstate__()
        state["_highs_configuration"] = {
            "continuous_lp_backend": self.configuration.continuous_lp_backend,
            "tolerances": self.configuration.tolerances.to_dict(),
            "mip_node_limit": self.configuration.mip_node_limit,
            "verbosity": self.configuration.verbosity,
            "timeout": self.configuration.timeout,
            "presolve": self.configuration.presolve,
        }
        return state

    def __setstate__(self, state):
        state = dict(state)
        config = state.pop("_highs_configuration", None)
        super(Model, self).__setstate__(state)
        if config is not None:
            self.configuration = Configuration(problem=self, **config)

    def _initialize_problem(self):
        self.problem = Problem()

    def _initialize_model_from_problem(self, problem):
        if isinstance(problem, Problem):
            self.problem = problem
        else:
            raise TypeError(
                "Problem must be an instance of optlang_highs_interface.Problem, "
                f"not {type(problem)!r}"
            )

    def _add_variables(self, variables):
        super(Model, self)._add_variables(variables)
        for variable in variables:
            self.problem.add_variable(variable.name, variable.type)
            self.problem.set_variable_bounds(
                variable.name, variable.lb, variable.ub
            )

    def _remove_variables(self, variables):
        variable_names = [variable.name for variable in variables]
        super(Model, self)._remove_variables(variables)
        for name in variable_names:
            self.problem.remove_variable(name)

    def _add_constraints(self, constraints, sloppy=False):
        if sloppy is False:
            for constraint in constraints:
                if not constraint.is_Linear:
                    raise ValueError(
                        "HiGHS only supports linear constraints: "
                        f"{constraint} is not linear."
                    )
        super(Model, self)._add_constraints(constraints, sloppy=sloppy)
        for constraint in constraints:
            self.problem.add_constraint(
                constraint.name,
                constraint.coefficient_dict(),
                lb=constraint.lb,
                ub=constraint.ub,
            )

    def _remove_constraints(self, constraints):
        constraint_names = [constraint.name for constraint in constraints]
        super(Model, self)._remove_constraints(constraints)
        for name in constraint_names:
            self.problem.remove_constraint(name)

    def _optimize(self):
        return self.problem.optimize(
            verbosity=self.configuration.verbosity,
            timeout=self.configuration.timeout,
            presolve=self.configuration.presolve,
            mip_rel_gap=self.configuration.tolerances.optimality,
            node_limit=self.configuration.mip_node_limit,
            continuous_lp_backend=self.configuration.continuous_lp_backend,
            feasibility_tolerance=self.configuration.tolerances.feasibility,
        )

    @interface.Model.objective.setter
    def objective(self, value):
        super(Model, Model).objective.fset(self, value)
        value.problem = None
        if value is None:
            self.problem.objective = {}
            self.problem.offset = 0.0
        else:
            offset, coefficients, _ = parse_optimization_expression(value)
            self.problem.objective = OrderedDict(
                (variable.name, coefficient)
                for variable, coefficient in coefficients.items()
            )
            self.problem.offset = offset
            self.problem.direction = value.direction
        value.problem = self


def register_highs_interface():
    """Register this module as cobra/optlang solver ``highs``; idempotent.

    cobra builds its solver lookup from the shared ``cobra.util.solver.solvers``
    dict (the same object backs ``cobra.core.configuration.SOLVERS``), so one
    insertion covers ``cobra.Configuration().solver`` and ``model.solver``.
    """
    import sys

    import cobra.util.solver

    cobra.util.solver.solvers["highs"] = sys.modules[__name__]
    return "highs"
