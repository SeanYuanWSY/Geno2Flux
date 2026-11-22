"""Formal reaction-pool audit fixtures use the author's emitted row format."""

import hashlib
import json
import shutil
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from unittest.mock import patch

import cobra
import pandas as pd
import pytest


BUNDLE_ROOT = Path(__file__).resolve().parents[1] / "third_party" / "alphagem"
sys.path.insert(0, str(BUNDLE_ROOT))

_tqdm_stub = ModuleType("tqdm")
_tqdm_stub.tqdm = lambda iterable: iterable
with patch.dict(sys.modules, {"tqdm": _tqdm_stub}):
    from src import add_reactions_based_pool as pool_builder  # noqa: E402
from src import pipeline_runner as runner_module  # noqa: E402
from src.pipeline_config import get_profile, resolve_profile_for_reference  # noqa: E402
from src.pipeline_runner import PipelineRun, PipelineRunner  # noqa: E402


def _write_pool(context, ids, monkeypatch):
    pool_path = context["pool_path"]
    pool_path.write_text(
        json.dumps({"reactions": [{"id": value} for value in ids]}),
        encoding="utf-8",
    )
    pool_record = context["runner"]._file_record(pool_path)
    context["runner"].preflight_by_stage.setdefault("reaction_pool", {})[
        "reaction_pool"
    ] = {
        "observed_path": str(pool_path),
        **pool_record,
    }
    monkeypatch.setattr(
        runner_module, "FROZEN_REACTION_POOL_SHA256", pool_record["sha256"]
    )


def _save_model(context, *, balanced=True, include_missing=False):
    model = cobra.Model("fixture")
    water = cobra.Metabolite("water_c", formula="H2O", charge=0, compartment="c")
    product = cobra.Metabolite(
        "product_c", formula="H2O" if balanced else "H2", charge=0,
        compartment="c",
    )
    added = cobra.Reaction("rhea_1_c")
    added.add_metabolites({water: -1, product: 1})
    model.add_reactions([added])
    cobra.io.save_yaml_model(model, str(context["draft_path"]))
    if include_missing:
        missing = cobra.Reaction("rhea_80727_c")
        missing.add_metabolites({water: -1, product: 1})
        model.add_reactions([missing])
    cobra.io.save_yaml_model(model, str(context["model_path"]))


def _write_audit(context, statuses):
    ids = ("rhea_1_c", "rhea_80727_c", "rhea_2_c")
    candidates = [
        {
            "input_row": index,
            "reaction_id": reaction_id,
            "status": status,
            "unmatched_metabolites": [],
        }
        for index, (reaction_id, status) in enumerate(zip(ids, statuses))
    ]
    counts = {}
    for status in statuses:
        counts[status] = counts.get(status, 0) + 1
    audit = {
        "input_candidate_count": len(candidates),
        "status_counts": counts,
        "candidates": candidates,
    }
    context["audit_path"].write_text(json.dumps(audit), encoding="utf-8")
    return audit


def _write_coverage(context, audit):
    coverage_path = context["working"] / "reaction_pool_coverage.json"
    coverage = pool_builder._build_reaction_pool_coverage(
        audit["candidates"],
        SimpleNamespace(reactions=[]),
        "author_skip",
        "reject",
    )
    coverage_path.write_text(
        json.dumps(coverage, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return coverage_path


@pytest.fixture
def context(tmp_path, monkeypatch):
    working = tmp_path / "working" / "sample_fixture"
    working.mkdir(parents=True)
    models = tmp_path / "models"
    models.mkdir()
    pool_path = models / "reaction_pool_final.json"
    score_path = working / "sample_fixturegpr_score.xlsx"
    model_path = working / "sample_fixture-GEM_withgaps.yml"
    draft_path = working / "model2sample_fixture.yml"
    audit_path = working / "reaction_pool_audit.json"
    pd.DataFrame({"rhea": [1, 80727, 2]}).to_excel(score_path, index=False)
    runner = PipelineRunner(tmp_path)
    run = PipelineRun(
        name="sample_fixture",
        refname="yeast",
        fasta=tmp_path / "fixture.faa",
        mode="plmsearch",
        profile=resolve_profile_for_reference(get_profile("paper"), "yeast"),
    )
    result = {
        "runner": runner,
        "run": run,
        "working": working,
        "pool_path": pool_path,
        "score_path": score_path,
        "model_path": model_path,
        "draft_path": draft_path,
        "audit_path": audit_path,
    }
    reference_path = models / "yeast-GEM.xml"
    reference_path.write_text("fixture reference\n", encoding="utf-8")
    runner.preflight_by_stage["reaction_pool"] = {
        "reference_model": {
            "observed_path": str(reference_path),
            **runner._file_record(reference_path),
        }
    }
    _write_pool(result, ("rhea_1_c", "rhea_2_c"), monkeypatch)
    _save_model(result)
    audit = _write_audit(result, ("added", "missing_from_reaction_pool", "rejected_mass_or_charge_balance"))
    _write_coverage(result, audit)
    return result


def test_author_skipped_absence_has_versioned_bound_receipt(context):
    runner, run = context["runner"], context["run"]
    artifacts = runner._verify_stage_outputs("reaction_pool", run)
    receipt_path = context["working"] / "reaction_pool_validation_v3.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    assert receipt["policy_version"] == 4
    assert receipt["reaction_pool_policies"] == {
        "already_present_policy": "author_skip",
        "conflict_policy": "reject",
    }
    assert receipt["missing_from_reaction_pool"] == [
        {"input_row": 1, "reaction_id": "rhea_80727_c"}
    ]
    assert receipt["status_counts"] == {
        "added": 1,
        "missing_from_reaction_pool": 1,
        "rejected_mass_or_charge_balance": 1,
    }
    assert receipt["candidate_audit"] == runner._file_record(context["audit_path"])
    assert artifacts[str(receipt_path.relative_to(runner.repo_root))] == runner._file_record(receipt_path)


def test_false_pool_absence_remains_fatal(context, monkeypatch):
    _write_pool(context, ("rhea_1_c", "rhea_2_c", "rhea_80727_c"), monkeypatch)
    with pytest.raises(RuntimeError, match="falsely marks a present"):
        context["runner"]._verify_stage_outputs("reaction_pool", context["run"])
    assert not (context["working"] / "reaction_pool_validation_v3.json").exists()


def test_changed_pool_without_matching_preflight_is_fatal(context):
    context["pool_path"].write_text('{"reactions": []}', encoding="utf-8")
    with pytest.raises(RuntimeError, match="frozen snapshot or preflight"):
        context["runner"]._verify_stage_outputs("reaction_pool", context["run"])


def test_candidate_conservation_checks_actual_spreadsheet(context):
    audit = json.loads(context["audit_path"].read_text(encoding="utf-8"))
    audit["candidates"][1]["reaction_id"] = "rhea_other_c"
    context["audit_path"].write_text(json.dumps(audit), encoding="utf-8")
    with pytest.raises(RuntimeError, match="candidate differs from input row"):
        context["runner"]._verify_stage_outputs("reaction_pool", context["run"])


def test_explicit_released_code_exclusion_stays_fatal(context):
    _write_audit(context, ("added", "explicitly_excluded_released_code", "rejected_mass_or_charge_balance"))
    with pytest.raises(RuntimeError, match="not evaluated by the published QC contract"):
        context["runner"]._verify_stage_outputs("reaction_pool", context["run"])
    assert not (context["working"] / "reaction_pool_validation_v3.json").exists()


def test_added_imbalanced_reaction_stays_fatal(context):
    _save_model(context, balanced=False)
    with pytest.raises(RuntimeError, match="added an imbalanced reaction"):
        context["runner"]._verify_stage_outputs("reaction_pool", context["run"])


def test_claimed_pool_absence_cannot_hide_model_addition(context):
    _save_model(context, include_missing=True)
    with pytest.raises(RuntimeError, match="was added to the model"):
        context["runner"]._verify_stage_outputs("reaction_pool", context["run"])


def test_real_five_reactions_are_rejected_before_model_insertion(tmp_path):
    pool = cobra.io.load_json_model(
        str(BUNDLE_ROOT / "models" / "reaction_pool_final.json")
    )
    model = cobra.Model("sample_draft_fixture")
    model.add_metabolites([
        cobra.Metabolite(
            "chebi57880_c", formula="C9H15O9P", charge=-1, compartment="c"
        ),
        cobra.Metabolite(
            "chebi18072_c", formula="C5H9NO3", charge=0, compartment="c"
        ),
    ])
    ids = (
        "rhea_11580_c", "rhea_14789_c", "rhea_19877_c",
        "rhea_21368_c", "rhea_23736_c",
    )
    records = []
    for reaction_id in ids:
        record = {"reaction_id": reaction_id, "candidate_source": "frozen_reaction_pool"}
        pool_builder.admit_candidate_reaction(
            pool.reactions.get_by_id(reaction_id).copy(), model, record
        )
        records.append(record)
    assert [record["reaction_id"] for record in records] == list(ids)
    assert all(
        record["status"] == "rejected_postmerge_metabolite_identity_conflict"
        and record["rejection_reason"] == "postmerge_metabolite_identity_conflict"
        and len(record["identity_conflicts"]) == 1
        for record in records
    )
    assert [
        record["identity_conflicts"][0]["metabolite_id"]
        for record in records
    ] == ["chebi57880_c"] * 3 + ["chebi18072_c"] * 2
    saved = tmp_path / "no_conflicting_reactions.yml"
    cobra.io.save_yaml_model(model, str(saved))
    reloaded = cobra.io.load_yaml_model(str(saved))
    assert not {reaction.id for reaction in reloaded.reactions} & set(ids)


def test_source_imbalance_keeps_original_rejection_before_identity_check():
    pool = cobra.io.load_json_model(
        str(BUNDLE_ROOT / "models" / "reaction_pool_final.json")
    )
    model = cobra.Model("sample_draft_fixture")
    for metabolite in pool.reactions.get_by_id("rhea_10464_c").metabolites:
        model.add_metabolites([cobra.Metabolite(
            metabolite.id, formula="C1", charge=0,
            compartment=metabolite.compartment,
        )])
    record = {"reaction_id": "rhea_10464_c"}
    pool_builder.admit_candidate_reaction(
        pool.reactions.get_by_id("rhea_10464_c").copy(), model, record
    )
    assert record["status"] == "rejected_mass_or_charge_balance"
    assert "identity_conflicts" not in record
    assert "rhea_10464_c" not in model.reactions


def test_real_five_conflicts_are_bound_in_formal_receipt(tmp_path):
    ids = (11580, 14789, 19877, 21368, 23736)
    name = "sample_real_conflicts"
    working = tmp_path / "working" / name
    working.mkdir(parents=True)
    models = tmp_path / "models"
    models.mkdir()
    pool_path = models / "reaction_pool_final.json"
    shutil.copyfile(BUNDLE_ROOT / "models" / "reaction_pool_final.json", pool_path)
    reference_path = models / "yeast-GEM.xml"
    reference_path.write_text("fixture reference\n", encoding="utf-8")
    score_path = working / f"{name}gpr_score.xlsx"
    pd.DataFrame({"rhea": ids}).to_excel(score_path, index=False)
    draft = cobra.Model(name)
    draft.add_metabolites([
        cobra.Metabolite(
            "chebi57880_c", formula="C9H15O9P", charge=-1, compartment="c"
        ),
        cobra.Metabolite(
            "chebi18072_c", formula="C5H9NO3", charge=0, compartment="c"
        ),
    ])
    draft_path = working / f"model2{name}.yml"
    saved_path = working / f"{name}-GEM_withgaps.yml"
    cobra.io.save_yaml_model(draft, str(draft_path))
    pool = cobra.io.load_json_model(str(pool_path))
    candidates = []
    for index, rhea_id in enumerate(ids):
        reaction_id = f"rhea_{rhea_id}_c"
        record = {
            "input_row": index,
            "reaction_id": reaction_id,
            "candidate_source": "frozen_reaction_pool",
            "unmatched_metabolites": [],
        }
        pool_builder.admit_candidate_reaction(
            pool.reactions.get_by_id(reaction_id).copy(), draft, record
        )
        candidates.append(record)
    cobra.io.save_yaml_model(draft, str(saved_path))
    audit_path = working / "reaction_pool_audit.json"
    audit_path.write_text(json.dumps({
        "input_candidate_count": 5,
        "status_counts": {"rejected_postmerge_metabolite_identity_conflict": 5},
        "candidates": candidates,
    }), encoding="utf-8")
    coverage = pool_builder._build_reaction_pool_coverage(
        candidates,
        SimpleNamespace(reactions=list(pool.reactions)),
        "author_skip",
        "reject",
    )
    (working / "reaction_pool_coverage.json").write_text(
        json.dumps(coverage, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    runner = PipelineRunner(tmp_path)
    runner.preflight_by_stage["reaction_pool"] = {
        "reaction_pool": {
            "observed_path": str(pool_path), **runner._file_record(pool_path),
        },
        "reference_model": {
            "observed_path": str(reference_path),
            **runner._file_record(reference_path),
        },
    }
    run = PipelineRun(
        name=name,
        refname="yeast",
        fasta=tmp_path / "fixture.faa",
        mode="plmsearch",
        profile=resolve_profile_for_reference(get_profile("paper"), "yeast"),
    )
    artifacts = runner._verify_stage_outputs("reaction_pool", run)
    receipt_path = working / "reaction_pool_validation_v3.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    assert receipt["postmerge_identity_conflict_count"] == 5
    assert [row["reaction_id"] for row in receipt["postmerge_identity_conflicts"]] == [
        f"rhea_{value}_c" for value in ids
    ]
    assert receipt["reaction_pool"] == runner._file_record(pool_path)
    assert receipt["candidate_spreadsheet"] == runner._file_record(score_path)
    assert receipt["draft_model"] == runner._file_record(draft_path)
    assert receipt["reference_model"] == runner._file_record(reference_path)
    assert receipt["saved_model"] == runner._file_record(saved_path)
    assert artifacts[str(receipt_path.relative_to(tmp_path))] == runner._file_record(receipt_path)
    saved_model = cobra.io.load_yaml_model(str(saved_path))
    assert not saved_model.reactions


def test_strict_admission_persists_only_balanced_reaction(tmp_path):
    model = cobra.Model("balanced_fixture")
    a = cobra.Metabolite("a_c", formula="H2O", charge=0, compartment="c")
    b = cobra.Metabolite("b_c", formula="H2O", charge=0, compartment="c")
    reaction = cobra.Reaction("rhea_balanced_c")
    reaction.add_metabolites({a: -1, b: 1})
    record = {"reaction_id": reaction.id}
    pool_builder.admit_candidate_reaction(reaction, model, record)
    assert record["status"] == "added"
    path = tmp_path / "balanced.yml"
    cobra.io.save_yaml_model(model, str(path))
    saved = cobra.io.load_yaml_model(str(path))
    assert saved.reactions.get_by_id("rhea_balanced_c").check_mass_balance() == {}


def test_packaged_pool_hash_is_the_policy_anchor():
    pool_path = BUNDLE_ROOT / "models" / "reaction_pool_final.json"
    assert hashlib.sha256(pool_path.read_bytes()).hexdigest() == (
        runner_module.FROZEN_REACTION_POOL_SHA256
    )
