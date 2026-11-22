import json
from pathlib import Path
import subprocess
import sys

import pytest

import scripts.prepare_alphagem_runtime as bundle


def test_bundled_alphagem_source_is_complete_and_regular():
    manifest, records = bundle.source_records()
    assert manifest["source_commit"] == "788a00b21ff1e516065166e2577dd8c499ef872d"
    assert len(records) == 92
    assert {item["path"] for item in records} >= {
        "AlphaGEM.py",
        "LICENSE",
        "src/optlang_highs_interface.py",
        "src/pipeline_config.py",
        "src/pipeline_runner.py",
        "models/yeast-GEM.xml",
        "models/reaction_pool_final.json",
    }
    help_result = subprocess.run(
        [sys.executable, str(bundle.BUNDLE_ROOT / "AlphaGEM.py"), "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--profile" in help_result.stdout
    assert "--stop-after" in help_result.stdout


def test_runtime_install_uses_real_files_and_never_overwrites(tmp_path):
    _, records = bundle.source_records()
    runtime = tmp_path / "runtime"
    installed = bundle.install_runtime(runtime, records)
    assert len(installed["installed"]) == len(records)
    assert not any(path.is_symlink() for path in runtime.rglob("*"))
    assert not bundle.inspect_runtime(runtime, records)["missing"]
    assert installed["runtime_source_commit"] == bundle.runtime_source_commit(runtime, records)
    assert subprocess.run(
        ["git", "status", "--porcelain=v1"], cwd=runtime,
        capture_output=True, text=True, check=True,
    ).stdout == ""
    (runtime / "data_available").mkdir()
    (runtime / "data_available" / "yeast.fasta").write_text(">fixture\nA\n")
    (runtime / "models" / "esm.pt").write_text("fixture\n")
    identity_result = subprocess.run(
        [
            sys.executable, "-c",
            "import json; from src.pipeline_runner import PipelineRunner; "
            "print(json.dumps(PipelineRunner('.')._source_code_identity()))",
        ],
        cwd=runtime, capture_output=True, text=True, check=True,
    )
    identity = json.loads(identity_result.stdout)
    assert identity["git_repository"] is True
    assert identity["commit"] == installed["runtime_source_commit"]
    assert identity["dirty"] is False
    again = bundle.install_runtime(runtime, records)
    assert again["installed"] == []

    protected = runtime / "src" / "pipeline_config.py"
    protected.write_text("locally changed\n", encoding="utf-8")
    with pytest.raises(ValueError, match="conflicts"):
        bundle.install_runtime(runtime, records)
    assert protected.read_text(encoding="utf-8") == "locally changed\n"


def test_runtime_install_rejects_symlinked_source_path(tmp_path):
    _, records = bundle.source_records()
    runtime = tmp_path / "runtime"
    runtime.mkdir()
    (runtime / "src").symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        bundle.install_runtime(runtime, records)
    assert not (tmp_path / "pipeline_config.py").exists()


def test_runtime_install_rejects_symlinked_parent_path(tmp_path):
    _, records = bundle.source_records()
    actual_parent = tmp_path / "actual"
    actual_parent.mkdir()
    linked_parent = tmp_path / "linked"
    linked_parent.symlink_to(actual_parent, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        bundle.install_runtime(linked_parent / "runtime", records)
    assert not (actual_parent / "runtime").exists()


def test_runtime_install_rejects_source_repository_destination():
    with pytest.raises(ValueError, match="outside the source repository"):
        bundle.require_real_runtime_path(bundle.REPOSITORY_ROOT / "runtime")


def test_cli_canonicalizes_symlinked_parent_before_preparing(tmp_path, capsys):
    actual_parent = tmp_path / "actual"
    actual_parent.mkdir()
    linked_parent = tmp_path / "linked"
    linked_parent.symlink_to(actual_parent, target_is_directory=True)

    assert bundle.main([
        "--prepare", "--runtime-root", str(linked_parent / "runtime")
    ]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["runtime_root"] == str(actual_parent / "runtime")
    assert (actual_parent / "runtime" / "AlphaGEM.py").is_file()
    assert not any(path.is_symlink() for path in (actual_parent / "runtime").rglob("*"))


def test_documented_yeast9_target_is_real_and_optimizable():
    import cobra

    model = cobra.io.read_sbml_model(
        str(bundle.BUNDLE_ROOT / "models" / "yeast-GEM.xml")
    )
    config = json.loads(
        (bundle.REPOSITORY_ROOT / "configs" / "unified_pipeline.example.json").read_text()
    )
    target = model.reactions.get_by_id(config["target_reaction"])
    assert target.id == "r_2056"
    assert config["biomass"] in model.reactions
    with model:
        model.objective = target
        solution = model.optimize()
    assert solution.status == "optimal"
    assert solution.objective_value > 0
