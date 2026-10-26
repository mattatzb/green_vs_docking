"""Packaging, relocation, command dispatch, and public imports."""
import importlib
import os
import pkgutil
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from green_vs_docking import paths
from green_vs_docking.cli import main as cli

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("command", ["run", "score", "plot", "table"])
def test_shared_entrypoint_dispatches_without_changing_arguments(command, monkeypatch):
    captured = []
    def load(name, package):
        assert name == "." + command
        return SimpleNamespace(main=lambda argv: captured.append(argv))
    monkeypatch.setattr(cli.importlib, "import_module", load)
    cli.main([command, "--help"])
    assert captured == [["--help"]]


def test_default_command_is_complete_pipeline(monkeypatch):
    captured = []
    def load(name, package):
        assert name == ".run"
        return SimpleNamespace(main=lambda argv: captured.append(argv))
    monkeypatch.setattr(cli.importlib, "import_module", load)
    cli.main([])
    cli.main(["--output-dir", "custom"])
    assert captured == [[], ["--output-dir", "custom"]]


def test_all_source_modules_import():
    import green_vs_docking
    for module in pkgutil.walk_packages(green_vs_docking.__path__, green_vs_docking.__name__ + "."):
        importlib.import_module(module.name)


@pytest.mark.parametrize("script,args", [
    ("main.py", []), ("main.py", ["run"]), ("main.py", ["score"]),
    ("main.py", ["plot"]), ("main.py", ["table"]),
])
def test_source_launchers_work_from_another_directory(tmp_path, script, args):
    result = subprocess.run([sys.executable, str(ROOT / script), *args, "--help"],
                            cwd=tmp_path, check=True, capture_output=True, text=True)
    assert "usage:" in result.stdout


def test_paths_survive_checkout_rename_and_spaces(tmp_path):
    checkout = tmp_path / "merged repository"
    shutil.copytree(ROOT / "src/green_vs_docking", checkout / "src/green_vs_docking")
    (checkout / "data").mkdir()
    shutil.copyfile(ROOT / "data/blocks.smi", checkout / "data/blocks.smi")
    shutil.copyfile(ROOT / "environment.yml", checkout / "environment.yml")
    env = dict(os.environ, PYTHONPATH=str(checkout / "src") + os.pathsep + os.environ.get("PYTHONPATH", ""))
    env.pop("GREEN_VS_DOCKING_ROOT", None)
    result = subprocess.run([sys.executable, "-c", "from green_vs_docking.paths import REPO_ROOT; print(REPO_ROOT)"],
                            cwd=tmp_path, env=env, check=True, capture_output=True, text=True)
    assert result.stdout.strip() == str(checkout)


def test_new_pipeline_defaults_keep_reference_outputs_separate():
    assert paths.DEFAULT_GENERATED_ROOT == ROOT / "outputs/9KQ3"
    assert paths.REFERENCE_OUTPUTS_ROOT == ROOT / "9KQ3_outputs"
    assert paths.DEFAULT_GENERATED_ROOT != paths.REFERENCE_OUTPUTS_ROOT


def test_explicit_repository_root_override(tmp_path, monkeypatch):
    monkeypatch.setenv("GREEN_VS_DOCKING_ROOT", str(tmp_path))
    assert paths.repository_root() == tmp_path
