"""Check CLI input and output paths without launching analysis jobs."""

import sys
from pathlib import Path

import pytest

from green_vs_docking.cli import score as scoring
from green_vs_docking.cli import plot as plots
from green_vs_docking.cli import table as table_summary
from green_vs_docking.paths import DEFAULT_GENERATED_ROOT


ROOT = Path(__file__).resolve().parents[1]


def test_table_cli_defaults_use_checkout_paths(tmp_path, monkeypatch):
    captured = {}
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["table"])
    monkeypatch.setattr(table_summary, "write_table_summary", lambda **args: captured.update(args))
    table_summary.main()
    assert captured["outputs_root"] == ROOT / "9KQ3_outputs"
    assert captured["oracle_root"] == ROOT / "data"
    assert captured["output_path"] == DEFAULT_GENERATED_ROOT / "latex_table_summary.csv"


def test_table_cli_custom_root_and_explicit_destination(tmp_path, monkeypatch):
    captured = {}
    monkeypatch.setattr(table_summary, "write_table_summary", lambda **args: captured.update(args))
    arguments = ["table", "--outputs-root", str(tmp_path)]
    monkeypatch.setattr(sys, "argv", arguments)
    table_summary.main()
    assert captured["output_path"] == DEFAULT_GENERATED_ROOT / "latex_table_summary.csv"
    custom = tmp_path / "explicit.csv"
    monkeypatch.setattr(sys, "argv", arguments + ["--output", str(custom)])
    table_summary.main()
    assert captured["output_path"] == custom


@pytest.mark.parametrize("blocks", [False, True])
def test_plot_cli_default_destinations_are_separate_from_inputs(tmp_path, monkeypatch, blocks):
    captured = {}
    arguments = ["plot", "--root", str(tmp_path)]
    if blocks:
        arguments.append("--prefer-blocks-file")
    else:
        arguments.append("--no-prefer-blocks-file")
    monkeypatch.setattr(sys, "argv", arguments)
    monkeypatch.setattr(plots, "load_experiment_points", lambda *args, **kwargs: None)
    monkeypatch.setattr(plots, "plot_points", lambda df, path: captured.update(plot=path))
    monkeypatch.setattr(plots, "build_resume_table", lambda root, path: None)
    monkeypatch.setattr(plots, "build_resume_table_from_comparisons",
                        lambda **args: captured.update(resume=args["resume_path"]))
    plots.main()
    suffix = "_blocks" if blocks else ""
    assert captured["plot"] == DEFAULT_GENERATED_ROOT / f"green_vs_docking_summary{suffix}.pdf"
    assert captured["resume"] == DEFAULT_GENERATED_ROOT / "resume.csv"


def test_plot_cli_default_root_is_checkout(tmp_path, monkeypatch):
    captured = {}
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["plot"])
    monkeypatch.setattr(plots, "load_experiment_points", lambda root, **kw: captured.update(root=root))
    monkeypatch.setattr(plots, "plot_points", lambda frame, path: captured.update(plot=path))
    monkeypatch.setattr(plots, "build_resume_table", lambda *args: None)
    monkeypatch.setattr(plots, "build_resume_table_from_comparisons",
                        lambda **args: captured.update(resume=args["resume_path"]))
    plots.main()
    assert captured["root"] == ROOT / "9KQ3_outputs"
    assert captured["plot"] == DEFAULT_GENERATED_ROOT / "green_vs_docking_summary_blocks.pdf"
    assert captured["resume"] == DEFAULT_GENERATED_ROOT / "resume.csv"


def test_explicit_cli_paths_are_preserved(tmp_path, monkeypatch):
    captured = {}
    custom = tmp_path / "explicit.pdf"
    resume = tmp_path / "explicit.csv"
    monkeypatch.setattr(sys, "argv", ["plot", "--root", str(tmp_path),
                                      "--output", str(custom), "--resume-output", str(resume),
                                      "--prefer-blocks-file"])
    monkeypatch.setattr(plots, "load_experiment_points", lambda *args, **kwargs: None)
    monkeypatch.setattr(plots, "plot_points", lambda df, path: captured.update(plot=path))
    monkeypatch.setattr(plots, "build_resume_table_from_comparisons",
                        lambda **args: captured.update(resume=args["resume_path"]))
    plots.main()
    assert captured["plot"] == custom
    assert captured["resume"] == resume


def test_joint_kde_default_destination_is_separate_from_inputs(tmp_path, monkeypatch):
    captured = {}
    monkeypatch.setattr(plots, "load_experiment_rows", lambda *args, **kwargs: None)
    monkeypatch.setattr(plots, "sample_rows_by_label", lambda *args: None)
    monkeypatch.setattr(plots, "plot_joint_kde", lambda frame, path, **kw: captured.update(path=path))
    plots.main(["--root", str(tmp_path), "--joint-kde-only"])
    assert captured["path"] == DEFAULT_GENERATED_ROOT / "green_vs_docking_joint_kde_blocks.png"


def test_namerxn_is_found_on_path(tmp_path, monkeypatch):
    binary = tmp_path / "namerxn"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(0o755)
    monkeypatch.delenv("NAMERXN_BINARY", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path))
    parser = scoring.build_parser()
    args = parser.parse_args(["--experiment", f"baseline={ROOT / 'data'}"])
    assert args.namerxn_binary == binary


def test_namerxn_environment_and_explicit_override(tmp_path, monkeypatch):
    binary = tmp_path / "custom-namerxn"
    monkeypatch.setenv("NAMERXN_BINARY", str(binary))
    parser = scoring.build_parser()
    arguments = ["--experiment", f"baseline={ROOT / 'data'}"]
    assert parser.parse_args(arguments).namerxn_binary == binary
    assert parser.parse_args(arguments + ["--namerxn-binary", "override"]).namerxn_binary == Path("override")


def test_analysis_environment_matches_shared_versions():
    entries = {line.strip().removeprefix("- ") for line in (ROOT / "environment.yml").read_text().splitlines()
               if line.lstrip().startswith("- ")}
    assert {"python=3.9.7", "pip=25.0", "numpy=1.26.4", "pandas=2.1.4", "scipy=1.13.1",
            "rdkit==2025.3.5", "matplotlib==3.9.2", "pillow==11.1.0", "seaborn==0.13.2",
            "pytest==8.4.1"} <= entries
