"""The default CLI reproduces the paper; optional modes keep plot/CSV aligned."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from green_vs_docking.cli import plot as plots
from green_vs_docking.analysis import summary as statistics
from green_vs_docking.paths import DEFAULT_GENERATED_ROOT


ROOT = Path(__file__).resolve().parents[1]


def test_default_cli_uses_manuscript_filename(tmp_path, monkeypatch):
    captured = {}
    monkeypatch.setattr(sys, "argv", ["plot", "--root", str(tmp_path)])
    monkeypatch.setattr(plots, "load_experiment_points", lambda *args, **kwargs: None)
    monkeypatch.setattr(plots, "plot_points", lambda frame, path: captured.update(path=path))
    monkeypatch.setattr(plots, "build_resume_table", lambda *args, **kwargs: None)
    monkeypatch.setattr(plots, "build_resume_table_from_comparisons", lambda **kwargs: None)
    plots.main()
    assert captured["path"] == DEFAULT_GENERATED_ROOT / "green_vs_docking_summary_blocks.pdf"




@pytest.fixture
def comparison_root(tmp_path):
    experiment = tmp_path / "pesticide_9KQ3_blocks"
    experiment.mkdir()
    rows = []
    for seed, score, block in [(0, 0.2, True), (0, 0.8, False), (1, 0.4, True),
                                (1, 0.6, True), (1, 0.9, False)]:
        row = dict.fromkeys(statistics.POOLED_RESUME_METRICS, score)
        row.update(label=experiment.name, seed=seed, uses_block=block,
                   quickvina2_gpu_raw_values=-5 - 2 * score)
        rows.append(row)
    frame = pd.DataFrame(rows)
    frame.to_csv(experiment / "green_vs_quickvina.csv", index=False)
    frame.loc[frame.uses_block].to_csv(experiment / "green_vs_quickvina_blocks.csv", index=False)
    return tmp_path


@pytest.mark.parametrize("seed_first", [False, True])
@pytest.mark.parametrize("filtered", [False, True])
def test_statistics_options_apply_to_both_plot_and_resume(
    comparison_root, tmp_path, monkeypatch, seed_first, filtered
):
    captured = {}
    destination = tmp_path / "result.csv"
    arguments = ["plot", "--root", str(comparison_root),
                 "--output", str(tmp_path / "result.pdf"), "--resume-output", str(destination)]
    if seed_first:
        arguments.append("--aggregate-by-seed")
    if not filtered:
        arguments.append("--no-prefer-blocks-file")
    monkeypatch.setattr(sys, "argv", arguments)
    monkeypatch.setattr(plots, "plot_points", lambda frame, path: captured.update(points=frame))
    plots.main()

    values = [0.2, 0.4, 0.6] if filtered else [0.2, 0.8, 0.4, 0.6, 0.9]
    if seed_first:
        values = [0.2, 0.5] if filtered else [0.5, (0.4 + 0.6 + 0.9) / 3]
    mean, std = np.mean(values), np.std(values, ddof=1)
    resume = pd.read_csv(destination).iloc[0]
    point = captured["points"].iloc[0]
    assert resume.green_score == pytest.approx(mean)
    assert resume.green_score_std == pytest.approx(std)
    assert point.green_score_mean == pytest.approx(mean)
    assert point.green_score_std == pytest.approx(std)
    assert point.docking_mean == pytest.approx(-5 - 2 * mean)
    assert point.docking_std == pytest.approx(2 * std)
    assert resume.type == ("average across seed means" if seed_first else "pooled molecule average")
    assert resume.num_molecules == ((1.5 if filtered else 2.5) if seed_first else (3 if filtered else 5))


def test_historical_summary_export_is_an_explicit_option(tmp_path, monkeypatch):
    captured = {}
    destination = tmp_path / "historical.csv"
    monkeypatch.setattr(sys, "argv", ["plot", "--root", str(tmp_path),
                                      "--resume-output", str(destination), "--resume-from-summaries"])
    monkeypatch.setattr(plots, "load_experiment_points", lambda *args, **kwargs: None)
    monkeypatch.setattr(plots, "plot_points", lambda *args: None)
    monkeypatch.setattr(plots, "build_resume_table", lambda root, path: captured.update(path=path))
    monkeypatch.setattr(plots, "build_resume_table_from_comparisons", lambda **kw: pytest.fail("wrong exporter"))
    plots.main()
    assert captured["path"] == destination


def test_summary_export_cannot_silently_override_seed_aggregation(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["plot", "--aggregate-by-seed",
                                      "--resume-from-summaries"])
    with pytest.raises(SystemExit) as error:
        plots.main()
    assert error.value.code == 2
