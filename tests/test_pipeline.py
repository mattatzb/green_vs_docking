"""Small orchestration tests that do not depend on manuscript datasets."""

from pathlib import Path
import json

import pandas as pd
import pytest


def test_pipeline_runs_every_condition_then_summarizes(tmp_path, monkeypatch):
    from green_vs_docking import pipeline

    data = tmp_path / "data"
    baseline = "pesticide_9KQ3_basic"
    names = ["pesticide_9KQ3_blocks", "pesticide_9KQ3_solvent"]
    for name in [baseline, *names]:
        (data / name).mkdir(parents=True)
    output = tmp_path / "result"
    calls = []
    monkeypatch.setattr(pipeline, "run_analysis", lambda **kw: calls.append(("score", kw)))
    monkeypatch.setattr(pipeline, "load_experiment_points", lambda *a, **kw: calls.append(("load", a, kw)))
    monkeypatch.setattr(pipeline, "plot_points", lambda *a: calls.append(("plot", a)))
    monkeypatch.setattr(pipeline, "build_resume_table_from_comparisons", lambda **kw: calls.append(("resume", kw)))
    monkeypatch.setattr(pipeline, "write_table_summary", lambda **kw: calls.append(("table", kw)))
    artifacts = pipeline.run_pipeline(data_root=data, output_dir=output)
    assert [call[0] for call in calls] == ["score", "score", "load", "plot", "resume", "table"]
    for call, name in zip(calls, names):
        assert call[1]["experiments"] == {"baseline": data / baseline, name: data / name}
        assert call[1]["output_dir"] == output / name
    assert calls[-1][1]["oracle_root"] == data
    assert calls[-1][1]["outputs_root"] == output
    assert artifacts["plot"] == output / "green_vs_docking_summary_blocks.pdf"
    assert artifacts["resume"] == output / "resume.csv"
    assert artifacts["table"] == output / "latex_table_summary.csv"


def test_existing_comparisons_are_read_without_rescoring(tmp_path, monkeypatch):
    from green_vs_docking import pipeline

    source = tmp_path / "saved"
    source.mkdir()
    destination = tmp_path / "fresh"
    captured = {}
    monkeypatch.setattr(pipeline, "run_analysis", lambda **kw: pytest.fail("unexpected rescoring"))
    monkeypatch.setattr(pipeline, "load_experiment_points", lambda root, **kw: captured.update(source=root))
    monkeypatch.setattr(pipeline, "plot_points", lambda *a: None)
    monkeypatch.setattr(pipeline, "build_resume_table_from_comparisons", lambda **kw: captured.update(resume=kw))
    monkeypatch.setattr(pipeline, "write_table_summary", lambda **kw: captured.update(table=kw))
    pipeline.run_pipeline(data_root=tmp_path, output_dir=destination, comparison_root=source)
    assert captured["source"] == source
    assert captured["table"]["outputs_root"] == source
    assert captured["table"]["output_path"].parent == destination
    assert captured["resume"]["resume_path"].parent == destination
    assert not list(source.iterdir())


def test_pipeline_requires_conditions_for_raw_scoring(tmp_path):
    from green_vs_docking import pipeline

    with pytest.raises(ValueError, match="condition"):
        pipeline.run_pipeline(data_root=tmp_path, output_dir=tmp_path / "output")


def test_raw_pipeline_preserves_pre_refactor_csvs(tiny_data_root, tmp_path, monkeypatch):
    from green_vs_docking import pipeline, workflow

    monkeypatch.setattr(workflow, "save_green_vs_quickvina_plot", lambda *a, **kw: None)
    output = tmp_path / "fresh"
    artifacts = pipeline.run_pipeline(data_root=tiny_data_root, output_dir=output)
    expected = json.loads((Path(__file__).parent / "fixtures/pre_refactor_scores.json").read_text())
    for name, frame in expected.items():
        reference = pd.DataFrame(frame["data"], columns=frame["columns"])
        actual = pd.read_csv(output / "pesticide_9KQ3_solvent" / name)
        pd.testing.assert_frame_equal(actual, reference, check_dtype=False, rtol=1e-12, atol=1e-12)
    table = pd.read_csv(artifacts["table"])
    assert table.solved_mean.tolist() == [2.0, 2.0]
    assert table.non_solved_mean.tolist() == [1.0, 1.0]
    assert table.available_replicates.tolist() == [2, 2]
    assert table.pooled_solved_molecules.tolist() == [4, 4]
    assert artifacts["plot"].read_bytes().startswith(b"%PDF-")
    assert pd.read_csv(artifacts["resume"]).num_molecules.tolist() == [4, 4]


def test_kde_style_does_not_leak_into_summary(tmp_path, monkeypatch):
    from green_vs_docking import pipeline

    expected = pipeline.plt.rcParams["xtick.labelsize"]
    def score(**kwargs):
        pipeline.plt.rcParams["xtick.labelsize"] = 42
    def draw(*args):
        assert pipeline.plt.rcParams["xtick.labelsize"] == expected
    monkeypatch.setattr(pipeline, "run_analysis", score)
    monkeypatch.setattr(pipeline, "load_experiment_points", lambda *a, **kw: None)
    monkeypatch.setattr(pipeline, "plot_points", draw)
    monkeypatch.setattr(pipeline, "build_resume_table_from_comparisons", lambda **kw: None)
    monkeypatch.setattr(pipeline, "write_table_summary", lambda **kw: None)
    pipeline.run_pipeline(data_root=tmp_path, output_dir=tmp_path / "out", conditions=["condition"])
    assert pipeline.plt.rcParams["xtick.labelsize"] == expected
