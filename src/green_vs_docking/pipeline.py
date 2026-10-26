"""The manuscript workflow: route scoring, aggregate plot/resume, final table."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence
import matplotlib.pyplot as plt
from .analysis.summary import build_resume_table, build_resume_table_from_comparisons, load_experiment_points
from .analysis.table_summary import DEFAULT_BASELINE_EXPERIMENT, DEFAULT_SEEDS, write_table_summary
from .paths import DEFAULT_DATA_ROOT, DEFAULT_GENERATED_ROOT, DEFAULT_SOLVENT_GUIDE_PATH, DEFAULT_BLOCK_FILE_PATH
from .plots.summary import plot_points
from .workflow import run_analysis


def run_pipeline(
    data_root: Path = DEFAULT_DATA_ROOT,
    output_dir: Path = DEFAULT_GENERATED_ROOT,
    comparison_root: Path | None = None,
    baseline_experiment: str = DEFAULT_BASELINE_EXPERIMENT,
    conditions: Sequence[str] | None = None,
    solvent_guide: Path = DEFAULT_SOLVENT_GUIDE_PATH,
    block_file: Path = DEFAULT_BLOCK_FILE_PATH,
    namerxn_binary: Path | None = None,
    aggregate_by_seed: bool = False,
    prefer_block_file: bool = True,
    resume_from_summaries: bool = False,
) -> dict[str, Path]:
    """Write fresh artifacts while reading the original inputs without mutation."""
    if aggregate_by_seed and resume_from_summaries:
        raise ValueError("Choose seed aggregation or the historical summary export.")
    data_root = data_root.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    if comparison_root is None:
        if conditions is None:
            conditions = [p.name for p in sorted(data_root.glob("pesticide_9KQ3_*"))
                          if p.is_dir() and p.name != baseline_experiment]
        if not conditions:
            raise ValueError(f"No condition experiments found under {data_root}")
        for condition in conditions:
            # Old scoring commands ran separately; KDE themes must not leak
            # between stages when the same work is done in one process.
            with plt.rc_context():
                run_analysis(
                    experiments={"baseline": data_root / baseline_experiment,
                                 condition: data_root / condition},
                    output_dir=output_dir / condition,
                    solvent_guide=solvent_guide, block_file=block_file,
                    namerxn_binary=namerxn_binary,
                )
        comparison_root = output_dir
    else:
        comparison_root = comparison_root.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = "_blocks" if prefer_block_file else ""
    artifacts = {"plot": output_dir / f"green_vs_docking_summary{suffix}.pdf",
                 "resume": output_dir / "resume.csv",
                 "table": output_dir / "latex_table_summary.csv"}
    points = load_experiment_points(
        comparison_root, aggregate_by_seed=aggregate_by_seed, prefer_block_file=prefer_block_file
    )
    plot_points(points, artifacts["plot"])
    if resume_from_summaries:
        build_resume_table(comparison_root, artifacts["resume"])
    else:
        build_resume_table_from_comparisons(
            root=comparison_root, resume_path=artifacts["resume"],
            prefer_block_file=prefer_block_file, aggregate_by_seed=aggregate_by_seed,
        )
    write_table_summary(
        outputs_root=comparison_root, oracle_root=data_root,
        output_path=artifacts["table"], baseline_experiment=baseline_experiment,
        seeds=DEFAULT_SEEDS,
    )
    return artifacts
