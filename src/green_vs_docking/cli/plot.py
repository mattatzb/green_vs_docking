"""Plot and resume CLI; manuscript statistics are the defaults."""

from __future__ import annotations

import argparse
from pathlib import Path
from ..paths import DEFAULT_GENERATED_ROOT, REFERENCE_OUTPUTS_ROOT
from ..analysis.summary import (load_experiment_points, load_experiment_rows, sample_rows_by_label, build_resume_table, build_resume_table_from_comparisons)
from ..plots.summary import plot_points, plot_joint_kde

def main(argv: list[str] | None = None) -> None:
    """Parse CLI options, create the aggregate plot, and write `resume.csv`."""
    parser = argparse.ArgumentParser(
        description="Reproduce the manuscript plot and resume CSV using block-filtered, pooled molecule statistics."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=REFERENCE_OUTPUTS_ROOT,
        help="Directory that contains experiment subfolders.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="PDF destination (default: outputs/9KQ3/green_vs_docking_summary_blocks.pdf in the repository; omit _blocks with --no-prefer-blocks-file).",
    )
    parser.add_argument(
        "--resume-output",
        type=Path,
        default=None,
        help="CSV destination (default: outputs/9KQ3/resume.csv in the repository).",
    )
    statistics_group = parser.add_mutually_exclusive_group()
    statistics_group.add_argument(
        "--aggregate-by-seed",
        action="store_true",
        help="Use mean/std across seed means for both plot and resume instead of the default pooled molecule statistics.",
    )
    statistics_group.add_argument(
        "--resume-from-summaries",
        action="store_true",
        help="Export historical per-seed summary-file averages (both full and block-only rows) to the resume CSV; the plot retains pooled statistics.",
    )
    parser.add_argument(
        "--prefer-blocks-file",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use block-filtered comparison files for block experiments (default); --no-prefer-blocks-file selects unfiltered comparisons.",
    )
    parser.add_argument(
        "--joint-kde-output",
        type=Path,
        default=None,
        help="Optional PNG path for one joint KDE overlay containing all experiment molecule clouds.",
    )
    parser.add_argument(
        "--joint-kde-sample",
        type=int,
        default=None,
        help="Optional maximum number of molecules sampled per experiment label for faster KDE rendering.",
    )
    parser.add_argument(
        "--joint-kde-only",
        action="store_true",
        help="Only create the joint KDE plot in outputs/9KQ3 (or --joint-kde-output); skip the average scatter plot and resume CSV.",
    )
    parser.add_argument(
        "--joint-kde-docking-min",
        type=float,
        default=-14.0,
        help="Minimum docking value shown on the joint KDE x-axis.",
    )
    parser.add_argument(
        "--joint-kde-docking-max",
        type=float,
        default=0.0,
        help="Maximum docking value shown on the joint KDE x-axis.",
    )
    args = parser.parse_args(argv)

    root_dir = args.root.expanduser().resolve()
    if args.output is None:
        suffix = "_blocks" if args.prefer_blocks_file else ""
        plot_path = DEFAULT_GENERATED_ROOT / f"green_vs_docking_summary{suffix}.pdf"
    else:
        plot_path = args.output.expanduser().resolve()

    resume_path = args.resume_output
    if resume_path is None:
        resume_path = DEFAULT_GENERATED_ROOT / "resume.csv"
    else:
        resume_path = resume_path.expanduser().resolve()

    if args.joint_kde_output is not None or args.joint_kde_only:
        if args.joint_kde_output is None:
            suffix = "_blocks" if args.prefer_blocks_file else ""
            joint_kde_path = DEFAULT_GENERATED_ROOT / f"green_vs_docking_joint_kde{suffix}.png"
        else:
            joint_kde_path = args.joint_kde_output.expanduser().resolve()
        molecule_rows = load_experiment_rows(root_dir, prefer_block_file=args.prefer_blocks_file)
        molecule_rows = sample_rows_by_label(molecule_rows, args.joint_kde_sample)
        plot_joint_kde(
            molecule_rows,
            joint_kde_path,
            docking_min=args.joint_kde_docking_min,
            docking_max=args.joint_kde_docking_max,
        )
        print(f"Wrote joint KDE plot to: {joint_kde_path}")

    if args.joint_kde_only:
        return

    df = load_experiment_points(root_dir, aggregate_by_seed=args.aggregate_by_seed, prefer_block_file=args.prefer_blocks_file)
    plot_points(df, plot_path)
    print(f"Wrote summary plot to: {plot_path}")

    if args.resume_from_summaries:
        resume_file = build_resume_table(root_dir, resume_path)
    else:
        resume_file = build_resume_table_from_comparisons(
            root=root_dir,
            resume_path=resume_path,
            prefer_block_file=args.prefer_blocks_file,
            aggregate_by_seed=args.aggregate_by_seed,
        )
    if resume_file is not None:
        print(f"Wrote resume table to: {resume_file}")


