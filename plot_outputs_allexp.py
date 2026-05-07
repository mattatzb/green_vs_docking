#!/usr/bin/env python3

"""Aggregate existing experiment outputs into one summary plot and resume CSV."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

# Use plasma palette for all plots by default.
_PLASMA = plt.cm.get_cmap("plasma")
plt.rcParams["image.cmap"] = "plasma"
plt.rcParams["axes.prop_cycle"] = plt.cycler(
    color=_PLASMA(np.linspace(0.1, 0.9, 10))
)
plt.rcParams["xtick.labelsize"] = 18
plt.rcParams["ytick.labelsize"] = 18

RESUME_COLUMNS = [
    "num_molecules",
    "green_score",
    "green_score_std",
    "atom_economy",
    "atom_economy_std",
    "solvent_score",
    "solvent_score_std",
    "temperature_score",
    "temperature_score_std",
    "num_steps",
    "num_steps_std",
    "biomass_utilization",
    "biomass_utilization_std",
    "block_usage_pct"
]

POOLED_RESUME_METRICS = [
    "green_score",
    "atom_economy",
    "solvent_score",
    "temperature_score",
    "num_steps",
    "biomass_utilization",
]

def _clean_label(raw: str) -> str:
    """ Cleans an experiment label by splitting on colons and taking the last non-empty segment, or returning the original string if it cannot be processed. This is useful for normalizing labels that may contain hierarchical information separated by colons. """
    if not isinstance(raw, str):
        return str(raw)
    parts = [segment.strip() for segment in raw.split(":") if segment.strip()]
    return parts[-1] if parts else raw.strip()

BLOCK_KEYWORDS = ("blocks",)

def _comparison_csv_for_experiment(experiment_dir: Path, prefer_block_file: bool) -> Path | None:
    """ 
    Determines the appropriate comparison CSV file for an experiment directory based on its name and the prefer_block_file flag. 
    If prefer_block_file is True and the experiment directory name contains any of the BLOCK_KEYWORDS, it will look for a "green_vs_quickvina_blocks.csv" file. Otherwise, it will look for a "green_vs_quickvina.csv" file. If the chosen file exists, its path is returned; otherwise, None is returned. 
    """
    default_path = experiment_dir / "green_vs_quickvina.csv"
    blocks_path = experiment_dir / "green_vs_quickvina_blocks.csv"

    # Check if we should prefer the blocks file based on the presence of block-related keywords in the experiment directory name and the existence of the blocks file. 
    # If so, we choose the blocks file; otherwise, we choose the default comparison file. We then check if the chosen file exists and return its path if it does, or None if it does not.
    use_blocks = (
        prefer_block_file
        and any(keyword in experiment_dir.name.lower() for keyword in BLOCK_KEYWORDS)
        and blocks_path.exists()
    )
    chosen = blocks_path if use_blocks else default_path
    return chosen if chosen.exists() else None

def aggregate_label_stats(comparison_df: pd.DataFrame, aggregate_by_seed: bool) -> pd.DataFrame:
    """ 
    Aggregates the comparison DataFrame by experiment label, optionally first averaging within each seed. 
    """

    # If aggregate_by_seed is True and the DataFrame contains a "seed" column, we first group by both "normalized_label" and "seed" to compute the mean green score and docking reward for each seed. We then group by "normalized_label" again to compute the overall mean and standard deviation across seeds for each label.
    if aggregate_by_seed and "seed" in comparison_df.columns:
        per_seed = (
            comparison_df.groupby(["normalized_label", "seed"])
            .agg(
                green_score_mean=("green_score", "mean"),
                docking_mean=("quickvina2_gpu_raw_values", "mean"),
            )
            .reset_index()
        )
        aggregated = (
            per_seed.groupby("normalized_label")
            .agg(
                green_score_mean=("green_score_mean", "mean"),
                green_score_std=("green_score_mean", "std"),
                docking_mean=("docking_mean", "mean"),
                docking_std=("docking_mean", "std"),
            )
            .reset_index()
        )
    
    # If aggregate_by_seed is False or there is no "seed" column, we directly group by "normalized_label" to compute the mean and standard deviation of the green score and docking reward across all data points for each label.
    else:
        aggregated = (
            comparison_df.groupby("normalized_label")[["green_score", "quickvina2_gpu_raw_values"]]
            .agg(["mean", "std"])
            .reset_index()
        )
        aggregated.columns = [
            "normalized_label",
            "green_score_mean",
            "green_score_std",
            "docking_mean",
            "docking_std",
        ]

    aggregated = aggregated.rename(columns={"normalized_label": "label"})
    return aggregated


def load_experiment_points(root: Path, aggregate_by_seed: bool, prefer_block_file: bool) -> pd.DataFrame:
    """
    Loads and aggregates experiment data from the specified root directory. For each experiment subdirectory, it looks for a comparison CSV file (either "green_vs_quickvina.csv" or "green_vs_quickvina_blocks.csv" based on the prefer_block_file flag and the presence of block-related keywords in the directory name). 
    It reads the comparison data, normalizes the labels, and aggregates the green score and docking reward statistics by label (optionally first averaging within each seed).
    The resulting DataFrame contains one row per unique label with mean and standard deviation values for both metrics, along with the experiment directory name for reference.
    """

    rows: list[dict[str, object]] = []
    baseline_seen = False
    
    # We iterate through each subdirectory in the root directory, checking if it is a directory. For each valid experiment directory, we determine the appropriate comparison CSV file to load based on the prefer_block_file flag and the presence of block-related keywords in the directory name. 
    for experiment_dir in sorted(root.iterdir()):
        if not experiment_dir.is_dir():
            continue
        comparison_path = _comparison_csv_for_experiment(experiment_dir, prefer_block_file)
        if comparison_path is None:
            continue


        comparison_df = pd.read_csv(comparison_path)
        if comparison_df.empty:
            continue
        if not {"label", "green_score", "quickvina2_gpu_raw_values"}.issubset(comparison_df.columns):
            continue
        comparison_df = comparison_df.copy()
        comparison_df["normalized_label"] = comparison_df["label"].map(_clean_label)

        grouped = aggregate_label_stats(comparison_df, aggregate_by_seed)
        grouped["experiment_dir"] = experiment_dir.name
        for record in grouped.to_dict(orient="records"):
            record["label"] = _clean_label(record["label"])
            label_lower = record["label"].lower()
            if label_lower == "baseline":
                if baseline_seen:
                    continue
                baseline_seen = True
            rows.append(record)

    return pd.DataFrame(rows)


def load_experiment_rows(root: Path, prefer_block_file: bool) -> pd.DataFrame:
    """
    Load molecule-level rows from the same comparison CSVs used by the plots.

    Baseline is present in every experiment folder, so it is kept only once to
    avoid giving it repeated weight in the combined KDE.
    """
    rows: list[pd.DataFrame] = []
    baseline_seen = False
    required_columns = {"label", "green_score", "quickvina2_gpu_raw_values"}

    for experiment_dir in sorted(root.iterdir()):
        if not experiment_dir.is_dir():
            continue
        comparison_path = _comparison_csv_for_experiment(experiment_dir, prefer_block_file)
        if comparison_path is None:
            continue

        comparison_df = pd.read_csv(comparison_path)
        if comparison_df.empty or not required_columns.issubset(comparison_df.columns):
            continue

        comparison_df = comparison_df.copy()
        comparison_df["normalized_label"] = comparison_df["label"].map(_clean_label)
        for label, label_df in comparison_df.groupby("normalized_label", sort=False):
            label = _clean_label(label)
            if label.lower() == "baseline":
                if baseline_seen:
                    continue
                baseline_seen = True

            selected = label_df.copy()
            selected["label"] = label
            selected["experiment_dir"] = experiment_dir.name
            rows.append(selected)

    if not rows:
        return pd.DataFrame(columns=["label", "green_score", "quickvina2_gpu_raw_values"])
    return pd.concat(rows, ignore_index=True)


def sample_rows_by_label(df: pd.DataFrame, max_rows_per_label: int | None) -> pd.DataFrame:
    """Optionally down-sample each experiment label to speed up KDE rendering."""
    if max_rows_per_label is None or max_rows_per_label <= 0 or df.empty:
        return df

    sampled: list[pd.DataFrame] = []
    for _, label_df in df.groupby("label", sort=False):
        if len(label_df) > max_rows_per_label:
            sampled.append(label_df.sample(n=max_rows_per_label, random_state=0))
        else:
            sampled.append(label_df)
    return pd.concat(sampled, ignore_index=True)


def plot_points(df: pd.DataFrame, output_path: Path) -> None:
    """Save a scatter/error-bar plot of mean green score versus docking."""
    if df.empty:
        raise RuntimeError("No experiment data found under the specified directory.")

    df["green_score_std"] = df["green_score_std"].fillna(0.0)
    df["docking_std"] = df["docking_std"].fillna(0.0).abs()

    fig, ax = plt.subplots(figsize=(9, 6))
    handles = []
    labels = []

    for _, row in df.iterrows():
        is_baseline = row["label"].lower() == "baseline"
        marker = "s" if is_baseline else "o"
        label = row["label"]
        container = ax.errorbar(
            row["docking_mean"],
            row["green_score_mean"],
            xerr=row["docking_std"],
            yerr=row["green_score_std"],
            fmt=marker,
            capsize=4,
            label=label,
        )
        handles.append(container)
        labels.append(label)

    ax.set_xlabel("QuickVina2 GPU Reward (mean ± std)")
    ax.set_ylabel("Green Score (mean ± std)")
    ax.set_title("Average Green Score vs Docking")
    ax.invert_xaxis()
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    ax.legend(handles, labels, fontsize=8, loc="best", frameon=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_joint_kde(
    df: pd.DataFrame,
    output_path: Path,
    docking_min: float = -14.0,
    docking_max: float = 0.0,
) -> None:
    """Save one joint KDE overlay containing all experiment molecule clouds."""
    try:
        import seaborn as sns
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "The joint KDE plot requires seaborn in the active Python environment."
        ) from exc

    required_columns = {"label", "green_score", "quickvina2_gpu_raw_values"}
    if df.empty or not required_columns.issubset(df.columns):
        raise RuntimeError("No molecule-level comparison data found for the joint KDE plot.")

    plot_df = df.dropna(subset=["label", "green_score", "quickvina2_gpu_raw_values"]).copy()
    if plot_df.empty:
        raise RuntimeError("All molecule-level rows have missing green or docking values.")

    # Negating the docking score preserves the desired visual convention:
    # more negative docking values appear on the right, while tick labels show
    # the original docking score.
    plot_df["_quickvina_inverted"] = -plot_df["quickvina2_gpu_raw_values"]
    labels = list(dict.fromkeys(plot_df["label"].astype(str)))
    colors = sns.color_palette("plasma", n_colors=len(labels))

    grid = sns.JointGrid(
        data=plot_df,
        x="_quickvina_inverted",
        y="green_score",
        height=8,
        ratio=5,
        space=0.08,
    )

    legend_handles: list[Line2D] = []
    for label, color in zip(labels, colors):
        label_df = plot_df.loc[plot_df["label"] == label]
        if len(label_df) < 2:
            continue

        sns.kdeplot(
            data=label_df,
            x="_quickvina_inverted",
            y="green_score",
            ax=grid.ax_joint,
            levels=8,
            thresh=0.05,
            fill=False,
            linewidths=0.8,
            color=color,
            warn_singular=False,
        )
        sns.kdeplot(
            data=label_df,
            x="_quickvina_inverted",
            ax=grid.ax_marg_x,
            color=color,
            linewidth=1.0,
            fill=False,
            warn_singular=False,
        )
        sns.kdeplot(
            data=label_df,
            y="green_score",
            ax=grid.ax_marg_y,
            color=color,
            linewidth=1.0,
            fill=False,
            warn_singular=False,
        )
        legend_handles.append(Line2D([0], [0], color=color, lw=1.5, label=label))

    if not legend_handles:
        raise RuntimeError("Not enough rows per experiment to draw a joint KDE plot.")

    x_min = min(docking_min, docking_max)
    x_max = max(docking_min, docking_max)
    grid.ax_joint.set_xlim(-x_max, -x_min)
    grid.ax_joint.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{-value:g}"))
    grid.ax_marg_x.set_xlim(grid.ax_joint.get_xlim())

    y_min = plot_df["green_score"].min()
    y_max = plot_df["green_score"].max()
    y_pad = max((y_max - y_min) * 0.05, 0.02)
    grid.ax_joint.set_ylim(max(0.0, y_min - y_pad), min(1.0, y_max + y_pad))
    grid.ax_marg_y.set_ylim(grid.ax_joint.get_ylim())

    grid.ax_joint.set_xlabel("QuickVina2 GPU Raw Value")
    grid.ax_joint.set_ylabel("Green Score")
    grid.ax_joint.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
    grid.ax_joint.legend(
        handles=legend_handles,
        fontsize=7,
        loc="lower left",
        frameon=True,
        title="Experiment",
    )
    grid.figure.suptitle("Green Score vs Docking KDE Across Experiments", y=1.02)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    grid.figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(grid.figure)


def build_resume_table(root: Path, resume_path: Path) -> Path | None:
    """
    Combine experiment-average rows from all output folders into one CSV.

    Both the standard summary and the optional block-only summary are included
    when present. Baseline rows from block-only summaries are skipped to avoid
    duplicate baseline entries.
    """
    rows: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()

    def append_rows(summary_path: Path, type_label: str, *, skip_baseline: bool = False) -> None:
        if not summary_path.exists():
            return
        df = pd.read_csv(summary_path)
        if df.empty or "summary_level" not in df.columns:
            return
        avg_rows = df.loc[df["summary_level"] == "experiment_average"]
        if avg_rows.empty:
            return
        for _, row in avg_rows.iterrows():
            raw_label = row.get("label", "")
            label = _clean_label(raw_label)
            if not label:
                continue
            if skip_baseline and label.lower() == "baseline":
                continue
            key = (type_label, label.lower())
            if key in seen:
                continue
            entry = {
                "label": label,
                "type": type_label,
            }
            for col in RESUME_COLUMNS:
                entry[col] = row.get(col)
            rows.append(entry)
            seen.add(key)

    for experiment_dir in sorted(root.iterdir()):
        if not experiment_dir.is_dir():
            continue
        append_rows(experiment_dir / "summary_scores.csv", "per seed average")
        append_rows(
            experiment_dir / "summary_scores_blocks_only.csv",
            "per seed average (only blocks)",
            skip_baseline=True,
        )

    if not rows:
        return None

    resume_df = pd.DataFrame(rows)
    resume_df.to_csv(resume_path, index=False)
    return resume_path


def _summarize_comparison_rows(
    df: pd.DataFrame,
    label: str,
    aggregate_by_seed: bool,
) -> dict[str, object]:
    if aggregate_by_seed and "seed" in df.columns:
        seed_metric_means = df.groupby("seed")[POOLED_RESUME_METRICS].mean()
        summary: dict[str, object] = {
            "label": label,
            "type": "average across seed means",
            "num_molecules": df.groupby("seed").size().mean(),
        }
        for metric in POOLED_RESUME_METRICS:
            summary[metric] = seed_metric_means[metric].mean()
            summary[f"{metric}_std"] = seed_metric_means[metric].std()
        if "uses_block" in df.columns:
            summary["block_usage_pct"] = df.groupby("seed")["uses_block"].mean().mean() * 100
        else:
            summary["block_usage_pct"] = pd.NA
        return summary

    summary = {
        "label": label,
        "type": "pooled molecule average",
        "num_molecules": len(df),
    }
    for metric in POOLED_RESUME_METRICS:
        summary[metric] = df[metric].mean()
        summary[f"{metric}_std"] = df[metric].std()
    if "uses_block" in df.columns:
        summary["block_usage_pct"] = df["uses_block"].mean() * 100
    else:
        summary["block_usage_pct"] = pd.NA
    return summary


def build_resume_table_from_comparisons(
    root: Path,
    resume_path: Path,
    prefer_block_file: bool,
    aggregate_by_seed: bool,
) -> Path | None:
    """
    Build a resume CSV from the same comparison files used for plotting.

    This mode keeps the resume statistics consistent with `plot_points`: block
    folders can use `green_vs_quickvina_blocks.csv`, and the standard deviation
    is either molecule-pooled or seed-first depending on `aggregate_by_seed`.
    """
    rows: list[dict[str, object]] = []
    baseline_seen = False
    required_columns = {"label", *POOLED_RESUME_METRICS}

    for experiment_dir in sorted(root.iterdir()):
        if not experiment_dir.is_dir():
            continue
        comparison_path = _comparison_csv_for_experiment(experiment_dir, prefer_block_file)
        if comparison_path is None:
            continue
        comparison_df = pd.read_csv(comparison_path)
        if comparison_df.empty or not required_columns.issubset(comparison_df.columns):
            continue

        comparison_df = comparison_df.copy()
        comparison_df["normalized_label"] = comparison_df["label"].map(_clean_label)
        for label, label_df in comparison_df.groupby("normalized_label", sort=False):
            label = _clean_label(label)
            if label.lower() == "baseline":
                if baseline_seen:
                    continue
                baseline_seen = True
            rows.append(_summarize_comparison_rows(label_df, label, aggregate_by_seed))

    if not rows:
        return None

    resume_df = pd.DataFrame(rows)
    final_columns = ["label", "type", *RESUME_COLUMNS]
    resume_df = resume_df.reindex(columns=final_columns)
    resume_df.to_csv(resume_path, index=False)
    return resume_path


def main() -> None:
    """Parse CLI options, create the aggregate plot, and write `resume.csv`."""
    parser = argparse.ArgumentParser(description="Plot mean green score vs docking for outputs.")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("/work/liac/tatzber/green_vs_docking/outputs/outputs_bae"),
        help="Directory that contains experiment subfolders.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("/work/liac/tatzber/green_vs_docking/outputs/outputs_bae/green_vs_docking_summary.png"),
        help="Path of the PNG file to create.",
    )
    parser.add_argument(
        "--resume-output",
        type=Path,
        default=None,
        help="Optional CSV path for an aggregated summary (defaults to <root>/resume.csv).",
    )
    parser.add_argument(
        "--aggregate-by-seed",
        action="store_true",
        help="Average docking/green per seed first, then aggregate across seeds (std reflects variability between seeds).",
    )
    parser.add_argument(
        "--prefer-blocks-file",
        action="store_true",
        help="For experiments whose folder names contain 'blocks', load green_vs_quickvina_blocks.csv instead of the default comparison file.",
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
        help="Only create the joint KDE plot; skip the average scatter plot and resume CSV.",
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
    args = parser.parse_args()

    root_dir = args.root.expanduser().resolve()
    plot_path = args.output.expanduser().resolve()
    parser_default_output = parser.get_default("output").expanduser().resolve()
    if args.prefer_blocks_file and plot_path == parser_default_output:
        plot_path = plot_path.with_name(f"{plot_path.stem}_blocks{plot_path.suffix}")

    resume_path = args.resume_output
    if resume_path is None:
        resume_path = root_dir / "resume.csv"
    else:
        resume_path = resume_path.expanduser().resolve()

    if args.joint_kde_output is not None or args.joint_kde_only:
        if args.joint_kde_output is None:
            suffix = "_blocks" if args.prefer_blocks_file else ""
            joint_kde_path = root_dir / f"green_vs_docking_joint_kde{suffix}.png"
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

    if args.prefer_blocks_file:
        resume_file = build_resume_table_from_comparisons(
            root=root_dir,
            resume_path=resume_path,
            prefer_block_file=args.prefer_blocks_file,
            aggregate_by_seed=args.aggregate_by_seed,
        )
    else:
        resume_file = build_resume_table(root_dir, resume_path)
    if resume_file is not None:
        print(f"Wrote resume table to: {resume_file}")


if __name__ == "__main__":
    main()
