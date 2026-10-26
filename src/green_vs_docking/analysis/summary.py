"""Comparison selection and experiment summary statistics."""
from __future__ import annotations
from pathlib import Path
import pandas as pd

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


BLOCK_KEYWORDS = ("blocks",)


def _clean_label(raw: str) -> str:
    """ Cleans an experiment label by splitting on colons and taking the last non-empty segment, or returning the original string if it cannot be processed. This is useful for normalizing labels that may contain hierarchical information separated by colons. """
    if not isinstance(raw, str):
        return str(raw)
    parts = [segment.strip() for segment in raw.split(":") if segment.strip()]
    return parts[-1] if parts else raw.strip()


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


def load_experiment_points(
    root: Path, aggregate_by_seed: bool = False, prefer_block_file: bool = True
) -> pd.DataFrame:
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


def load_experiment_rows(root: Path, prefer_block_file: bool = True) -> pd.DataFrame:
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
    prefer_block_file: bool = True,
    aggregate_by_seed: bool = False,
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


