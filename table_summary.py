"""Build compact CSV inputs for the 9KQ3 LaTeX results table."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from functools import lru_cache
from typing import Iterable
from rdkit.Chem.Scaffolds import MurckoScaffold

import pandas as pd
from rdkit import Chem
from rdkit.Chem import QED

DEFAULT_OUTPUTS_ROOT = Path("/work/liac/tatzber/green_vs_docking/outputs/9KQ3")
DEFAULT_ORACLE_ROOT = Path("/work/liac/vsabanzagil/condition_optimization")
DEFAULT_OUTPUT = DEFAULT_OUTPUTS_ROOT / "latex_table_summary.csv"
DEFAULT_BASELINE_EXPERIMENT = "pesticide_9KQ3_basic"
DEFAULT_SEEDS = tuple(range(5))


@dataclass(frozen=True)
class ExperimentSpec:
    """Mapping between one output folder and its original oracle experiment."""

    table_label: str
    output_label: str
    oracle_name: str
    output_dir: Path
    is_baseline: bool = False


@lru_cache(maxsize=None)
def canonicalize_smiles(smiles: str) -> str | None:
    """Return RDKit canonical SMILES, or None when parsing fails."""
    if not isinstance(smiles, str) or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


@lru_cache(maxsize=None)
def qed_from_smiles(smiles: str) -> float:
    """ Returns the QED score for a given SMILES string, or NaN if the SMILES string is invalid. """
    if not isinstance(smiles, str) or not smiles:
        return math.nan
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return math.nan
    return float(QED.qed(mol))


@lru_cache(maxsize=None)
def bms_from_smiles(smiles: str) -> str | None:
    """ Returns the Bemis-Murcko scaffold for a given SMILES string, or None if the SMILES string is invalid or has no scaffold. """
    if not isinstance(smiles, str) or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        return None
    return Chem.MolToSmiles(scaffold)


def sample_std(values: Iterable[float]) -> float:
    """ Returns the sample standard deviation of the given values, or NaN if there are fewer than 2 values. """
    series = pd.Series(list(values), dtype="float64")
    return float(series.std()) if len(series) > 1 else math.nan


def mean_value(values: Iterable[float]) -> float:
    """ Returns the mean of the given values, or NaN if there are no values. """
    series = pd.Series(list(values), dtype="float64")
    return float(series.mean()) if not series.empty else math.nan


def clean_table_label(name: str) -> str:
    """ Converts an experiment directory name into a more human-readable table label. """
    if name == "baseline":
        return "Baseline"
    label = name
    label = label.removeprefix("pesticide_9KQ3_")
    label = label.replace("no_deprotection", "No protection/deprotection")
    label = label.replace("solvent", "Solvent")
    label = label.replace("blocks", "Blocks")
    label = label.replace("_", " + ")
    return label


def discover_experiments(outputs_root: Path, baseline_experiment: str) -> list[ExperimentSpec]:
    """ Discovers experiment specifications by looking for green_vs_quickvina.csv files under the outputs root directory. The baseline experiment is expected to be present and will be included as the first entry in the returned list. """

    # The baseline experiment is included as a special case to ensure it is always present and listed first in the table summary, even if it does not have a green_vs_quickvina.csv file. 
    # This allows the baseline to be included in the summary statistics and comparisons, while still allowing for flexibility in how the baseline experiment is structured within the outputs directory.
    specs = [
        ExperimentSpec(
            table_label="Baseline",
            output_label="baseline",
            oracle_name=baseline_experiment,
            output_dir=next_output_dir_with_baseline(outputs_root),
            is_baseline=True,
        )
    ]
    
    # Discover additional experiments by looking for green_vs_quickvina.csv files under the outputs root directory. 
    # Each experiment is expected to be in its own subdirectory, and the presence of a green_vs_quickvina.csv file is used as an indicator that the subdirectory contains a valid experiment to include in the summary.
    # The table label for each experiment is derived from the subdirectory name using the clean_table_label function.
    for experiment_dir in sorted(path for path in outputs_root.iterdir() if path.is_dir()):
        comparison_path = experiment_dir / "green_vs_quickvina.csv"
        if not comparison_path.exists():
            continue
        specs.append(
            ExperimentSpec(
                table_label=clean_table_label(experiment_dir.name),
                output_label=experiment_dir.name,
                oracle_name=experiment_dir.name,
                output_dir=experiment_dir,
            )
        )
    return specs


def next_output_dir_with_baseline(outputs_root: Path) -> Path:
    """ Finds the next output directory under the outputs root that contains a green_vs_quickvina.csv file, which is expected to correspond to the baseline experiment. This is used to ensure that the baseline experiment is included in the summary statistics and comparisons, even if it does not have a green_vs_quickvina.csv file in its own subdirectory. """
    for experiment_dir in sorted(path for path in outputs_root.iterdir() if path.is_dir()):
        comparison_path = experiment_dir / "green_vs_quickvina.csv"
        if comparison_path.exists():
            return experiment_dir
    raise FileNotFoundError(f"No green_vs_quickvina.csv found under {outputs_root}")


def load_green_data(spec: ExperimentSpec) -> pd.DataFrame:
    """ Loads the green_vs_quickvina.csv file for the given experiment specification and returns a DataFrame containing the relevant columns, filtered to only include rows matching the output label specified in the experiment specification. The SMILES strings are canonicalized, and any rows with invalid SMILES strings are dropped. """
    comparison_path = spec.output_dir / "green_vs_quickvina.csv"
    if not comparison_path.exists():
        raise FileNotFoundError(comparison_path)
    columns = [
        "smiles",
        "label",
        "seed",
        "green_score",
        "quickvina2_gpu_raw_values",
        "uses_block",
        "has_protection",
    ]
    df = pd.read_csv(comparison_path, usecols=lambda col: col in columns)
    df = df.loc[df["label"] == spec.output_label].copy()
    # Green outputs are usually already canonical, but canonicalizing here makes
    # the merge robust to different SMILES forms in oracle_history.csv.
    df["canonical_smiles"] = df["smiles"].map(canonicalize_smiles)
    df = df.dropna(subset=["canonical_smiles"])
    return df


def load_oracle_seed(oracle_root: Path, oracle_name: str, seed: int) -> tuple[int, pd.DataFrame] | None:
    """ 
    Loads the oracle history for a given seed from the specified oracle experiment and returns a tuple containing the total number of molecules in the oracle history and a DataFrame containing the relevant columns
    for rows with non-zero syntheseus reward, with canonicalized SMILES strings. If the oracle history file is missing or cannot be loaded, returns None. 
    """
    oracle_path = oracle_root / oracle_name / f"seed{seed}" / "oracle_history.csv"
    if not oracle_path.exists():
        print(f"[warning] Missing oracle history: {oracle_path}")
        return None
    columns = [
        "smiles",
        "syntheseus_reward",
        "quickvina2_gpu_raw_values",
    ]
    df = pd.read_csv(oracle_path, usecols=lambda col: col in columns)
    total_molecules = len(df)
    # Only syntheseus-positive rows can become solved, so only these need RDKit canonicalization.
    df = df.loc[df["syntheseus_reward"] != 0].copy()
    df["seed"] = seed
    df["canonical_smiles"] = df["smiles"].map(canonicalize_smiles)
    return total_molecules, df.dropna(subset=["canonical_smiles"])


def apply_enforced_filters(merged: pd.DataFrame, experiment_name: str) -> pd.DataFrame:
    """
    Applies additional filters to the merged DataFrame based on the experiment name. For experiments that include "blocks" in their name, only rows where uses_block is True are included. 
    For experiments that include "no_deprotection" in their name, only rows where has_protection is False are included. 
    This ensures that the summary statistics for each experiment are calculated based on the appropriate subset of molecules that meet the expected conditions for that experiment.
    """
    filtered = merged.copy()
    lower_name = experiment_name.lower()
    if "blocks" in lower_name:
        filtered = filtered.loc[filtered["uses_block"] == True].copy()
    if "no_deprotection" in lower_name:
        filtered = filtered.loc[filtered["has_protection"] == False].copy()
    return filtered


def summarize_solved_pool(df: pd.DataFrame) -> dict[str, object]:
    """
    Given a DataFrame containing information about solved molecules, calculates summary statistics including the number of solved molecules, mean and standard deviation of docking scores, 
    mean and standard deviation of QED scores, mean and standard deviation of green scores, and the number of unique Bemis-Murcko scaffolds. Returns a dictionary containing these summary statistics.
    """
    solved_df = df.copy()

    # Only calculate QED and BMS for solved molecules, since these are only needed for the summary statistics of the solved pool. This also avoids unnecessary RDKit calculations for molecules that are not part of the solved pool.
    if not solved_df.empty:
        solved_df["qed"] = solved_df["smiles"].map(qed_from_smiles)
        solved_df["bms"] = solved_df["smiles"].map(bms_from_smiles)

    return {
        "pooled_solved_molecules": len(solved_df),
        "docking_mean": solved_df["quickvina2_gpu_raw_values"].mean() if not solved_df.empty else math.nan,
        "docking_std": solved_df["quickvina2_gpu_raw_values"].std() if len(solved_df) > 1 else math.nan,
        "qed_mean": solved_df["qed"].mean() if not solved_df.empty else math.nan,
        "qed_std": solved_df["qed"].std() if len(solved_df) > 1 else math.nan,
        "green_score_mean": solved_df["green_score"].mean() if not solved_df.empty else math.nan,
        "green_score_std": solved_df["green_score"].std() if len(solved_df) > 1 else math.nan,
        "bms_unique": solved_df["bms"].dropna().nunique() if not solved_df.empty else 0,
    }


def summarize_experiment(spec: ExperimentSpec, oracle_root: Path, seeds: Iterable[int], expected_seed_count: int) -> list[dict[str, object]]:
    """
    Summarizes the results for a given experiment specification by loading the green_vs_quickvina.csv data for the experiment and the oracle history for each seed, merging the data, 
    applying enforced filters based on the experiment name, and calculating summary statistics for the solved molecules. 
    Returns a list of dictionaries containing the summary statistics for each seed, as well as pooled summary statistics for all solved molecules across seeds.
    """
    green_df = load_green_data(spec)
    solved_frames: list[pd.DataFrame] = []
    solved_counts: list[int] = []
    non_solved_counts: list[int] = []
    available_seeds: list[int] = []

    merge_columns = [
        "seed",
        "canonical_smiles",
        "green_score",
        "uses_block",
        "has_protection",
    ]

    # If the green_vs_quickvina.csv file contains quickvina2_gpu_raw_values, include this in the merge columns so that it can be used in the summary statistics. If it is not present, the summary statistics will be calculated based on the oracle history data alone.
    if "quickvina2_gpu_raw_values" in green_df.columns:
        merge_columns.append("quickvina2_gpu_raw_values")
    green_merge = green_df[merge_columns].drop_duplicates(["seed", "canonical_smiles"])

    # We iterate over the seeds and attempt to load the oracle history for each seed. If the oracle history is missing for a seed, we skip that seed and do not include it in the summary statistics. 
    # For seeds with available oracle history, we merge the oracle data with the green data, apply enforced filters based on the experiment name, and calculate summary statistics for the solved molecules. We also keep track of which seeds had available oracle history so that we can report this in the final summary.
    for seed in seeds:
        oracle_data = load_oracle_seed(oracle_root, spec.oracle_name, seed)
        if oracle_data is None:
            continue
        total_molecules, oracle_df = oracle_data
        available_seeds.append(seed)
        # Merge on canonical SMILES rather than raw strings because the oracle
        # and green outputs can represent the same molecule differently.
        merged = oracle_df.merge(
            green_merge,
            on=["seed", "canonical_smiles"],
            how="left",
            suffixes=("", "_green"),
        )
        if "quickvina2_gpu_raw_values_green" in merged.columns:
            # Prefer the oracle docking value, but fill gaps from the green
            # comparison file when the column was present there.
            merged["quickvina2_gpu_raw_values"] = merged["quickvina2_gpu_raw_values"].fillna(
                merged["quickvina2_gpu_raw_values_green"]
            )
        solved = apply_enforced_filters(merged.dropna(subset=["green_score"]), spec.output_label)
        solved_frames.append(solved)
        solved_count = len(solved)
        solved_counts.append(solved_count)
        non_solved_counts.append(total_molecules - solved_count)

    # After processing all seeds, we concatenate the solved DataFrames for each seed into a single DataFrame representing the pooled solved molecules across all seeds.
    # We then calculate summary statistics for the pooled solved molecules and include these in the final summary dictionary that is returned for this experiment specification.
    pooled_solved = pd.concat(solved_frames, ignore_index=True) if solved_frames else pd.DataFrame()
    successful_replicates = sum(count > 0 for count in solved_counts)
    common = {
        "experiment": spec.table_label,
        "source_output_label": spec.output_label,
        "source_oracle_experiment": spec.oracle_name,
        "expected_replicates": expected_seed_count,
        "available_replicates": len(available_seeds),
        "successful_replicates_N": successful_replicates,
        "available_seeds": ";".join(str(seed) for seed in available_seeds),
        "solved_mean": mean_value(solved_counts),
        "solved_std": sample_std(solved_counts),
        "non_solved_mean": mean_value(non_solved_counts),
        "non_solved_std": sample_std(non_solved_counts),
    }

    row = common.copy()
    row.update(summarize_solved_pool(pooled_solved))
    return [row]



def build_table_summary(
    outputs_root: Path,
    oracle_root: Path,
    baseline_experiment: str = DEFAULT_BASELINE_EXPERIMENT,
    seeds: Iterable[int] = DEFAULT_SEEDS,
) -> pd.DataFrame:
    
    """
    Builds a summary DataFrame for all experiments discovered under the outputs root directory, using the oracle history data from the specified oracle root directory. 
    The baseline experiment is included as a special case to ensure it is always present and listed first in the table summary. The summary includes statistics for each experiment based on the available seeds, as well as pooled statistics for all solved molecules across seeds.
    """

    resolved_outputs_root = outputs_root.expanduser().resolve()
    resolved_oracle_root = oracle_root.expanduser().resolve()
    seed_tuple = tuple(seeds)

    rows: list[dict[str, object]] = []
    for spec in discover_experiments(resolved_outputs_root, baseline_experiment):
        rows.extend(
            summarize_experiment(
                spec,
                resolved_oracle_root,
                seed_tuple,
                expected_seed_count=len(seed_tuple),
            )
        )
    return pd.DataFrame(rows)


def write_table_summary(
    outputs_root: Path,
    oracle_root: Path,
    output_path: Path,
    baseline_experiment: str = DEFAULT_BASELINE_EXPERIMENT,
    seeds: Iterable[int] = DEFAULT_SEEDS,
) -> Path:
    """
    Builds a summary DataFrame for all experiments discovered under the outputs root directory, using the oracle history data from the specified oracle root directory.
    """
    summary_df = build_table_summary(
        outputs_root=outputs_root,
        oracle_root=oracle_root,
        baseline_experiment=baseline_experiment,
        seeds=seeds,
    )
    resolved_output = output_path.expanduser().resolve()
    resolved_output.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(resolved_output, index=False)
    return resolved_output


def main() -> None:
    """CLI entry point for writing the LaTeX-table summary CSV."""
    parser = argparse.ArgumentParser(
        description="Build CSV inputs for the 9KQ3 LaTeX docking/QED/green-score table."
    )
    parser.add_argument("--outputs-root", type=Path, default=DEFAULT_OUTPUTS_ROOT)
    parser.add_argument("--oracle-root", type=Path, default=DEFAULT_ORACLE_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--baseline-experiment", default=DEFAULT_BASELINE_EXPERIMENT)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(DEFAULT_SEEDS))
    args = parser.parse_args()

    output_path = write_table_summary(
        outputs_root=args.outputs_root,
        oracle_root=args.oracle_root,
        output_path=args.output,
        baseline_experiment=args.baseline_experiment,
        seeds=tuple(args.seeds),
    )
    print(f"Wrote table summary CSV to: {output_path}")


if __name__ == "__main__":
    main()
