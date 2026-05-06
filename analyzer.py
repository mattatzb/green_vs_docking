"""Experiment-level orchestration for green metric scoring."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .scoring import DEFAULT_WEIGHTS, score_tracker_file, summarize_scores
from .solvent_score import build_reagent_score_map
from .utils import find_seed_dirs


def analyze_experiments(
    experiments: dict[str, Path],
    reagent_guide_path: Path,
    weights: dict[str, float] | None = None,
    default_reagent_score: float = 0.8,
    block_smiles: set[str] | None = None,
    block_mass_lookup: dict[str, float] | None = None,
    namerxn_binary: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Score every valid seed from each experiment and return molecule-level and
    summary-level tables.

    Each experiment directory is expected to contain `seed*` subdirectories with
    both `oracle_history.csv` and `syntheseus_results/smiles_rxn_tracker.json`.
    Missing seeds are skipped, while malformed tracker JSON files produce a
    warning and do not stop the full analysis.
    """
    active_weights = weights or DEFAULT_WEIGHTS
    reagent_score_map, reagent_default = build_reagent_score_map(
        reagent_guide_path=reagent_guide_path,
        default_score=default_reagent_score,
    )

    all_scores: list[pd.DataFrame] = []
    for label, experiment_dir in experiments.items():
        for seed, seed_dir in find_seed_dirs(experiment_dir).items():
            tracker_path = seed_dir / "syntheseus_results" / "smiles_rxn_tracker.json"
            oracle_history_path = seed_dir / "oracle_history.csv"
            if not tracker_path.exists() or not oracle_history_path.exists():
                continue

            try:
                all_scores.append(
                    score_tracker_file(
                        tracker_path=tracker_path,
                        oracle_history_path=oracle_history_path,
                        reagent_score_map=reagent_score_map,
                        default_reagent_score=reagent_default,
                        weights=active_weights,
                        label=label,
                        seed=seed,
                        block_smiles=block_smiles,
                        block_mass_lookup=block_mass_lookup,
                        namerxn_binary=namerxn_binary,
                    )
                )
            except json.JSONDecodeError as exc:
                print(
                    f"[warning] Skipping experiment '{label}' seed {seed} because {tracker_path} contains invalid JSON ({exc})."
                )
                continue

    if not all_scores:
        empty_scores = pd.DataFrame(
            columns=[
                "smiles",
                "label",
                "seed",
                "num_steps",
                "atom_economy",
                "solvent_score",
                "temperature_score",
                "green_score",
                "biomass_utilization",
                "oracle_calls",
                "rxn_steps_from_tracker",
            ]
        )
        empty_summary = pd.DataFrame(
            columns=[
                "label",
                "seed",
                "summary_level",
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
                "block_usage_pct",
            ]
        )
        return empty_scores, empty_summary

    scores_df = pd.concat(all_scores, ignore_index=True)
    summary_df = summarize_scores(scores_df)
    return scores_df, summary_df
