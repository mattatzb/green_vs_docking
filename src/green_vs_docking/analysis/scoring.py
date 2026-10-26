"""Per-route scoring, filtering flags, and summary-statistics utilities."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .atom_economy import route_atom_economy
from .solvent_score import score_agents
from .temperature import temperature_score_from_range
from .utils import canonicalize_smiles, load_allowed_smiles, load_json
from .biomass_utilization import compute_biomass_utilization


DEFAULT_WEIGHTS = {
    "atom_economy": 0.25,
    "solvent_score": 0.25,
    "temperature_score": 0.25,
    "biomass_utilization": 0.25,
}

PROTECTION_KEYWORD = "protection"


@dataclass(frozen=True)
class RunSpec:  
    """Paths and labels for one experiment seed run."""

    label: str
    seed: int
    tracker_path: Path
    oracle_history_path: Path


def reaction_nodes(route_entry: dict[str, Any]) -> list[dict[str, Any]]:
    """ Returns a list of reaction nodes from the route entry, sorted by their depth in the route. Reaction nodes are identified as nodes that have the "is_rxn" key set to True. """
    nodes = [value for key, value in route_entry.items() if key.startswith("node_") and isinstance(value, dict)]
    reactions = [node for node in nodes if node.get("is_rxn")]
    return sorted(reactions, key=lambda node: node.get("depth", 9999))

def molecule_nodes(route_entry: dict[str, Any]) -> list[dict[str, Any]]:
    """ Returns a list of molecule nodes from the route entry, sorted by their depth in the route. Molecule nodes are identified as nodes that have the "is_mol" key set to True. """
    nodes = [value for key, value in route_entry.items() if key.startswith("node_") and isinstance(value, dict)]
    molecules = [node for node in nodes if node.get("is_mol")]
    return sorted(molecules, key=lambda node: node.get("depth", 9999))

def route_uses_block(route_entry: dict[str, Any], block_smiles: set[str]) -> bool:
    """ Checks if any molecule node in the route entry has a SMILES string that is present in the block_smiles set. Returns True if at least one molecule node uses a block, and False otherwise. """
    for node in molecule_nodes(route_entry):
        smiles = node.get("mol_smiles")
        canonical_mol = canonicalize_smiles(smiles) if smiles else None
        if canonical_mol is not None and canonical_mol in block_smiles:
            return True

    return False


def route_has_protection(route_entry: dict[str, Any]) -> bool:
    """ Checks if any reaction node in the route entry has a reaction class or name that contains the PROTECTION_KEYWORD. Returns True if at least one reaction node indicates the use of protection, and False otherwise. """
    keyword = PROTECTION_KEYWORD
    for node in reaction_nodes(route_entry):
        for field in ("rxn_class", "rxn_name"):
            value = node.get(field)
            if isinstance(value, str) and keyword in value.lower():
                return True
    return False

def route_component_scores(
    route_entry: dict[str, Any],
    reagent_score_map: dict[str, float],
    default_reagent_score: float,
) -> dict[str, float]:
    """ 
    Computes the component scores for a given route entry, including the number of steps, atom economy, solvent score, and temperature score. 
    The function extracts the reaction nodes from the route entry, calculates the atom economy using the route_atom_economy function, computes the solvent score 
    by scoring the agents in each reaction node using the provided reagent score map and default reagent score, and calculates the temperature score by scoring the temperature ranges in each reaction node using the temperature_score_from_range function. 
    The function returns a dictionary containing all the computed component scores. If there are no reaction nodes in the route entry, it returns default values for each component score.
    """

    # We first extract the reaction nodes from the route entry using the reaction_nodes function. 
    # If there are no reaction nodes, we return a dictionary with default values for each component score, including 0 for the number of steps and NaN for the atom economy, solvent score, and temperature score.
    rxn_nodes = reaction_nodes(route_entry)
    if not rxn_nodes:
        return {
            "num_steps": 0,
            "atom_economy": np.nan,
            "solvent_score": np.nan,
            "temperature_score": np.nan,
        }

    # We calculate the atom economy for the route using the route_atom_economy function, which takes the entire route entry as input.
    atom_economy = route_atom_economy(route_entry)

    # We compute the solvent score by iterating over each reaction node, scoring the agents in each node using the score_agents function with the provided reagent score map and default reagent score, and taking the minimum solvent score found among all reaction nodes as the overall solvent score for the route. 
    # If there are no solvent scores computed (e.g., if no reaction nodes have agents), we use the default reagent score as the solvent score.
    solvent_scores = [
        score_agents(node.get("agents") or [], reagent_score_map, default_reagent_score)
        for node in rxn_nodes
    ]
    solvent_score = float(min(solvent_scores)) if solvent_scores else default_reagent_score

    # We calculate the temperature score by iterating over each reaction node, scoring the temperature range in each node using the temperature_score_from_range function, and taking the mean of all valid temperature scores found among the reaction nodes as the overall temperature score for the route. 
    # If there are no valid temperature scores computed (e.g., if no reaction nodes have valid temperature ranges), we set the temperature score to NaN.
    temperature_scores = [temperature_score_from_range(node.get("temperature")) for node in rxn_nodes]
    valid_temperature_scores = [score for score in temperature_scores if score is not None]
    temperature_score = float(np.mean(valid_temperature_scores)) if valid_temperature_scores else np.nan

    return {
        "num_steps": len(rxn_nodes),
        "atom_economy": atom_economy,
        "solvent_score": solvent_score,
        "temperature_score": temperature_score,
    }


def aggregate_green_score(component_scores: dict[str, float], weights: dict[str, float]) -> float:
    """Aggregate defined component scores with weight renormalization."""

    weighted_values: list[float] = []
    valid_weights: list[float] = []

    # We iterate over each component and its corresponding weight in the weights dictionary. For each component, we retrieve its score from the component_scores dictionary. 
    # If the score is valid (i.e., not None and not NaN), we multiply it by its weight and add it to the list of weighted values, and we also add the weight to the list of valid weights.
    for key, weight in weights.items():
        value = component_scores.get(key, np.nan)
        if value is None or (isinstance(value, float) and math.isnan(value)):
            continue
        weighted_values.append(float(value) * float(weight))
        valid_weights.append(float(weight))

    # If there are no valid weighted values or if the sum of valid weights is zero, we return NaN as the overall green score. Otherwise, we calculate the green score as the sum of weighted values divided by the sum of valid weights and return it as a float.
    if not weighted_values or sum(valid_weights) == 0:
        return np.nan
    return float(sum(weighted_values) / sum(valid_weights))


def score_tracker_file(
    tracker_path: Path,
    oracle_history_path: Path,
    reagent_score_map: dict[str, float],
    default_reagent_score: float,
    weights: dict[str, float],
    label: str,
    seed: int,
    block_smiles: set[str] | None = None,
    block_mass_lookup: dict[str, float] | None = None,
    namerxn_binary: Path | None = None,
) -> pd.DataFrame:
    """
    Scores the routes in a given tracker file by computing the component scores (number of steps, atom economy, solvent score, temperature score, and biomass utilization) for each route and aggregating them into an overall green score using the provided weights.
    """
    tracker = load_json(tracker_path)
    allowed_smiles = load_allowed_smiles(oracle_history_path)

    rows: list[dict[str, Any]] = []

    # We iterate over each route entry in the tracker, where each entry is identified by a SMILES string and contains information about the route. 
    # For each route entry, we first canonicalize the SMILES string and check if it is present in the allowed_smiles set (which contains the canonical SMILES strings of molecules with non-zero syntheseus reward in the oracle history). 
    # If the canonical SMILES string is not valid or not allowed, we skip scoring that route entry.
    for smiles, route_entry in tracker.items():
        canonical_smiles = canonicalize_smiles(smiles)
        if canonical_smiles is None or canonical_smiles not in allowed_smiles:
            continue
        
        # We compute the component scores for the route entry using the route_component_scores function, which calculates the number of steps, atom economy, solvent score, and temperature score based on the reaction nodes in the route.
        components = route_component_scores(route_entry, reagent_score_map, default_reagent_score)

        # We check if the route uses any blocks by calling the route_uses_block function with the route entry and the block_smiles set. 
        # If the route uses a block, we compute the biomass utilization score using the compute_biomass_utilization function, which requires the route entry, block_smiles set, block_mass_lookup dictionary, and namerxn_binary path. 
        # If the route does not use any blocks, we set the biomass utilization score to 0.0.
        use_block = route_uses_block(route_entry, block_smiles) if block_smiles else False
        if use_block:
            if block_mass_lookup is None or namerxn_binary is None:
                raise ValueError("Biomass utilization requires block masses and a NameRxn binary.")
            biomass_result = compute_biomass_utilization(
                route_entry=route_entry,
                block_smiles=block_smiles,
                block_mass_lookup=block_mass_lookup,
                namerxn_binary=namerxn_binary,
            )
            biomass_util = float(biomass_result.utilization)
        else:
            biomass_util = 0.0
        components["biomass_utilization"] = biomass_util

        # We check if the route has any protection steps by calling the route_has_protection function with the route entry. 
        # This function checks if any reaction node in the route indicates the use of protection based on its reaction class or name. The presence of protection steps is recorded as a boolean value.
        has_protection = route_has_protection(route_entry)

        # We aggregate the component scores into an overall green score for the route using the aggregate_green_score function, which takes the component scores and the weights as input and calculates a weighted average of the valid component scores to produce the final green score.
        green_score = aggregate_green_score(components, weights)

        rows.append(
            {
                "smiles": canonical_smiles,
                "label": label,
                "seed": seed,
                "num_steps": components["num_steps"],
                "atom_economy": components["atom_economy"],
                "solvent_score": components["solvent_score"],
                "temperature_score": components["temperature_score"],
                "green_score": green_score,
                "biomass_utilization": biomass_util,
                "oracle_calls": route_entry.get("oracle_calls"),
                "rxn_steps_from_tracker": route_entry.get("rxn_steps"),
                "uses_block": use_block,
                "has_protection": has_protection,
            }
        )

    return pd.DataFrame(rows)


def summarize_scores(scores_df: pd.DataFrame) -> pd.DataFrame:
    """ Summarizes the scores by computing the mean and standard deviation of each metric for each label and seed, as well as the average metrics across seeds for each label. The function returns a summary DataFrame with the computed statistics. """

    metric_columns = ["green_score", "atom_economy", "solvent_score", "temperature_score", "num_steps", "biomass_utilization"]

    # We first compute the mean and standard deviation of each metric for each combination of label and seed by grouping the scores_df DataFrame by "label" and "seed" and applying the mean and std aggregation functions to the specified metric columns. 
    mean_df = (
        scores_df.groupby(["label", "seed"], as_index=False)[metric_columns]
        .mean()
    )
    std_df = (
        scores_df.groupby(["label", "seed"], as_index=False)[metric_columns]
        .std()
        .rename(columns={col: f"{col}_std" for col in metric_columns})
    )

    # We then merge the mean and standard deviation DataFrames on the "label" and "seed" columns to create a per-seed summary DataFrame. 
    # We also compute the count of molecules for each label and seed by grouping the scores_df DataFrame and merging the counts into the per-seed summary. 
    per_seed = mean_df.merge(std_df, on=["label", "seed"], how="left")
    counts = (
        scores_df.groupby(["label", "seed"], as_index=False)
        .size()
        .rename(columns={"size": "num_molecules"})
    )

    per_seed = per_seed.merge(counts, on=["label", "seed"], how="left")

    # If the "uses_block" column is present in the scores_df DataFrame, we calculate the average block usage percentage for each label and seed and merge that into the per-seed summary as well.
    if "uses_block" in scores_df.columns:
        block_usage = (
            scores_df.groupby(["label", "seed"], as_index=False)["uses_block"]
            .mean()
            .rename(columns={"uses_block": "block_usage_pct"})
        )
        block_usage["block_usage_pct"] = block_usage["block_usage_pct"] * 100
        per_seed = per_seed.merge(block_usage, on=["label", "seed"], how="left")
    else:
        per_seed["block_usage_pct"] = np.nan

    per_seed["summary_level"] = "seed"

    # We compute the average metrics across seeds for each label by grouping the per-seed summary DataFrame by "label" and applying the mean aggregation function to the metric columns, as well as the count of molecules and block usage percentage. 
    # We also compute the standard deviation of the metrics across seeds for each label by applying the std aggregation function to the metric columns. 
    experiment_metric_means = (
        per_seed.groupby("label", as_index=False)[metric_columns]
        .mean()
    )

    experiment_metric_stds = (
        per_seed.groupby("label", as_index=False)[metric_columns]
        .std()
        .rename(columns={col: f"{col}_std" for col in metric_columns})
    )
    experiment_other_means = (
        per_seed.groupby("label", as_index=False)[["num_molecules", "block_usage_pct"]]
        .mean()
    )

    # We then merge the average metrics and their standard deviations into an experiment-level summary DataFrame, which contains the average metrics across seeds for each label, as well as the number of molecules and block usage percentage. 
    experiment_avg = (
        experiment_metric_means
        .merge(experiment_metric_stds, on="label", how="left")
        .merge(experiment_other_means, on="label", how="left")
    )

    # We set the "seed" column to "average" and the "summary_level" column to "experiment_average" for the experiment-level summary.
    experiment_avg["seed"] = "average"
    experiment_avg["summary_level"] = "experiment_average"

    summary_df = pd.concat([per_seed, experiment_avg], ignore_index=True, sort=False)
    final_columns = [
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
    return summary_df[final_columns]
