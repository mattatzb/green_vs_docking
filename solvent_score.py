"""Solvent and reagent scoring from a canonical-SMILES score guide."""

from __future__ import annotations

from pathlib import Path

from .utils import canonicalize_smiles, load_json


def build_reagent_score_map(
    reagent_guide_path: Path,
    manual_scores: dict[str, float] | None = None,
    default_score: float = 0.8,
) -> tuple[dict[str, float], float]:
    """
    Builds a reagent score map based on the provided reagent guide JSON file or manual scores dictionary. The reagent guide JSON file is expected to have a structure where it contains 
    a "categories" key mapping category names to lists of SMILES strings, and a "weights" key mapping category names to their corresponding scores. 
    The manual scores dictionary is expected to map SMILES strings directly to their scores. If the reagent guide file exists, it takes precedence over the manual scores. 
    The function returns a tuple containing the reagent score map (a dictionary mapping canonical SMILES strings to their scores) and the default score to be used for any reagents not found in the map.
    """
    score_map: dict[str, float] = {}

    # If the reagent guide file exists, we load the categories and weights from the file and populate the score map with canonical SMILES strings and their corresponding scores. 
    if reagent_guide_path.exists():
        guide = load_json(reagent_guide_path)
        categories = guide.get("categories", {})
        weights = guide.get("weights", {})

        # We iterate over the categories and their associated SMILES lists, canonicalize each SMILES string, and assign the corresponding score from the weights (or the default score if the category is not found in the weights) to each canonical SMILES string in the score map.
        for category, smiles_list in categories.items():
            score = float(weights.get(category, default_score))
            for smiles in smiles_list:
                canonical = canonicalize_smiles(smiles)
                if canonical is not None:
                    score_map[canonical] = score
    # If the reagent guide file does not exist but manual scores are provided, we populate the score map directly from the manual scores dictionary by canonicalizing each SMILES string and assigning the provided score to each canonical SMILES string in the score map.
    elif manual_scores:
        for smiles, score in manual_scores.items():
            canonical = canonicalize_smiles(smiles)
            if canonical is not None:
                score_map[canonical] = float(score)

    return score_map, float(default_score)


def score_agents(
    agents: list[str] | None,
    reagent_score_map: dict[str, float],
    default_reagent_score: float,
) -> float:
    """
    Scores a list of agent SMILES strings based on the provided reagent score map and default reagent score. 
    The function canonicalizes each SMILES string in the agents list, looks up its score in the reagent score map, and returns the minimum score found among the agents. 
    If an agent's canonical SMILES string is not found in the reagent score map, the default reagent score is used for that agent. If the agents list is empty or None, the function returns the default reagent score.
    """
    if not agents:
        return default_reagent_score

    scores: list[float] = []
    for smiles in agents:
        canonical = canonicalize_smiles(smiles)
        if canonical is None:
            continue
        scores.append(reagent_score_map.get(canonical, default_reagent_score))
    return float(min(scores)) if scores else default_reagent_score
