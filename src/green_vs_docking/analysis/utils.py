"""Shared file, seed-directory, and SMILES utilities for green_vs_docking."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem


SEED_DIR_PATTERN = re.compile(r"seed(\d+)$")


def canonicalize_smiles(smiles: str | None) -> str | None:
    """Return RDKit canonical SMILES, or None for empty/unparseable input."""
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def load_json(path: Path) -> dict[str, Any]:
    """Load a JSON file into a dictionary."""
    with path.open("r") as handle:
        return json.load(handle)


def load_allowed_smiles(oracle_history_path: Path) -> set[str]:
    """Return canonical SMILES with non-zero Syntheseus reward."""
    oracle_df = pd.read_csv(oracle_history_path)
    return {
        canonical
        for smiles in oracle_df.loc[oracle_df["syntheseus_reward"] != 0, "smiles"]
        if isinstance(smiles, str)
        for canonical in [canonicalize_smiles(smiles)]
        if canonical is not None
    }

def load_block_smiles(blocks_path: Path) -> set[str]:
    """ Returns the set of canonical SMILES strings for molecules that are listed in the blocks file. """
    if not blocks_path.exists():
        return set()
    with blocks_path.open("r") as handle:
        # The blocks file is expected to be a CSV file where each line contains one or more SMILES strings separated by commas. We will read the file line by line, split each line by commas, and canonicalize each SMILES string. We will then return the set of unique canonical SMILES strings.
        return {
            canonical
            for line in handle
            for smiles in line.split(",")
            if isinstance(smiles, str)
            for canonical in [canonicalize_smiles(smiles)]
            if canonical is not None
        }

def find_seed_dirs(experiment_dir: Path) -> dict[int, Path]:
    """ Returns a dictionary mapping seed numbers to their corresponding directories within the experiment directory. """
    seed_dirs: dict[int, Path] = {}
    for path in sorted(experiment_dir.glob("seed*")):
        if not path.is_dir():
            continue
        match = SEED_DIR_PATTERN.match(path.name)
        if match is None:
            continue
        seed_dirs[int(match.group(1))] = path
    return seed_dirs
