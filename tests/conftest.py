"""Small synthetic fixtures; src imports are configured in pyproject.toml."""
import json
import os
from pathlib import Path
import pytest
os.environ.setdefault("MPLBACKEND", "Agg")


def write_tiny_inputs(root):
    """Two seeds, two solved ester routes, and one unsolved target per seed."""
    import pandas as pd

    for experiment in ["pesticide_9KQ3_basic", "pesticide_9KQ3_solvent"]:
        for seed in range(2):
            directory = root / experiment / f"seed{seed}"
            (directory / "syntheseus_results").mkdir(parents=True)
            tracker = {}
            rows = []
            for index, (target, alcohol) in enumerate([("COC(C)=O", "CO"), ("CCOC(C)=O", "CCO")]):
                tracker[target] = {
                    "oracle_calls": 10 + index, "rxn_steps": 1,
                    "node_1": {"depth": 0, "is_mol": True, "mol_smiles": target},
                    "node_2": {"depth": 1, "is_rxn": True,
                               "rxn_smiles": f"CC(=O)O.{alcohol}>>{target}",
                               "agents": ["O"] if index == 0 else ["CCN(CC)CC"],
                               "temperature": "[20, 30)" if index == 0 else "[80, 90)",
                               "rxn_class": "Esterification" if index == 0 else "Protection"},
                    "node_3": {"depth": 2, "is_mol": True, "mol_smiles": "CC(=O)O"},
                    "node_4": {"depth": 2, "is_mol": True, "mol_smiles": alcohol},
                }
                rows.append({"smiles": target, "syntheseus_reward": 1,
                             "quickvina2_gpu_raw_values": -6.0 - index - seed / 10})
            rows.append({"smiles": "C", "syntheseus_reward": 0, "quickvina2_gpu_raw_values": -3.0})
            pd.DataFrame(rows).to_csv(directory / "oracle_history.csv", index=False)
            (directory / "syntheseus_results/smiles_rxn_tracker.json").write_text(json.dumps(tracker))
    return root


@pytest.fixture
def tiny_data_root(tmp_path):
    return write_tiny_inputs(tmp_path / "data")
