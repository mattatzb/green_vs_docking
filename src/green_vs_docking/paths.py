"""Shared checkout paths, independent of the shell's working directory."""

from __future__ import annotations

import os
import shutil
from pathlib import Path


def repository_root() -> Path:
    configured = os.environ.get("GREEN_VS_DOCKING_ROOT")
    if configured:
        return Path(configured).expanduser().resolve()
    for anchor in [Path(__file__).resolve().parent, Path.cwd()]:
        for candidate in [anchor, *anchor.parents]:
            if (candidate / "data/blocks.smi").is_file() and (candidate / "environment.yml").is_file():
                return candidate
    return Path.cwd().resolve()


REPO_ROOT = repository_root()
DEFAULT_DATA_ROOT = REPO_ROOT / "data"
REFERENCE_OUTPUTS_ROOT = REPO_ROOT / "9KQ3_outputs"
DEFAULT_GENERATED_ROOT = REPO_ROOT / "outputs" / "9KQ3"
DEFAULT_SOLVENT_GUIDE_PATH = DEFAULT_DATA_ROOT / "solvent_guide.json"
DEFAULT_BLOCK_FILE_PATH = DEFAULT_DATA_ROOT / "blocks.smi"


def default_namerxn_binary() -> Path:
    return Path(os.environ.get("NAMERXN_BINARY") or shutil.which("namerxn") or "namerxn")
