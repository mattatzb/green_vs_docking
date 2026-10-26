"""Arguments for the complete manuscript pipeline."""

from __future__ import annotations

import argparse
from pathlib import Path
from ..analysis.table_summary import DEFAULT_BASELINE_EXPERIMENT
from ..paths import (DEFAULT_DATA_ROOT, DEFAULT_GENERATED_ROOT, DEFAULT_SOLVENT_GUIDE_PATH,
                     DEFAULT_BLOCK_FILE_PATH, default_namerxn_binary)
from ..pipeline import run_pipeline


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Score the supplied routes and generate all three manuscript artifacts.")
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_GENERATED_ROOT,
                        help="Fresh output directory (default: outputs/9KQ3; references stay separate).")
    parser.add_argument("--from-comparisons", dest="comparison_root", type=Path,
                        help="Reuse comparison CSVs from this directory, skipping route rescoring.")
    parser.add_argument("--baseline-experiment", default=DEFAULT_BASELINE_EXPERIMENT)
    parser.add_argument("--condition", dest="conditions", action="append",
                        help="Condition directory name; repeat as needed (default: all supplied 9KQ3 conditions).")
    parser.add_argument("--solvent-guide", type=Path, default=DEFAULT_SOLVENT_GUIDE_PATH)
    parser.add_argument("--block-file", type=Path, default=DEFAULT_BLOCK_FILE_PATH)
    parser.add_argument("--namerxn-binary", type=Path, default=default_namerxn_binary())
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--aggregate-by-seed", action="store_true")
    modes.add_argument("--resume-from-summaries", action="store_true")
    parser.add_argument("--prefer-blocks-file", dest="prefer_block_file",
                        action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args(argv)
    for name, path in run_pipeline(**vars(args)).items():
        print(f"Wrote {name}: {path}")
