"""CLI for one baseline/condition scoring comparison."""

from __future__ import annotations

import argparse
from pathlib import Path
from ..paths import DEFAULT_SOLVENT_GUIDE_PATH, DEFAULT_BLOCK_FILE_PATH, default_namerxn_binary
from ..workflow import run_analysis

def parse_experiment(value: str) -> tuple[str, Path]:
    """Parse a CLI experiment argument of the form `label=/path/to/run`."""
    if "=" not in value:
        raise argparse.ArgumentTypeError("Expected NAME=PATH for --experiment.")
    label, raw_path = value.split("=", 1)
    label = label.strip()
    path = Path(raw_path).expanduser().resolve()
    if not label:
        raise argparse.ArgumentTypeError("Experiment label cannot be empty.")
    return label, path


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser used by `run_green_analysis.py`."""
    parser = argparse.ArgumentParser(description="Analyze green metrics across experiment seeds.")
    parser.add_argument(
        "--experiment",
        action="append",
        type=parse_experiment,
        required=True,
        help="Experiment mapping in the form label=/path/to/experiment.",
    )
    parser.add_argument(
        "--solvent-guide",
        type=Path,
        default=DEFAULT_SOLVENT_GUIDE_PATH,
        help="Path to solvent guide JSON.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path.cwd() / "green_metric_outputs",
        help="Directory where CSV outputs will be written.",
    )

    parser.add_argument(
        "--block-file",
        type=Path,
        default=DEFAULT_BLOCK_FILE_PATH,
        help="Path to enforced block SMILES file.",
    )
    parser.add_argument(
        "--namerxn-binary",
        type=Path,
        default=default_namerxn_binary(),
        help="Path to NameRxn (defaults to NAMERXN_BINARY or namerxn on PATH).",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    experiments = {}
    for label, path in args.experiment:
        normalized_label = path.name if label == "conditions_no_block" else label
        experiments[normalized_label] = path
    run_analysis(
        experiments=experiments, output_dir=args.output_dir,
        solvent_guide=args.solvent_guide, block_file=args.block_file,
        namerxn_binary=args.namerxn_binary,
    )
