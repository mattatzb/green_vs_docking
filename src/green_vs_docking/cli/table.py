"""Manuscript table CLI."""

from __future__ import annotations

import argparse
from pathlib import Path
from ..analysis.table_summary import (DEFAULT_OUTPUTS_ROOT, DEFAULT_ORACLE_ROOT, DEFAULT_OUTPUT, DEFAULT_BASELINE_EXPERIMENT, DEFAULT_SEEDS, write_table_summary)

def main(argv: list[str] | None = None) -> None:
    """CLI entry point for writing the LaTeX-table summary CSV."""
    parser = argparse.ArgumentParser(
        description="Build CSV inputs for the 9KQ3 LaTeX docking/QED/green-score table."
    )
    parser.add_argument("--outputs-root", type=Path, default=DEFAULT_OUTPUTS_ROOT)
    parser.add_argument("--oracle-root", type=Path, default=DEFAULT_ORACLE_ROOT)
    parser.add_argument(
        "--output", type=Path, default=None,
        help="CSV destination (default: outputs/9KQ3/latex_table_summary.csv in the repository).",
    )
    parser.add_argument("--baseline-experiment", default=DEFAULT_BASELINE_EXPERIMENT)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(DEFAULT_SEEDS))
    args = parser.parse_args(argv)

    output_path = write_table_summary(
        outputs_root=args.outputs_root,
        oracle_root=args.oracle_root,
        output_path=args.output if args.output is not None else DEFAULT_OUTPUT,
        baseline_experiment=args.baseline_experiment,
        seeds=tuple(args.seeds),
    )
    print(f"Wrote table summary CSV to: {output_path}")


