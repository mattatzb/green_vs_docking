"""Command-line entry point for route scoring and green-vs-docking plots."""

from __future__ import annotations

import argparse
from pathlib import Path

from rdkit import RDLogger

from .analyzer import analyze_experiments
from .scoring import summarize_scores
from .plots import HAS_PLOTTING, build_quickvina_comparison_df, save_green_vs_quickvina_plot
from .biomass_utilization import build_block_mass_lookup
from .utils import load_block_smiles

DEFAULT_SOLVENT_GUIDE_PATH = Path(__file__).resolve().parent / "data" / "solvent_guide.json"
DEFAULT_BLOCK_FILE_PATH = Path(__file__).resolve().parent / "data" / "blocks.smi"

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
        default=Path("/work/liac/tatzber/HazELNut/namerxn"),
        help="Path to the NameRxn binary used for biomass utilization.",
    )
    return parser


def main() -> None:
    """Run scoring, write CSV summaries, and save available comparison plots."""
    RDLogger.DisableLog("rdApp.error")
    parser = build_parser()
    args = parser.parse_args()

    block_smiles = load_block_smiles(args.block_file)
    block_mass_lookup = build_block_mass_lookup(block_smiles)
    namerxn_binary = args.namerxn_binary.expanduser().resolve() if args.namerxn_binary else None
    experiments: dict[str, Path] = {}
    for label, path in args.experiment:
        normalized_label = path.name if label == "conditions_no_block" else label
        experiments[normalized_label] = path
    scores_df, summary_df = analyze_experiments(
        experiments=experiments,
        reagent_guide_path=args.solvent_guide.expanduser().resolve(),
        block_smiles=block_smiles,
        block_mass_lookup=block_mass_lookup,
        namerxn_binary=namerxn_binary,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    scores_path = args.output_dir / "molecule_scores.csv"
    summary_path = args.output_dir / "summary_scores.csv"
    comparison_path = args.output_dir / "green_vs_quickvina.csv"
    plot_path = args.output_dir / "green_vs_quickvina_kde.png"
    scores_df.to_csv(scores_path, index=False)
    summary_df.to_csv(summary_path, index=False)

    block_summary_path = args.output_dir / "summary_scores_blocks_only.csv"
    block_summary_df = None
    if "uses_block" in scores_df.columns:
        block_labels = {
            label
            for label in scores_df["label"].unique()
            if isinstance(label, str) and "blocks" in label.lower()
        }
        if block_labels:
            block_filtered_scores = scores_df.copy()
            uses_block_col = block_filtered_scores["uses_block"].fillna(False)
            block_mask = block_filtered_scores["label"].isin(block_labels)
            block_filtered_scores = block_filtered_scores.loc[
                (~block_mask) | (uses_block_col == True)
            ]
            block_summary_df = summarize_scores(block_filtered_scores)
            block_summary_df.to_csv(block_summary_path, index=False)

    default_df = build_quickvina_comparison_df(
        scores_df=scores_df,
        experiments=experiments,
    )
    block_df = build_quickvina_comparison_df(
        scores_df=scores_df,
        experiments=experiments,
        require_block=True,
    )
    protection_free_df = build_quickvina_comparison_df(
        scores_df=scores_df,
        experiments=experiments,
        exclude_protection=True,
    )

    default_df.to_csv(comparison_path, index=False)
    saved_default_plot_path = save_green_vs_quickvina_plot(comparison_df=default_df, output_path=plot_path)
    block_comparison_path = args.output_dir / "green_vs_quickvina_blocks.csv"
    block_plot_path = args.output_dir / "green_vs_quickvina_blocks_kde.png"
    saved_block_plot_path = None
    protection_free_comparison_path = args.output_dir / "green_vs_quickvina_no_protection.csv"
    protection_free_plot_path = args.output_dir / "green_vs_quickvina_no_protection_kde.png"
    saved_protection_free_plot_path = None

    if not block_df.empty:
        block_df.to_csv(block_comparison_path, index=False)
        saved_block_plot_path = save_green_vs_quickvina_plot(block_df, block_plot_path)

    if not protection_free_df.empty:
        protection_free_df.to_csv(protection_free_comparison_path, index=False)
        saved_protection_free_plot_path = save_green_vs_quickvina_plot(protection_free_df, protection_free_plot_path)

    print("Summary per seed and experiment average:")
    if summary_df.empty:
        print("No valid experiment data found.")
    else:
        print(summary_df.to_string(index=False))

    print(f"\nWrote molecule-level scores to: {scores_path}")
    print(f"Wrote summary scores to: {summary_path}")
    print(f"Wrote green/docking comparison data to: {comparison_path}")

    if block_summary_df is not None:
        print(f"Wrote block-only summary scores to: {block_summary_path}")

    if saved_default_plot_path is not None:
        print(f"Wrote default green/docking KDE plot to: {saved_default_plot_path}")
    
    if saved_block_plot_path is not None:
        print(f"Wrote block-only green/docking KDE plot to: {saved_block_plot_path}")
    if saved_protection_free_plot_path is not None:
        print(f"Wrote protection-free green/docking KDE plot to: {saved_protection_free_plot_path}")

    if saved_default_plot_path is None and saved_block_plot_path is None and saved_protection_free_plot_path is None:
        if HAS_PLOTTING:
            print("Skipped … no comparison rows …")
        else:
            print("Skipped … matplotlib/seaborn not available.")


if __name__ == "__main__":
    main()
