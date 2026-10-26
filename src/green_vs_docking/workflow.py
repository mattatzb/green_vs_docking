"""Score one comparison and write its existing CSV and KDE artifacts."""

from __future__ import annotations

from pathlib import Path
from rdkit import RDLogger
from .analysis.analyzer import analyze_experiments
from .analysis.scoring import summarize_scores
from .analysis.comparisons import build_quickvina_comparison_df
from .analysis.biomass_utilization import build_block_mass_lookup
from .analysis.utils import load_block_smiles
from .plots.molecule import HAS_PLOTTING, save_green_vs_quickvina_plot
from .paths import DEFAULT_SOLVENT_GUIDE_PATH, DEFAULT_BLOCK_FILE_PATH, default_namerxn_binary


def run_analysis(
    experiments: dict[str, Path], output_dir: Path,
    solvent_guide: Path = DEFAULT_SOLVENT_GUIDE_PATH,
    block_file: Path = DEFAULT_BLOCK_FILE_PATH,
    namerxn_binary: Path | None = None,
) -> None:
    RDLogger.DisableLog("rdApp.error")
    if namerxn_binary is None:
        namerxn_binary = default_namerxn_binary()
    block_smiles = load_block_smiles(block_file)
    block_mass_lookup = build_block_mass_lookup(block_smiles)
    namerxn_binary = namerxn_binary.expanduser().resolve() if namerxn_binary else None
    scores_df, summary_df = analyze_experiments(
        experiments=experiments,
        reagent_guide_path=solvent_guide.expanduser().resolve(),
        block_smiles=block_smiles,
        block_mass_lookup=block_mass_lookup,
        namerxn_binary=namerxn_binary,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    scores_path = output_dir / "molecule_scores.csv"
    summary_path = output_dir / "summary_scores.csv"
    comparison_path = output_dir / "green_vs_quickvina.csv"
    plot_path = output_dir / "green_vs_quickvina_kde.png"
    scores_df.to_csv(scores_path, index=False)
    summary_df.to_csv(summary_path, index=False)

    block_summary_path = output_dir / "summary_scores_blocks_only.csv"
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
    block_comparison_path = output_dir / "green_vs_quickvina_blocks.csv"
    block_plot_path = output_dir / "green_vs_quickvina_blocks_kde.png"
    saved_block_plot_path = None
    protection_free_comparison_path = output_dir / "green_vs_quickvina_no_protection.csv"
    protection_free_plot_path = output_dir / "green_vs_quickvina_no_protection_kde.png"
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

