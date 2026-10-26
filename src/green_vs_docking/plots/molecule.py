"""Molecule-level KDE rendering."""

from __future__ import annotations

from pathlib import Path
import pandas as pd
from ..analysis.comparisons import SPECIAL_SUFFIXES

try:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter
    import seaborn as sns

    HAS_PLOTTING = True
except Exception:
    HAS_PLOTTING = False


DEFAULT_PALETTE = {
    "baseline": "#6ea86a",
    "conditions_no_block": "#c1b46a",
}


def _resolve_palette_color(resolved_palette: dict[str, str], experiment_name: str) -> str:
    """ 
    Resolves the color for a given experiment name based on the provided palette and special suffixes. The function first checks if there is a direct match for the experiment name in the resolved palette. 
    If not, it checks if the experiment name ends with any of the defined special suffixes (e.g., " (blocks only)", " (no protection)") and attempts to find a base label by removing the suffix. If a base label is found in the palette, its color is used. 
    If no match is found, a default color is returned.
    """
    color = resolved_palette.get(experiment_name)
    if color:
        return color
    for suffix in SPECIAL_SUFFIXES:
        if experiment_name.endswith(suffix):
            base_label = experiment_name[: -len(suffix)]
            if base_label:
                color = resolved_palette.get(base_label)
                if color:
                    return color
    return "#4c72b0"


def save_green_vs_quickvina_plot(
    comparison_df: pd.DataFrame,
    output_path: Path,
    palette: dict[str, str] | None = None,
) -> Path | None:
    """
    Save the KDE comparison plot used by the analysis pipeline.

    The docking axis is displayed in the original QuickVina scale, but the
    internal plotting coordinate is inverted so more negative scores appear on
    the right.
    """
    if not HAS_PLOTTING or comparison_df.empty:
        return None

    resolved_palette = dict(DEFAULT_PALETTE)
    if palette:
        resolved_palette.update(palette)

    plot_df = comparison_df.copy()
    experiment_label_column = "plot_label" if "plot_label" in plot_df.columns else "label"
    plot_df["experiment"] = plot_df[experiment_label_column]
    plot_df = plot_df.dropna(subset=["experiment", "quickvina2_gpu_raw_values", "green_score"])
    if plot_df.empty:
        return None

    # Store inverted docking reward so that the most negative original values appear on the right.
    plot_df["_quickvina_inverted"] = -plot_df["quickvina2_gpu_raw_values"].astype(float)

    sns.set_theme(
        style="ticks",
        rc={
            "axes.facecolor": "#f2f2f2",
            "figure.facecolor": "#f2f2f2",
            "axes.edgecolor": "#444444",
            "axes.linewidth": 1.0,
            "grid.color": "#d9d9d9",
            "grid.linestyle": "-",
            "grid.linewidth": 0.6,
            "font.size": 12,
        },
    )

    g = sns.JointGrid(
        data=plot_df,
        x="_quickvina_inverted",
        y="green_score",
        height=7.2,
        ratio=5,
        space=0.08,
    )

    experiment_order = list(dict.fromkeys(plot_df["experiment"].tolist()))
    for experiment_name in experiment_order:
        subset = plot_df.loc[plot_df["experiment"] == experiment_name]
        if subset.empty:
            continue

        color = _resolve_palette_color(resolved_palette, experiment_name)
        sns.kdeplot(
            data=subset,
            x="_quickvina_inverted",
            y="green_score",
            ax=g.ax_joint,
            levels=10,
            thresh=0.05,
            bw_adjust=1.0,
            fill=False,
            linewidths=2,
            color=color,
        )
        sns.kdeplot(
            data=subset,
            x="_quickvina_inverted",
            ax=g.ax_marg_x,
            bw_adjust=1.0,
            fill=False,
            linewidth=2,
            color=color,
        )
        sns.kdeplot(
            data=subset,
            y="green_score",
            ax=g.ax_marg_y,
            bw_adjust=1.0,
            fill=False,
            linewidth=2,
            color=color,
        )

    # Bounds are computed in the original docking scale first because the user
    # sees original QuickVina values even though the plotted coordinate is negated.
    orig_x_values = plot_df["quickvina2_gpu_raw_values"].astype(float)
    orig_x_min = float(orig_x_values.min())
    orig_x_max = float(orig_x_values.max())
    x_upper = 0.0 if orig_x_max <= 0 else orig_x_max
    span = max(x_upper - orig_x_min, 1e-6)
    left_pad = 0.02 * span
    right_pad = 0.0 if x_upper == 0.0 else 0.02 * span
    orig_left_bound = orig_x_min - left_pad
    orig_right_bound = x_upper + right_pad

    # Convert the target original bounds into the inverted axis space.
    inverted_left = -orig_right_bound
    inverted_right = -orig_left_bound

    g.ax_joint.scatter(
        plot_df["_quickvina_inverted"].mean(),
        plot_df["green_score"].mean(),
        s=40,
        color="black",
        zorder=5,
    )
    g.ax_joint.set_xlabel("QuickVina2 GPU Reward")
    g.ax_joint.set_ylabel("Green Score")
    g.ax_joint.set_xlim(inverted_left, inverted_right)
    g.ax_marg_x.set_xlim(inverted_left, inverted_right)

    formatter = FuncFormatter(lambda value, _: f"{-value:g}")
    g.ax_joint.xaxis.set_major_formatter(formatter)
    g.ax_marg_x.xaxis.set_major_formatter(formatter)

    g.ax_joint.grid(False)
    g.ax_marg_x.grid(False)
    g.ax_marg_y.grid(False)

    legend_handles = []
    for name in experiment_order:
        legend_color = _resolve_palette_color(resolved_palette, name)
        legend_handles.append(
            Line2D([0], [0], color=legend_color, lw=2, label=name)
        )
    g.ax_joint.legend(handles=legend_handles, title="Experiment", loc="upper right", frameon=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    g.figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(g.figure)
    return output_path


