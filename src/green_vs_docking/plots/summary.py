"""Manuscript scatter/error-bar and joint KDE figures."""
from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

# Use plasma palette for all plots by default.
_PLASMA = plt.cm.get_cmap("plasma")
plt.rcParams["image.cmap"] = "plasma"
plt.rcParams["axes.prop_cycle"] = plt.cycler(
    color=_PLASMA(np.linspace(0.1, 0.9, 10))
)
plt.rcParams["xtick.labelsize"] = 18
plt.rcParams["ytick.labelsize"] = 18

def plot_points(df: pd.DataFrame, output_path: Path) -> None:
    """Save a scatter/error-bar plot of mean green score versus docking."""
    if df.empty:
        raise RuntimeError("No experiment data found under the specified directory.")

    df["green_score_std"] = df["green_score_std"].fillna(0.0)
    df["docking_std"] = df["docking_std"].fillna(0.0).abs()

    markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    labels = list(dict.fromkeys(df["label"].astype(str)))
    label_to_marker = {
        label: markers[i % len(markers)] for i, label in enumerate(labels)
    }

    fig, ax = plt.subplots(figsize=(9, 6))
    plotted_labels: set[str] = set()

    for _, row in df.iterrows():
        label = row["label"]
        marker = label_to_marker[label]
        label_text = label if label not in plotted_labels else None
        ax.errorbar(
            row["docking_mean"],
            row["green_score_mean"],
            xerr=row["docking_std"],
            yerr=row["green_score_std"],
            fmt=marker,
            markersize=8,
            capsize=4,
            elinewidth=1.0,
            markeredgewidth=1.0,
            label=label_text,
        )
        plotted_labels.add(label)

    ax.set_xlabel("Docking score")
    ax.set_ylabel("Green chemistry score")
    ax.set_title("Green chemistry-constrained experiments")
    ax.invert_xaxis()
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    ax.legend(fontsize=8, loc="best", frameon=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_joint_kde(
    df: pd.DataFrame,
    output_path: Path,
    docking_min: float = -14.0,
    docking_max: float = 0.0,
) -> None:
    """Save one joint KDE overlay containing all experiment molecule clouds."""
    try:
        import seaborn as sns
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "The joint KDE plot requires seaborn in the active Python environment."
        ) from exc

    required_columns = {"label", "green_score", "quickvina2_gpu_raw_values"}
    if df.empty or not required_columns.issubset(df.columns):
        raise RuntimeError("No molecule-level comparison data found for the joint KDE plot.")

    plot_df = df.dropna(subset=["label", "green_score", "quickvina2_gpu_raw_values"]).copy()
    if plot_df.empty:
        raise RuntimeError("All molecule-level rows have missing green or docking values.")

    # Negating the docking score preserves the desired visual convention:
    # more negative docking values appear on the right, while tick labels show
    # the original docking score.
    plot_df["_quickvina_inverted"] = -plot_df["quickvina2_gpu_raw_values"]
    labels = list(dict.fromkeys(plot_df["label"].astype(str)))
    colors = sns.color_palette("plasma", n_colors=len(labels))

    grid = sns.JointGrid(
        data=plot_df,
        x="_quickvina_inverted",
        y="green_score",
        height=8,
        ratio=5,
        space=0.08,
    )

    legend_handles: list[Line2D] = []
    for label, color in zip(labels, colors):
        label_df = plot_df.loc[plot_df["label"] == label]
        if len(label_df) < 2:
            continue

        sns.kdeplot(
            data=label_df,
            x="_quickvina_inverted",
            y="green_score",
            ax=grid.ax_joint,
            levels=8,
            thresh=0.05,
            fill=False,
            linewidths=0.8,
            color=color,
            warn_singular=False,
        )
        sns.kdeplot(
            data=label_df,
            x="_quickvina_inverted",
            ax=grid.ax_marg_x,
            color=color,
            linewidth=1.0,
            fill=False,
            warn_singular=False,
        )
        sns.kdeplot(
            data=label_df,
            y="green_score",
            ax=grid.ax_marg_y,
            color=color,
            linewidth=1.0,
            fill=False,
            warn_singular=False,
        )
        legend_handles.append(Line2D([0], [0], color=color, lw=1.5, label=label))

    if not legend_handles:
        raise RuntimeError("Not enough rows per experiment to draw a joint KDE plot.")

    x_min = min(docking_min, docking_max)
    x_max = max(docking_min, docking_max)
    grid.ax_joint.set_xlim(-x_max, -x_min)
    grid.ax_joint.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{-value:g}"))
    grid.ax_marg_x.set_xlim(grid.ax_joint.get_xlim())

    y_min = plot_df["green_score"].min()
    y_max = plot_df["green_score"].max()
    y_pad = max((y_max - y_min) * 0.05, 0.02)
    grid.ax_joint.set_ylim(max(0.0, y_min - y_pad), min(1.0, y_max + y_pad))
    grid.ax_marg_y.set_ylim(grid.ax_joint.get_ylim())

    grid.ax_joint.set_xlabel("Docking score")
    grid.ax_joint.set_ylabel("Green chemistry score")
    grid.ax_joint.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
    grid.ax_joint.legend(
        handles=legend_handles,
        fontsize=7,
        loc="lower left",
        frameon=True,
        title="Experiment",
    )
    grid.figure.suptitle("Green Score vs Docking KDE Across Experiments", y=1.02)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    grid.figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(grid.figure)


