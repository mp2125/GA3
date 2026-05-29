import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def plot_latex(
    x,
    y,
    xlabel="x",
    ylabel="y",
    title=None,
    scale=1.0,
    font_size=11,
    save_path=None,
    design_x=None,
    design_y=None,
    weight_limit_x=None,
):
    """
    Plot x vs y in a style suitable for LaTeX reports.

    Uses Times New Roman font throughout. The `scale` parameter controls the
    physical size of the figure (and line/tick thicknesses so they look
    proportionally correct at that size), while font sizes stay fixed so they
    match your document's body text regardless of figure size.

    Parameters
    ----------
    x : array-like
        Horizontal data values.
    y : array-like
        Vertical data values.
    xlabel : str
        Label for the x-axis.
    ylabel : str
        Label for the y-axis.
    title : str or None
        Optional figure title.
    scale : float
        Controls the physical size of the figure.
          1.0  → single-column width (~3.5 in)
          1.5  → larger / 1.5-column width (~5.25 in)
          2.0  → full double-column width (~7 in)
        Line widths, tick sizes, and spine thicknesses scale with this so
        the figure looks correct at every size. Font sizes do NOT scale.
    font_size : float
        Base font size in points for all text (tick labels, axis labels,
        legend). Set this to match your LaTeX document's body font size
        (e.g. 10 for a 10 pt document, 11 or 12 otherwise).
        Title (if used) is rendered at font_size + 1 pt.
    save_path : str or None
        If given, save the figure to this path. Use .pdf for lossless
        vector embedding in LaTeX.
    design_x : float or None
        x-value of the design point. A marker is drawn on the curve at
        this location. design_y must also be provided.
    design_y : float or None
        y-value of the design point.
    weight_limit_x : float or None
        x-value at which the weight crosses 1.1 kg. A vertical dashed line
        labelled "Weight Limit" is drawn at this position.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes objects.
    """

    # ── Figure size scales with `scale` ───────────────────────────────────────
    FIG_WIDTH  = 3.5 * scale   # inches
    FIG_HEIGHT = 2.8 * scale   # inches

    # ── Line/tick geometry scales with `scale` ────────────────────────────────
    LINE_WIDTH      = 1.5 * scale
    TICK_WIDTH      = 0.8 * scale
    TICK_LENGTH     = 4.0 * scale
    AXIS_LINE_WIDTH = 0.8 * scale
    MARKER_SIZE     = 6.0 * scale

    # ── Font sizes are FIXED — set by font_size only ──────────────────────────
    FONT_SIZE  = font_size
    LABEL_SIZE = font_size
    TITLE_SIZE = font_size + 1

    # ── Global rcParams ───────────────────────────────────────────────────────
    mpl.rcParams.update({
        # Font
        "font.family":        "serif",
        "font.serif":         ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset":   "stix",   # STIX matches Times New Roman in equations
        "font.size":          FONT_SIZE,
        "axes.labelsize":     LABEL_SIZE,
        "axes.titlesize":     TITLE_SIZE,
        "xtick.labelsize":    FONT_SIZE,
        "ytick.labelsize":    FONT_SIZE,
        "legend.fontsize":    FONT_SIZE,
        # Lines & spines
        "lines.linewidth":    LINE_WIDTH,
        "axes.linewidth":     AXIS_LINE_WIDTH,
        # Ticks
        "xtick.major.width":  TICK_WIDTH,
        "ytick.major.width":  TICK_WIDTH,
        "xtick.minor.width":  TICK_WIDTH * 0.7,
        "ytick.minor.width":  TICK_WIDTH * 0.7,
        "xtick.major.size":   TICK_LENGTH,
        "ytick.major.size":   TICK_LENGTH,
        "xtick.minor.size":   TICK_LENGTH * 0.6,
        "ytick.minor.size":   TICK_LENGTH * 0.6,
        "xtick.direction":    "out",
        "ytick.direction":    "out",
        "xtick.top":          False,
        "ytick.right":        False,
    })

    # ── Figure & axes ─────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))

    ax.plot(x, y, linewidth=LINE_WIDTH, color="cornflowerblue")

    # ── Design point marker ───────────────────────────────────────────────────
    if design_x is not None and design_y is not None:
        ax.plot(
            design_x, design_y,
            marker="o",
            markersize=MARKER_SIZE,
            color="black",
            markerfacecolor="white",
            markeredgewidth=LINE_WIDTH * 0.8,
            linestyle="none",
            label="Design point",
            zorder=5,
        )
        ax.legend(frameon=False)

    # ── Weight limit vertical line ────────────────────────────────────────────
    if weight_limit_x is not None:
        ax.axvline(
            weight_limit_x,
            linestyle="--",
            linewidth=LINE_WIDTH * 0.8,
            color="black",
        )
        ax.text(
            weight_limit_x, ax.get_ylim()[1],
            "  Weight Limit",
            va="top",
            ha="left",
            fontsize=FONT_SIZE,
            fontstyle="italic",
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    fig.tight_layout()

    # ── Save ──────────────────────────────────────────────────────────────────
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"Figure saved to: {save_path}")
    else: fig.show()

    return fig, ax


# ── Quick demo ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    x = np.linspace(0, 2 * np.pi, 300)
    y = np.sin(x)

    fig, ax = plot_latex(
        x, y,
        xlabel=r"$x$ (rad)",
        ylabel=r"$\sin(x)$",
        scale=1.0,
        font_size=11,
        design_x=np.pi / 2,
        design_y=1.0,
        weight_limit_x=4.5,
    )
    plt.show()