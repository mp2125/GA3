import numpy as np
from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX
from optimiser import fullSolver
from parameters import *
from weight_limits import getWeight
from plotterFuncs import plot_latex
import matplotlib.pyplot as plt
import matplotlib as mpl

WEIGHT_LIMIT = 1.1  # kg

# ── Nominal (fixed) values ────────────────────────────────────────────────────
NOMINAL = {
    "length":       0.281,
    "tubes":        8,
    "baffles":      12,
    "shell_passes": 2,
}

# ── Sweep ranges ──────────────────────────────────────────────────────────────
SWEEPS = {
    "length":       np.linspace(0.15, 0.35, 101),
    "tubes":        np.arange(4, 12),
    "baffles":      np.arange(5, 25),
    #"shell_passes": np.arange(1, 4),
}

LABELS = {
    "length":       "Tube Length",
    "tubes":        "No. Tubes",
    "baffles":      "No. Baffles",
    "shell_passes": "Shell Passes",
}

COLORS = {
    "length":       "cornflowerblue",
    "tubes":        "tomato",
    "baffles":      "mediumseagreen",
    "shell_passes": "darkorange",
}


def compute_pitch(tubes, ds, do):
    phi = tubes * (do / ds) ** 2
    return do / 2 * (2 * np.pi / (3 ** 0.5 * phi)) ** 0.5


def plot_q_vs_weight(use_eNTU=True, scale=1.0, font_size=11,
                     save_path="Report-Files/q_vs_weight.pdf"):
    """
    Plot heat transfer vs weight for each parameter sweep as a separate
    labelled curve on a single axes. The design point and weight limit
    are also marked.

    Parameters
    ----------
    use_eNTU : bool
        If True, use Q from eNTU method; otherwise LMTD.
    scale : float
        Controls figure size (see plot_latex).
    font_size : float
        Fixed font size in points.
    save_path : str
        Output PDF path.
    """
    method = "eNTU" if use_eNTU else "LMTD"

    # ── Pre-compute nominal Q and W ───────────────────────────────────────────
    nominal_pitch = compute_pitch(int(NOMINAL["tubes"]), ds, do)
    nominal_hx = HX(
        int(NOMINAL["tubes"]), int(NOMINAL["baffles"]),
        NOMINAL["length"], nominal_pitch,
        False, int(NOMINAL["shell_passes"]), int(NOMINAL["shell_passes"])
    )
    nominal_Q = fullSolver(nominal_hx)[1 if use_eNTU else 0]
    nominal_W = getWeight(nominal_hx)

    # ── Style (mirrors plot_latex) ────────────────────────────────────────────
    FIG_WIDTH       = 3.5 * scale
    FIG_HEIGHT      = 2.8 * scale
    LINE_WIDTH      = 1.5 * scale
    TICK_WIDTH      = 0.8 * scale
    TICK_LENGTH     = 4.0 * scale
    AXIS_LINE_WIDTH = 0.8 * scale
    MARKER_SIZE     = 6.0 * scale
    FONT_SIZE       = font_size

    mpl.rcParams.update({
        "font.family":        "serif",
        "font.serif":         ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset":   "stix",
        "font.size":          FONT_SIZE,
        "axes.labelsize":     FONT_SIZE,
        "axes.titlesize":     FONT_SIZE + 1,
        "xtick.labelsize":    FONT_SIZE,
        "ytick.labelsize":    FONT_SIZE,
        "legend.fontsize":    FONT_SIZE,
        "lines.linewidth":    LINE_WIDTH,
        "axes.linewidth":     AXIS_LINE_WIDTH,
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

    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))

    # ── Plot each parameter sweep ─────────────────────────────────────────────
    for param, sweep_values in SWEEPS.items():
        print(f"\nSweeping {param}...")
        w_vals = []
        q_vals = []

        for val in sweep_values:
            p = dict(NOMINAL)
            p[param] = val

            tubes        = int(p["tubes"])
            baffles      = int(p["baffles"])
            length       = p["length"]
            shell_passes = int(p["shell_passes"])
            pitch        = compute_pitch(tubes, ds, do)

            hx = HX(tubes, baffles, length, pitch,
                    False, shell_passes, shell_passes)

            Q = fullSolver(hx)[1 if use_eNTU else 0]
            W = getWeight(hx)

            w_vals.append(W)
            q_vals.append(Q)
            print(f"  {param} = {val:.4g}  →  Q = {Q:.4f} kW,  W = {W:.4f} kg")

        ax.plot(w_vals, q_vals, '-',
                linewidth=LINE_WIDTH,
                color=COLORS[param],
                label=LABELS[param])

    # ── Design point ──────────────────────────────────────────────────────────
    ax.plot(
        nominal_W, nominal_Q,
        marker="o",
        markersize=MARKER_SIZE,
        color="black",
        markerfacecolor="white",
        markeredgewidth=LINE_WIDTH * 0.8,
        linestyle="none",
        label="Design point",
        zorder=5,
    )

    # ── Weight limit vertical line ────────────────────────────────────────────
    ax.axvline(
        WEIGHT_LIMIT,
        linestyle="--",
        linewidth=LINE_WIDTH * 0.8,
        color="black",
    )
    ylims = ax.get_ylim()
    ax.text(
        WEIGHT_LIMIT, ylims[0] + 0.6 * (ylims[1] - ylims[0]),
        "  Weight Limit",
        va="top",
        ha="left",
        fontsize=FONT_SIZE,
        fontstyle="italic",
    )

    ax.set_xlabel("Weight (kg)")
    ax.set_ylabel(f"Heat Transfer, $Q$ (kW)")
    ax.legend(frameon=False, edgecolor="black", fancybox=False, loc="upper left")

    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        print(f"\nSaved → {save_path}")
    else:
        fig.show()

    return fig, ax


if __name__ == "__main__":
    plot_q_vs_weight(use_eNTU=True, scale=1.2, font_size=11)