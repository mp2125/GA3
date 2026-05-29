import numpy as np
from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX
from optimiser import fullSolver
from parameters import *
from weight_limits import getWeight
from plotterFuncs import plot_latex

WEIGHT_LIMIT = 1.1  # kg

# ── Nominal (fixed) values ────────────────────────────────────────────────────
NOMINAL = {
    "length":       0.281,
    "tubes":        8,
    "baffles":      12,
    "shell_passes": 2,
}

# ── Sweep ranges for each parameter ──────────────────────────────────────────
SWEEPS = {
    "length":       np.linspace(0.25, 0.80, 101),
    "tubes":        np.arange(4, 21),       # 4 to 20 inclusive
    "baffles":      np.arange(8, 21),       # 8 to 20 inclusive
    "shell_passes": np.arange(1, 4),        # 1 to 3 inclusive
}

XLABELS = {
    "length":       "Tube Length (m)",
    "tubes":        "Number of Tubes",
    "baffles":      "Number of Baffles",
    "shell_passes": "Number of Shell Passes",
}


def compute_pitch(tubes, ds, do):
    """Approximate pitch from tube count and geometry."""
    phi = tubes * (do / ds) ** 2
    return do / 2 * (2 * np.pi / (3 ** 0.5 * phi)) ** 0.5


def find_weight_limit_x(x_vals, weights):
    """
    Interpolate the x value where weight crosses WEIGHT_LIMIT.
    Returns None if the limit is never crossed within the sweep range.
    """
    for i in range(len(weights) - 1):
        w0, w1 = weights[i], weights[i + 1]
        if (w0 - WEIGHT_LIMIT) * (w1 - WEIGHT_LIMIT) < 0:
            # Linear interpolation between the two bracketing points
            t = (WEIGHT_LIMIT - w0) / (w1 - w0)
            return x_vals[i] + t * (x_vals[i + 1] - x_vals[i])
    return None


def run_sensitivity(use_eNTU=True, scale=1.0, font_size=11,
                    save_dir="Report-Files"):
    """
    Run a one-at-a-time sensitivity analysis over length, tubes, baffles,
    and shell_passes.  Each parameter is swept across its range while the
    others are held at their nominal values.  One PDF plot is saved per
    parameter.

    Parameters
    ----------
    use_eNTU : bool
        If True, plot Q from the eNTU method; otherwise use LMTD.
    scale : float
        Passed to plot_latex to control figure size.
    font_size : float
        Passed to plot_latex; set to match your LaTeX document font size.
    save_dir : str
        Directory to save PDFs into (must already exist).
    """
    method = "eNTU" if use_eNTU else "LMTD"

    # Pre-compute the design point Q using nominal values
    nominal_pitch = compute_pitch(int(NOMINAL["tubes"]), ds, do)
    nominal_hx = HX(
        int(NOMINAL["tubes"]), int(NOMINAL["baffles"]),
        NOMINAL["length"], nominal_pitch,
        False, int(NOMINAL["shell_passes"]), int(NOMINAL["shell_passes"])
    )
    nominal_Q = fullSolver(nominal_hx)[1 if use_eNTU else 0]

    for param, sweep_values in SWEEPS.items():
        print(f"\nSweeping {param} ({len(sweep_values)} points)...")
        x_vals   = []
        q_vals   = []
        w_vals   = []   # weight at each sweep point

        for val in sweep_values:
            p = dict(NOMINAL)
            p[param] = val

            length       = p["length"]
            tubes        = int(p["tubes"])
            baffles      = int(p["baffles"])
            shell_passes = int(p["shell_passes"])
            pitch        = compute_pitch(tubes, ds, do)

            hx = HX(tubes, baffles, length, pitch,
                    False, shell_passes, shell_passes)

            Q_LMTD, Q_eNTU = fullSolver(hx)
            Q = Q_eNTU if use_eNTU else Q_LMTD
            W = getWeight(hx)

            x_vals.append(val)
            q_vals.append(Q)
            w_vals.append(W)

            print(f"  {param} = {val:.4g}  →  Q = {Q:.4f} kW,  W = {W:.4f} kg")

        # Design point — x is the nominal value, y is the nominal Q
        design_x = NOMINAL[param]
        design_y = nominal_Q

        # Weight limit crossing x
        weight_limit_x = find_weight_limit_x(x_vals, w_vals)
        if weight_limit_x is None:
            print(f"  Weight limit ({WEIGHT_LIMIT} kg) not crossed in sweep range.")

        save_path = f"{save_dir}/{param}_{method}.pdf"
        plot_latex(
            x_vals, q_vals,
            xlabel=XLABELS[param],
            ylabel=f"Heat Transfer, $Q_\\mathrm{{{method}}}$ (kW)",
            scale=scale,
            font_size=font_size,
            design_x=design_x,
            design_y=design_y,
            weight_limit_x=weight_limit_x,
            save_path=save_path,
        )
        print(f"  Saved → {save_path}")


if __name__ == "__main__":
    run_sensitivity(use_eNTU=True, scale=1.0, font_size=11)