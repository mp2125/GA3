import numpy as np
from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX
from optimiser import fullSolver
from parameters import *
from plotterFuncs import plot_latex

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
    "shell_passes": "Number of Shell Passes (N-N)",
}


def compute_pitch(tubes, ds, do):
    """Approximate pitch from tube count and geometry."""
    phi = tubes * (do / ds) ** 2
    return do / 2 * (2 * np.pi / (3 ** 0.5 * phi)) ** 0.5


def run_sensitivity(use_eNTU=False, scale=1.0, font_size=11,
                    save_dir="Sensitivities"):
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
    q_index = 1 if use_eNTU else 0   # fullSolver returns (Q_LMTD, Q_eNTU)

    for param, sweep_values in SWEEPS.items():
        print(f"\nSweeping {param} ({len(sweep_values)} points)...")
        x_vals = []
        q_vals = []

        for val in sweep_values:
            # Build the parameter set for this evaluation
            p = dict(NOMINAL)
            p[param] = val

            length       = p["length"]
            tubes        = int(p["tubes"])
            baffles      = int(p["baffles"])
            shell_passes = int(p["shell_passes"])

            pitch = compute_pitch(tubes, ds, do)

            hx = HX(tubes, baffles, length, pitch,
                    False, shell_passes, shell_passes)

            Q_LMTD, Q_eNTU = fullSolver(hx)
            Q = Q_eNTU if use_eNTU else Q_LMTD

            x_vals.append(val)
            q_vals.append(Q)

            print(f"  {param} = {val:.4g}  →  Q_{method} = {Q:.4f} kW")

        save_path = f"{save_dir}/{param}_{method}.pdf"
        plot_latex(
            x_vals, q_vals,
            xlabel=XLABELS[param],
            ylabel=f"Heat Transfer, $Q$ (kW)",
            scale=scale,
            font_size=font_size,
            save_path=save_path,
        )
        print(f"  Saved → {save_path}")


if __name__ == "__main__":
    run_sensitivity(use_eNTU=True, scale=1.0, font_size=11)