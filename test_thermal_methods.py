import numpy as np
from parameters import *

from hydraulic_analysis import solve_mass_flows
from thermal_analysis import handoutThermalCoefficient, tempIteratorLMTD, eNTUProcessor

from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX

from previous_HXs import hxs, mass_flows, temperatures, heat_transfers


# ------------------------------------------------------------
# SOLVE OPERATING MASS FLOWS
# ------------------------------------------------------------

for hx, mdots, temps, Q_real in zip(hxs, mass_flows, temperatures, heat_transfers):
    mass_flow_cold , mass_flow_hot = mdots  # Compares using the actual mass flows from the tests to isolate validity of thermal analysis

    shell_velocity = hx.shell_side_velocity(mass_flow_cold)
    tube_velocity = hx.tube_side_velocity(mass_flow_hot)

    ReSh = hx.shell_side_reynolds_number(shell_velocity)
    ReTu = hx.tube_side_reynolds_number(tube_velocity)

    H = handoutThermalCoefficient(ReSh, ReTu, hx.baffle_spacing)

    print(
        f"mdot1={mass_flow_cold:.3f} kg/s  "
        f"mdot2={mass_flow_hot:.3f} kg/s  "
        f"ReSh={ReSh:.0f}  "
        f"ReTu={ReTu:.0f}"
    )

    # ------------------------------------------------------------
    # Test Results
    # ------------------------------------------------------------

    print(
        f"For Real          "
        f"T1_out={temps[1]:.1f}°C  "
        f"T2_out={temps[3]:.1f}°C  "
        f"Q={Q_real*1e3:.1f} W"
    )

    # ------------------------------------------------------------
    # LMTD METHOD
    # ------------------------------------------------------------

    T1_out, T2_out, Q = tempIteratorLMTD(
        hx.tube_length,
        hx.number_of_tubes,
        mass_flow_cold,
        mass_flow_hot,
        H,
        T1in=temps[0],
        T2in=temps[2],
        passes=hx.tube_passes
    )

    print(
        f"For LMTD          "
        f"T1_out={T1_out:.1f}°C  "
        f"T2_out={T2_out:.1f}°C  "
        f"Q={Q:.1f} W"
    )


    # ------------------------------------------------------------
    # eNTU METHOD
    # ------------------------------------------------------------

    T1_out, T2_out, Q, eps = eNTUProcessor(
        hx.tube_length,
        hx.number_of_tubes,
        mass_flow_cold,
        mass_flow_hot,
        H,
        T1_in=temps[0],
        T2_in=temps[2],
        N=hx.shell_passes
    )

    print(
        f"For eNTU ε={eps:.4f} "
        f"T1_out={T1_out:.1f}°C  "
        f"T2_out={T2_out:.1f}°C  "
        f"Q={Q:.1f} W"
    )
    print()