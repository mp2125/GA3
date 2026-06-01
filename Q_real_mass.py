"""
    number_of_tubes,
    number_of_baffles,
    tube_length,
    tube_pitch,
    is_square_layout,
    tube_passes=2,
    shell_passes=1

"""
from thermal_analysis import handoutThermalCoefficient, eNTUProcessor, tempIteratorLMTD
from hydraulic_analysis import solve_mass_flows
from parameters import *

def realMassSolver( hx, mdotcold, mdothot,
                Tin_cold = 20,
                Tin_hot = 60,
                ):
    
    vel_shell = hx.shell_side_velocity(mdotcold)
    vel_tube = hx.tube_side_velocity(mdothot)

    Re_shell = hx.shell_side_reynolds_number(vel_shell)
    Re_tube = hx.tube_side_reynolds_number(vel_tube)

    if hx.tube_passes == hx.shell_passes:
        cfg = 'counterflow'
        # cfg = 'parallelflow'
        tube_length = hx.tube_length * hx.shell_passes
        number_of_tubes = hx.number_of_tubes / hx.shell_passes
        shell_passes = 1
        tube_passes = 1
    else:
        cfg = 'N-2N'
        tube_length = hx.tube_length
        number_of_tubes = hx.number_of_tubes
        shell_passes = hx.shell_passes
        tube_passes = hx.tube_passes

    shape = 'square' if hx.is_square_layout else 'triangle'
    H = handoutThermalCoefficient(Re_shell, Re_tube, hx.baffle_spacing, shape)

    lengthCorrection = (hx.number_of_baffles+2)*0.0015
    Tout_cold_LMTD, Tout_hot_LMTD, Q_LMTD = tempIteratorLMTD((tube_length - lengthCorrection),
                                                             number_of_tubes,
                                                             mdotcold,
                                                             mdothot,
                                                             H,
                                                             Tin_cold,
                                                             Tin_hot,
                                                             shell_passes,
                                                             tube_passes,
                                                             )

    Tout_cold_eNTU, Tout_hot_eNTU, Q_eNTU, eps = eNTUProcessor((tube_length - lengthCorrection),
                                                          number_of_tubes,
                                                          mdotcold,
                                                          mdothot,
                                                          H,
                                                          Tin_cold,
                                                          Tin_hot,
                                                          config=cfg,
                                                          N=shell_passes)


    return Q_LMTD, Q_eNTU

if __name__ == "__main__":
    from previous_HXs import hxs, heat_transfers, temperatures, mass_flows

    lengths = []
    errors = []
    baffle_numbers = []
    tubes = []
    tubePasses = []
    shellPasses = []
    relPasses = []
    Qs = []

    for n in range(len(hxs)):

        Q_calc = realMassSolver(hxs[n], *mass_flows[n], temperatures[n][0],temperatures[n][2])[1]
        Q_exp = heat_transfers[n] * 1e3

        error = ((Q_calc-Q_exp)/Q_exp)*100
        if abs(error) > 30: print(n)
        Qs.append(Q_calc)

        errors.append(error)
        lengths.append(hxs[n].tube_length)
        baffle_numbers.append(hxs[n].number_of_baffles)
        tubes.append(hxs[n].number_of_tubes)
        tubePasses.append(hxs[n].tube_passes)
        shellPasses.append(hxs[n].shell_passes)
        relPasses.append(hxs[n].tube_passes/hxs[n].shell_passes)


    from plotterFuncs import plot_latex

    plot_latex(
        baffle_numbers, errors,
        xlabel="Number of Baffles",
        ylabel=f"Error in Heat Transfer (%)",
        scale=1.0,
        font_size=11,
        save_path="Report-Files/erros_real_mass.pdf",
        linestyle='none',
        marker='o'
    )

    import numpy as np
    print(np.sqrt(np.mean(errors)**2))



    