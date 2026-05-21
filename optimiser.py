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

def fullSolver( hx,
                Tin_cold = 20,
                Tin_hot = 60,
                ):
    
    shell_mass, tube_mass = solve_mass_flows(hx)
    
    vel_shell = hx.shell_side_velocity(shell_mass)
    vel_tube = hx.tube_side_velocity(tube_mass)

    Re_shell = hx.shell_side_reynolds_number(vel_shell)
    Re_tube = hx.tube_side_reynolds_number(vel_tube)

    shape = 'square' if hx.is_square_layout else 'triangle'
    H = handoutThermalCoefficient(Re_shell, Re_tube, hx.baffle_spacing, shape)

    if hx.tube_passes == hx.shell_passes:
        cfg = 'counterflow'
    else:
        cfg = 'N-2N'

    Tout_cold_LMTD, Tout_hot_LMTD, Q_LMTD = tempIteratorLMTD(hx.tube_length,
                                                             hx.number_of_tubes,
                                                             shell_mass,
                                                             tube_mass,
                                                             H,
                                                             Tin_cold,
                                                             Tin_hot,
                                                             hx.tube_passes,
                                                             )

    Tout_cold_eNTU, Tout_hot_eNTU, Q_eNTU, eps = eNTUProcessor(hx.tube_length,
                                                          hx.number_of_tubes,
                                                          shell_mass,
                                                          tube_mass,
                                                          H,
                                                          Tin_cold,
                                                          Tin_hot,
                                                          config=cfg,
                                                          N=hx.shell_passes)

    return Q_LMTD, Q_eNTU/1.7

if __name__ == "__main__":
    from previous_HXs import hxs, heat_transfers

    lengths = []
    errors = []
    baffle_numbers = []
    tubes = []
    tubePasses = []
    shellPasses = []
    relPasses = []

    for n in range(len(hxs)):

        Q_calc = fullSolver(hxs[n])[1]
        Q_exp = heat_transfers[n] * 1e3

        error = ((Q_calc-Q_exp)/Q_exp)
        if error > 1.0: print(n)

        errors.append(error)
        lengths.append(hxs[n].tube_length)
        baffle_numbers.append(hxs[n].number_of_baffles)
        tubes.append(hxs[n].number_of_tubes)
        tubePasses.append(hxs[n].tube_passes)
        shellPasses.append(hxs[n].shell_passes)
        relPasses.append(hxs[n].tube_passes/hxs[n].shell_passes)

    import matplotlib.pyplot as plt
    import numpy as np
    # plt.scatter(tubePasses, errors, label='tube passes')
    # plt.scatter(lengths, errors, label='lengths')
    # plt.scatter(tubes, errors, label='tubes')
    plt.scatter(baffle_numbers, errors, label='baffles')
    print(f'mean: {np.average(errors):.3f}, sd: {(np.var(errors))**0.5:.3f}')
    plt.ylabel('error')
    plt.legend()
    plt.show()



    