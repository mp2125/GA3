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
                                                          N=1)

    return Q_LMTD, Q_eNTU

# from previous_HXs import hxs, heat_transfers

# lengths = []
# errors = []
# baffle_numbers = []
# tubes = []

# for n in range(len(hxs)):
#     Q_calc = fullSolver(hxs[n])[1]
#     Q_exp = heat_transfers[n] * 1e3

#     error = ((Q_calc-Q_exp)/Q_exp)
#     if error > 1.0: print(n)

#     errors.append(error)
#     lengths.append(hxs[n].tube_length)
#     baffle_numbers.append(hxs[n].number_of_baffles)
#     tubes.append(hxs[n].number_of_tubes)

# import matplotlib.pyplot as plt
# plt.scatter(lengths, errors)
# plt.ylabel('error')
# plt.show()


import numpy as np
from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX
import matplotlib.pyplot as plt

length_varying = []
lengths = np.linspace(0.05,0.35,21)
for length in lengths:
    max_tubes = int(min(0.9069 * (ds/do)**2, totalCu/length))
tubes = np.arange(5,max_tubes)
print(tubes)
baffles = np.arange(1,10)
pitches = np.linspace(do*2,do*4,21)
shell_passes = [1,2,3,4]
shapes = [True,False]

Qmax_LMTD = 0
Qmax_eNTU = 0
params_LMTD = []
params_eNTU = []

output = []

total = len(tubes) * len(baffles) * len(lengths) * len(pitches) * len(shapes) * len(shell_passes)
i = 0

for tube in tubes:
    for baffle in baffles:
        for length in lengths:
            for pitch in pitches:
                for shape in shapes:
                    for shell_pass in shell_passes:
                        hx = HX(tube,baffle,length,pitch,shape,shell_pass*2,shell_pass)
                        Q_LMTD,Q_eNTU = fullSolver(hx)
                        if Q_LMTD > Qmax_LMTD:
                            Qmax_LMTD = Q_LMTD
                            params_LMTD = [tube,baffle,length,pitch,shape,shell_pass*2,shell_pass]
                        if Q_eNTU > Qmax_eNTU:
                            Qmax_eNTU = Q_eNTU
                            params_eNTU = [tube,baffle,length,pitch,shape,shell_pass*2,shell_pass]
                        i+=1
                        percent = i / total
                        bar = "#" * int(percent * 40)
                        spaces = " " * (40 - len(bar))
                        output.append([tube,baffle,length,pitch,shape,shell_pass*2,shell_pass,Q_LMTD,Q_eNTU])
                        print(f"\r[{bar}{spaces}] {percent:.0%}")

print(f'LMTD: {Qmax_LMTD}W with params: {params_LMTD}')
print(f'eNTU: {Qmax_eNTU}W with params: {params_eNTU}')

with open("outputOptimisation.txt", "w") as f:
    for item in output:
        f.write(str(item) + "\n")


    