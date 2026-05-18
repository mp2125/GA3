import numpy as np
from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX
import matplotlib.pyplot as plt
from optimiser import fullSolver
from parameters import *
from weight_limits import getWeight

length_varying = []
lengths = np.linspace(0.05,0.325,21)
tubeNumber = 0
for length in lengths:
    tubeNumber += int(min(0.9069 * (ds/do)**2, totalCu/length)) - 4

baffles = np.arange(1,10)
shell_passes = [1,2,3,4]
shapes = [False]

Qmax_LMTD = 0
Qmax_eNTU = 0
params_LMTD = []
params_eNTU = []

output = []

total = tubeNumber * len(baffles) * len(shapes) * len(shell_passes)
i = 0

for length in lengths:
    max_packing = 0.9069 * (ds/do)**2 
    max_tubes = int(min(max_packing, totalCu/length))
    tubes = np.arange(5,max_tubes)
    for baffle in baffles:
        for tube in tubes:
            # approximate pitch distance
            phi = tube * (do/ds)**2
            pitch = do/2 * (2*pi/(3**0.5 * phi))**0.5
            if pitch > 2*do: # TO DO: this *2 can probably be reduced, but reluctant to without pressures being more certain
                for shape in shapes:
                    for shell_pass in shell_passes:
                        i+=1
                        hx = HX(tube,baffle,length,pitch,shape,shell_pass*2,shell_pass)
                        if getWeight(hx) < maxWeight:
                            Q_LMTD,Q_eNTU = fullSolver(hx)
                            if Q_LMTD > Qmax_LMTD:
                                Qmax_LMTD = Q_LMTD
                                params_LMTD = [tube,baffle,length,pitch,shape,shell_pass*2,shell_pass]
                            if Q_eNTU > Qmax_eNTU:
                                Qmax_eNTU = Q_eNTU
                                params_eNTU = [tube,baffle,length,pitch,shape,shell_pass*2,shell_pass]

                            percent = i / total
                            bar = "#" * int(percent * 40)
                            spaces = " " * (40 - len(bar))
                            output.append([tube,baffle,length,pitch,shape,shell_pass*2,shell_pass,Q_LMTD,Q_eNTU])
                        print(f"\r[{bar}{spaces}] {percent:.0%}")
            else:
                i+= len(shapes) * len(shell_passes)
                percent = i / total
                bar = "#" * int(percent * 40)
                spaces = " " * (40 - len(bar))
                output.append([tube,baffle,length,pitch,shape,shell_pass*2,shell_pass,Q_LMTD,Q_eNTU])

print(f'LMTD: {Qmax_LMTD}W with params: {params_LMTD} and weight')
print(f'eNTU: {Qmax_eNTU}W with params: {params_eNTU}')
print(f'{total} values checked')

with open("outputOptimisation.txt", "w") as f:
    for item in output:
        f.write(str(item) + "\n")