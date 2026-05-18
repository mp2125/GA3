import numpy as np
from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX
import matplotlib.pyplot as plt
from optimiser import fullSolver
from parameters import *
from weight_limits import getWeight

length_varying = []
lengths = np.linspace(0.05,0.325,21)
for length in lengths:
    max_tubes = int(min(0.9069 * (ds/do)**2, totalCu/length))
tubes = np.arange(5,max_tubes)
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
                        if getWeight(hx) < maxWeight:
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