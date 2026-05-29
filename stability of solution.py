import numpy as np
from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX
import matplotlib.pyplot as plt
from optimiser import fullSolver
from parameters import *
from weight_limits import getWeight

length_varying = []
lengths = [0.281] #np.linspace(0.250,0.80,101)
tubeNumber = 8
baffles = [12] # np.arange(8,20)
shell_passes = [1,2,3]
shapes = [False]

# lengths = [0.2426]
# baffles = [11]
# shell_passes = [2]
# shapes = [False]
# tubeNumber = 10

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
    # tubes = np.arange(5,max_tubes)
    tubes = [8] #np.arange(4,20)
    for baffle in baffles:
        for tube in tubes:
            # approximate pitch distance
            phi = tube * (do/ds)**2
            pitch = do/2 * (2*pi/(3**0.5 * phi))**0.5
            if 1:
            # if pitch > do*2: # TO DO: this *2 can probably be reduced, but reluctant to without pressures being more certain
                for shape in shapes:
                    for shell_pass in shell_passes:
                        hx = HX(tube,baffle,length,pitch,shape,shell_pass*2,shell_pass)

                        # if getWeight(hx) < maxWeight:
                        if True:
                            Q_LMTD,Q_eNTU = fullSolver(hx)
                            if Q_LMTD > Qmax_LMTD:
                                Qmax_LMTD = Q_LMTD
                                params_LMTD = [tube,baffle,length,pitch,shape,shell_pass*2,shell_pass]
                            if Q_eNTU > Qmax_eNTU:
                                Qmax_eNTU = Q_eNTU
                                params_eNTU = [tube,baffle,length,pitch,shape,shell_pass*2,shell_pass]

                            output.append([tube,baffle,length,pitch,shape,shell_pass*2,shell_pass,Q_LMTD,Q_eNTU])
                        
                        hx = HX(tube,baffle,length,pitch,shape,shell_pass,shell_pass)
                        #if getWeight(hx) < maxWeight:
                        # if True:
                        #     Q_LMTD,Q_eNTU = fullSolver(hx)
                        #     if Q_LMTD > Qmax_LMTD:
                        #         Qmax_LMTD = Q_LMTD
                        #         params_LMTD = [tube,baffle,length,pitch,shape,shell_pass,shell_pass]
                        #     if Q_eNTU > Qmax_eNTU:
                        #         Qmax_eNTU = Q_eNTU
                        #         params_eNTU = [tube,baffle,length,pitch,shape,shell_pass,shell_pass]

                        #     output.append([tube,baffle,length,pitch,shape,shell_pass,shell_pass,Q_LMTD,Q_eNTU])

            else:
                output.append([tube,baffle,length,pitch,shape,shell_pass*2,shell_pass,Q_LMTD,Q_eNTU])

x = [hx[6] for hx in output]
y = [hx[8] for hx in output]

import matplotlib.pyplot as plt

plt.scatter(x,y)
plt.show()