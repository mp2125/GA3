import numpy as np
from parameters import *
from thermal_analysis import handoutThermalCoefficient

pipes = 13
passes = 1
dividers = 0
length = 0.35

segmentLength = length/(dividers+1)
Aheat = np.pi * di * segmentLength

T1in = Tcold_in
T2in = Thot_in

mdot1 = 1
mdot2 = 1

C1 = mdot1 * c_p
C2 = mdot2 * c_p

Cmin = min(C1,C2)
