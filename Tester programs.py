import numpy as np
from parameters import *

from hydraulic_analysis import get_Re_shell, get_Re_tube, solve_mass_flows
from thermal_analysis import handoutThermalCoefficient, tempIteratorLMTD, eNTUProcessor

config = '1-2'
length = 0.35
tubes = 13
pitch = 0.014
baffles = 9
square = False

mdot1,mdot2 = solve_mass_flows(length,pitch,tubes,baffles,square)
ReSh = get_Re_shell(mdot1,length,pitch,baffles)
ReTu = get_Re_tube(mdot2,length,tubes,baffles)
H = handoutThermalCoefficient(ReSh,ReTu)

print(f"mdot1={mdot1:.3f} kg/s  mdot2={mdot2:.3f} kg/s  ReSh={ReSh:.0f}  ReTu={ReTu:.0f}")

T1_out,T2_out,Q = tempIteratorLMTD(length,tubes,mdot1,mdot2,H)
print(f"For LMTD          T1_out={T1_out:.1f}°C  T2_out={T2_out:.1f}°C  Q={Q:.1f} W")

T1_out,T2_out,Q,eps = eNTUProcessor(length,tubes,mdot1,mdot2,H)
print(f"For eNTU ε={eps:.4f} T1_out={T1_out:.1f}°C  T2_out={T2_out:.1f}°C  Q={Q:.1f} W")