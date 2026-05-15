import numpy as np
from parameters import *

from thermal_analysis import handoutThermalCoefficient, tempIteratorLMTD, eNTUProcessor

config = '1-2'

length = 0.35
tubes = 13

T1_out,T2_out,Q,eps = eNTUProcessor(length,tubes,mdot1,mdot2)

print(f"{config:<30}  ε={eps:.4f}  T1_out={T1_out:.1f}°C  T2_out={T2_out:.1f}°C  Q={Q:.1f}W")



