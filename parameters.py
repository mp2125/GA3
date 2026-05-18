# parameters.py
from math import pi
# temperatures
Tcold_in = 20.8
Thot_in = 57.7

# physical parameters (for units, check handout)
c_p = 4179
mu = 6.51e-4
rho_w = 990
Pr = 4.31
k_w = 0.632
k_tube = 380

# part constraints (all m)
Lmax = 400e-3
di = 6e-3
do = 8e-3
ds = 64e-3
t_splitter = 1.5e-3

Ai_water = pi*di**2/4
Ao_water = pi*do**2/4

d_noz = 20e-3
A_noz = pi*d_noz**2/4

rho_Cu = 0.2 # kg/m
rho_Sh = 0.65 # kg/m
rho_splitter = 2.39 # kg/m2

smallOring = 0.0008 # kg
bigOring = 0.0053 # kg

totalCu = 3.5
totalSh = 0.5

hose_diameter = 0.025
hose_length = 0.75
hose_area = pi*hose_diameter**2/4