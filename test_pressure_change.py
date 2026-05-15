import numpy as np
from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger
from parameters import *


mass_flow_cold = np.array([
    0.492, 0.525, 0.400, 0.467, 0.442,
    0.442, 0.458, 0.417, 0.525, 0.500
])

dp_cold_bar = np.array([
    0.270, 0.227, 0.340, 0.276, 0.293,
    0.303, 0.292, 0.329, 0.221, 0.254
])

mass_flow_cold = mass_flow_cold/1000 * rho_w
dp_cold_pa = dp_cold_bar * 1e5


mass_flow_hot = np.array([
    0.392, 0.316, 0.389, 0.386, 0.378,
    0.378, 0.372, 0.389, 0.316, 0.339
])

dp_hot_bar = np.array([
    0.099, 0.226, 0.105, 0.149, 0.134,
    0.132, 0.149, 0.109, 0.229, 0.095
])

mass_flow_hot = mass_flow_hot/1000 * rho_w
dp_hot_pa = dp_hot_bar * 1e5


# initialise HX
number_of_tubes = [14,12,16,12,12]
number_of_baffles = [5,6,7,8,6]
tube_length = [0.33,0.34,0.27,0.34,0.34]
tube_pitch = [14,14,14,14,14]
is_square_layout = [False,False,False,True,True]
tube_outer_diameter = do
tube_inner_diameter = di
shell_inner_diameter = ds
nozzle_area_shell_side = A_noz
nozzle_area_tube_side = A_noz
fluid_density = rho_w
fluid_viscosity = mu
hose_diameter = hose_diameter
hose_length = hose_length

for i in range(5):
    hx = ShellAndTubeHeatExchanger(number_of_tubes[i],number_of_baffles[i],tube_length[i],tube_pitch[i],is_square_layout[i],tube_outer_diameter,tube_inner_diameter,shell_inner_diameter,
                               nozzle_area_shell_side, nozzle_area_tube_side, fluid_density, fluid_viscosity, hose_diameter, hose_length)


    dp_cold_pred = hx.cold_side_pressure_drop(mass_flow_cold[i])
    dp_hot_pred  = hx.hot_side_pressure_drop(mass_flow_hot[i])

    error_cold = dp_cold_pred - dp_cold_pa[i]
    error_hot  = dp_hot_pred - dp_hot_pa[i]

    rel_error_cold = 100*(error_cold / dp_cold_pa[i])
    rel_error_hot  = 100*(error_hot / dp_hot_pa[i])

    print("Rel error cold:" , rel_error_cold,"%. Rel error hot: ",rel_error_hot,"%.",sep="")
