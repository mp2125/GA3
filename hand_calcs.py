from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX
from previous_HXs import hxs, mass_flows

hx = hxs[3]
m_cold, m_hot = mass_flows[3]

# baffle spacing = 0.0288
print(hx.baffle_spacing)

# tube area = 2.83e-5
print(hx.tube_side_flow_area)

vtube = 2.16
print(hx.tube_side_velocity(m_hot))

# Re tube = 19700
print(hx.tube_side_reynolds_number(vtube))

hx.tube_side_pressure_drop(m_hot)

# vnoz_hot = 1.17
print(hx.nozzle_velocity_hot)

# del p tube = 