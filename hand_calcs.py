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

Re_tube = 19700
print(hx.tube_side_reynolds_number(vtube))

hx.tube_side_pressure_drop(m_hot)

# vnoz_hot = 1.17
print(hx.nozzle_velocity_hot)

# friction_factor tube
print(hx.friction_factor(Re_tube, 0))

# pipe friction loss = 5240
print(hx.pipe_friction_loss)

# sigma = 0.0525
print(hx.sigma)

#Kc = 0.45, ke = 0.8
print(hx.get_entrance_exit_coefficients(Re_tube))

# entrance_exit losses = 6230
print(hx.entrance_exit_loss)

print(hx.K_hose_cold, hx.K_hose_hot)
# dp hose = 2.65kpa each side
print(hx.hose_pressure_drop(m_hot, hx.K_hose_hot))

print(hx.hot_side_pressure_drop(m_hot))


# happy with hot side physicalness, now lets tune it