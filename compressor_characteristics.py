import numpy as np
from scipy.interpolate import interp1d
from parameters import rho_w, hose_area

m_dot_cold = np.array([0.6580, 0.6290, 0.5830, 0.5380, 0.4670,
                   0.3920, 0.3210, 0.2790, 0.2210, 0.001])*1e-3*rho_w

dp_cold = np.array([0.1584, 0.1958, 0.2493, 0.3127, 0.3723,
                    0.4436, 0.4950, 0.5318, 0.5739, 0.7077])*1e5*0.75

m_dot_hot = np.array([0.4360, 0.3870, 0.3520, 0.3110, 0.2600,
                  0.2290, 0.1670, 0.1180, 0.0690, 0.0010])*1e-3*rho_w

dp_hot = np.array([0.0932, 0.1688, 0.2209, 0.2871, 0.3554,
                   0.4041, 0.4853, 0.5260, 0.5665, 0.6239])*1e5*0.8


def make_interp(mdot, dp):
    idx = np.argsort(mdot)
    return interp1d(mdot[idx], dp[idx], kind='linear', fill_value="extrapolate")

cold_curve = make_interp(m_dot_cold, dp_cold)
hot_curve = make_interp(m_dot_hot, dp_hot)

min_cold_mass_flow = min(m_dot_cold)
max_cold_mass_flow = max(m_dot_cold)
min_hot_mass_flow = min(m_dot_hot)
max_hot_mass_flow = max(m_dot_hot)

k_hose_cold = dp_cold[0] / (0.5*rho_w*(max_cold_mass_flow/(rho_w*hose_area))**2) *0.5
k_hose_hot = dp_hot[0] / (0.5*rho_w*(max_hot_mass_flow/(rho_w*hose_area))**2) *0.5

# print(2*k_hose_hot)
