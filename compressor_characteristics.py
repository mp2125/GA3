import numpy as np
from scipy.interpolate import interp1d

Q_cold = np.array([0.6580, 0.6290, 0.5830, 0.5380, 0.4670,
                   0.3920, 0.3210, 0.2790, 0.2210, 0.0])*1e-3

dp_cold = np.array([0.1584, 0.1958, 0.2493, 0.3127, 0.3723,
                    0.4436, 0.4950, 0.5318, 0.5739, 0.7077])*1e5

Q_hot = np.array([0.4360, 0.3870, 0.3520, 0.3110, 0.2600,
                  0.2290, 0.1670, 0.1180, 0.0690, 0.0010])*1e-3

dp_hot = np.array([0.0932, 0.1688, 0.2209, 0.2871, 0.3554,
                   0.4041, 0.4853, 0.5260, 0.5665, 0.6239])*1e5


def make_interp(Q, dp):
    idx = np.argsort(Q)
    return interp1d(Q[idx], dp[idx], kind='linear', fill_value="extrapolate")

cold_curve = make_interp(Q_cold, dp_cold)
hot_curve = make_interp(Q_hot, dp_hot)