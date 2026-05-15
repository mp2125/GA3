# ---------------------------------------------------------------------------
# ε-NTU heat exchanger equations
# Configurations: 'counterflow', 'parallelflow', '1-2', 'N-2N', 'crossflow'
# Stream 1 = cold, Stream 2 = hot
# C_r  = C_min / C_max  where C = m_dot * c_p
# NTU  = H * A / C_min
# ---------------------------------------------------------------------------

import numpy as np

def effectiveness(NTU, C_r, config='counterflow', N=1):
    """
    Returns ε given NTU and C_r for the specified configuration.
    N is the number of shell passes (only used for 'N-2N').
    """
    match config:
        case 'counterflow':
            if C_r < 1e-6:  return 1 - np.exp(-NTU)                           # condenser/evaporator
            if C_r == 1:    return NTU / (1 + NTU)
            return (1 - np.exp(-NTU * (1 - C_r))) / (1 - C_r * np.exp(-NTU * (1 - C_r)))
 
        case 'parallelflow':
            return (1 - np.exp(-NTU * (1 + C_r))) / (1 + C_r)
 
        case '1-2':
            sq = np.sqrt(1 + C_r**2)
            return 2 / (1 + C_r + sq / np.tanh(NTU * sq / 2))
 
        case 'N-2N':
            e1   = effectiveness(NTU / N, C_r, config='1-2')                   # per-shell ε
            if abs(C_r - 1) < 1e-6: return N * e1 / (1 + (N - 1) * e1)       # limiting form
            ratio = (1 - e1 * C_r) / (1 - e1)
            return (ratio**N - 1) / (ratio**N - C_r)
 
        case 'crossflow-both-unmixed':                                          # approximation; no closed form exists
            return 1 - np.exp((NTU**0.22 / C_r) * (np.exp(-C_r * NTU**0.78) - 1))
 
        case 'crossflow-Cmax-mixed':
            return (1 / C_r) * (1 - np.exp(-C_r * (1 - np.exp(-NTU))))
 
        case 'crossflow-Cmin-mixed':
            return 1 - np.exp(-(1 / C_r) * (1 - np.exp(-C_r * NTU)))
 
 
def NTU_required(epsilon, C_r, config='counterflow'):
    """
    Returns NTU needed to achieve target ε.
    Only defined analytically for counterflow, parallelflow, and 1-2.
    """
    match config:
        case 'counterflow':
            if C_r < 1e-6:  return -np.log(1 - epsilon)
            if C_r == 1:    return epsilon / (1 - epsilon)
            return np.log((1 - epsilon * C_r) / (1 - epsilon)) / (1 - C_r)
 
        case 'parallelflow':
            return -np.log(1 - epsilon * (1 + C_r)) / (1 + C_r)
 
        case '1-2':
            sq = np.sqrt(1 + C_r**2)
            A  = (2 / epsilon - 1 - C_r) / sq
            return (2 / sq) * np.arctanh(1 / A)
 
 
def outlet_temperatures(epsilon, C_cold, C_hot, T1_in, T2_in):
    """Returns (T1_out, T2_out) given ε and inlet conditions."""
    C_min = min(C_cold, C_hot)
    Q     = epsilon * C_min * (T2_in - T1_in)
    return T1_in + Q / C_cold, T2_in - Q / C_hot
 
 
# ---------------------------------------------------------------------------
# Example
# ---------------------------------------------------------------------------
 
if __name__ == "__main__":
    C_cold, C_hot  = 1.0 * 4200, 0.8 * 4200    # W/K
    T1_in, T2_in   = 20.0, 80.0                 # °C
    NTU            = (500 * 2.0) / min(C_cold, C_hot)
    C_r            = min(C_cold, C_hot) / max(C_cold, C_hot)
 
    configs = ['counterflow', 'parallelflow', '1-2',
               'crossflow-both-unmixed', 'crossflow-Cmax-mixed', 'crossflow-Cmin-mixed']
 
    print(f"NTU = {NTU:.3f}, C_r = {C_r:.3f}\n")
    for cfg in configs:
        eps            = effectiveness(NTU, C_r, cfg)
        T1_out, T2_out = outlet_temperatures(eps, C_cold, C_hot, T1_in, T2_in)
        print(f"{cfg:<30}  ε={eps:.4f}  T1_out={T1_out:.1f}°C  T2_out={T2_out:.1f}°C")
 
    print("\nN-2N shells in series:")
    for N in [1, 2, 4]:
        eps            = effectiveness(NTU, C_r, 'N-2N', N=N)
        T1_out, T2_out = outlet_temperatures(eps, C_cold, C_hot, T1_in, T2_in)
        print(f"  N={N}  ε={eps:.4f}  T1_out={T1_out:.1f}°C  T2_out={T2_out:.1f}°C")