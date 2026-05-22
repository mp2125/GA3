from parameters import *
import numpy as np

# TODO update these to use hx as function inputs

# shell side
def handoutThermalCoefficient(ReSh, ReTu, baffle_spacing, shape='triangle'):
    # c = 0.2059 if shape == 'square' else 0.1579
    c = 0.0697 if shape == 'square' else 0.0639

    Nui = 0.023  * ReTu**0.8 * Pr**0.3
    Nuo = c * ReSh**0.6 * Pr**0.3 * (ds / baffle_spacing)

    hi = Nui * k_w / di
    ho = Nuo * k_w / do

    Hinv = 1/hi + (di * np.log(do/di)) / (2 * k_tube) + (di/do) / ho
    return 1 / Hinv

def tempIteratorLMTD(length, tubes, mdot1, mdot2, H, T1in=Tcold_in, T2in=Thot_in, shell_passes=1, tube_passes=2):
    Aheat = np.pi * di * length * tubes

    def LMTD(T1in, T1out, T2in, T2out):
        dT1 = T2in - T1out
        dT2 = T2out - T1in
        if abs(dT1 - dT2) < 1e-9:
            return dT1
        return (dT1 - dT2) / np.log(dT1 / dT2)

    def passCorrection(T1in, T1out, T2in, T2out, shell_passes, tube_passes):
        """
        F-factor correction for shell-and-tube heat exchangers.

        tube_passes=1        → pure counter-flow, F=1
        tube_passes even     → 1-2 formula applies (F identical for all even Nt)
        shell_passes>1       → N-shell-in-series combination applied first

        Raises ValueError for odd tube_passes > 1 (non-standard, rarely used).
        """
        if abs(T1out - T1in) < 1e-9:
            return 1.0

        # Single tube pass = pure counter-flow
        if tube_passes == 1:
            return 1.0

        # For even tube passes (2, 4, 6, ...) the F-factor formula
        # is identical to the 1-2 case — tube-pass count beyond the first
        # pair does not change F.
        P1 = (T1out - T1in) / (T2in - T1in)
        R  = (T2in - T2out) / (T1out - T1in)

        # For N shells in series, invert the combination formula to get
        # the per-shell effectiveness P, then apply the single-shell F formula.
        if shell_passes > 1:
            ratio = (1 - P1 * R) / (1 - P1)
            Pdash = ratio ** (1.0 / shell_passes)
            P = (Pdash - 1) / (Pdash - R)
        else:
            P = P1

        if abs(R - 1) < 1e-6:
            numerator   = np.sqrt(2) * P
            denominator = (1 - P) * np.log(
                (2 - P * (2 - np.sqrt(2))) /
                (2 - P * (2 + np.sqrt(2)))
            )
        else:
            numerator   = np.sqrt(R**2 + 1) * np.log((1 - P) / (1 - P * R))
            denominator = (R - 1) * np.log(
                (2 - P * (R + 1 - np.sqrt(R**2 + 1))) /
                (2 - P * (R + 1 + np.sqrt(R**2 + 1)))
            )

        return numerator / denominator

    T1out = 25.0
    T2out = T2in - (mdot1 / mdot2) * (T1out - T1in)
    Qe    = mdot1 * c_p * (T1out - T1in)
    Qlmtd = H * Aheat * passCorrection(T1in, T1out, T2in, T2out, shell_passes, tube_passes) * LMTD(T1in, T1out, T2in, T2out)

    while abs(Qe - Qlmtd) / abs(Qe) > 1e-4:
        T1out = T1in + Qlmtd / (mdot1 * c_p)
        T2out = T2in - (mdot1 / mdot2) * (T1out - T1in)
        Qe    = mdot1 * c_p * (T1out - T1in)
        Qlmtd = H * Aheat * passCorrection(T1in, T1out, T2in, T2out, shell_passes, tube_passes) * LMTD(T1in, T1out, T2in, T2out)

    return T1out, T2out, Qe

def eNTUProcessor(length, tubes, mdot1, mdot2, H, T1_in=Tcold_in, T2_in=Thot_in, config='N-2N', N=1):

    def effectiveness(NTU, C_r, config='1-2', N=1):
        """
        Returns ε given NTU and C_r for the specified configuration.
        N is the number of shell passes (only used for 'N-2N').
        """
        match config:
            case 'counterflow':
                if C_r < 1e-6: return 1 - np.exp(-NTU)                            # condenser/evaporator
                if C_r == 1:   return NTU / (1 + NTU)
                return (1 - np.exp(-NTU * (1 - C_r))) / (1 - C_r * np.exp(-NTU * (1 - C_r)))
            case 'parallelflow':
                return (1 - np.exp(-NTU * (1 + C_r))) / (1 + C_r)
            case '1-2':
                sq = np.sqrt(1 + C_r**2)
                return 2 / (1 + C_r + sq / np.tanh(NTU * sq / 2))
            case 'N-2N':
                e1 = effectiveness(NTU / N, C_r, config='1-2')                     # per-shell ε
                if abs(C_r - 1) < 1e-6: return N * e1 / (1 + (N - 1) * e1)       # limiting form
                ratio = (1 - e1 * C_r) / (1 - e1)
                return (ratio**N - 1) / (ratio**N - C_r) * 0.95
            case 'crossflow-both-unmixed':                                          # approximation; no closed form exists
                return 1 - np.exp((NTU**0.22 / C_r) * (np.exp(-C_r * NTU**0.78) - 1))
            case 'crossflow-Cmax-mixed':
                return (1 / C_r) * (1 - np.exp(-C_r * (1 - np.exp(-NTU))))
            case 'crossflow-Cmin-mixed':
                return 1 - np.exp(-(1 / C_r) * (1 - np.exp(-C_r * NTU)))

    def outlet_temperatures(epsilon, C_cold, C_hot, T1_in, T2_in):
        """Returns (T1_out, T2_out, Q) given ε and inlet conditions."""
        C_min = min(C_cold, C_hot)
        Q     = epsilon * C_min * (T2_in - T1_in)
        return T1_in + Q / C_cold, T2_in - Q / C_hot, Q

    Aheat        = np.pi * di * length * tubes
    C_cold, C_hot = mdot1 * c_p, mdot2 * c_p

    NTU = Aheat * H / min(C_cold, C_hot)
    C_r = min(C_cold, C_hot) / max(C_cold, C_hot)
    eps = effectiveness(NTU, C_r, config, N=N)

    T1_out, T2_out, Q = outlet_temperatures(eps, C_cold, C_hot, T1_in, T2_in)
    return T1_out, T2_out, Q, eps
# eNTU method (https://www.mathworks.com/help/hydro/ref/entuheattransfer.html)

