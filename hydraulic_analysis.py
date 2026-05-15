from parameters import *
from math import log10
from loss_coefficients import kc_turb, ke_turb, kc_lam, ke_lam
from compressor_characteristics import cold_curve, hot_curve
from scipy.optimize import root_scalar


f = lambda Re: (1.82*log10(Re) - 1.64)**(-2)  # TODO make use of Danny's Moody

# def get_pressure_changes(mdot1, mdot2, N, L, Nb, Y, square):
#     # determine mass flow rates for energy balance
#     B = L/(Nb+1)
#     mdot_tube = mdot2 / N
#     v_tube = mdot_tube / (rho_w*Ai_water)
#     Re_tube = rho_w*v_tube*di/mu
#     v_noz2 = mdot2/(rho_w*A_noz)

#     # determine pressure drops in hot stream
#     del_p_tube = 0.5*rho_w*v_tube**2*(f(Re_tube)*L/di)


#     sigma = N*(di/ds)**2
#     if Re_tube > 3000:
#         kc = kc_turb(Re_tube, sigma)
#         ke = ke_turb(Re_tube, sigma)
#     else:
#         kc = kc_lam(L/di, Re_tube, sigma)

#     del_p_end = 0.5*rho_w*v_tube**2*(kc+ke)
#     del_p_noz2 = 2*0.5*rho_w*v_noz2**2
#     del_p2 = del_p_tube+del_p_end+del_p_noz2

#     # determine the pressure drop in the cold stream

#     A_shell = ds/Y * (Y-do)*B
#     v_shell = mdot1/(rho_w*A_shell)
#     A_pipe = pi*ds**2/4
#     dprime_shell = ds*(A_shell/A_pipe)
#     Re_shell = rho_w*v_shell*dprime_shell/mu
#     pitch_const = 0.34 if square else 0.2
#     del_p_shell = 4*pitch_const*Re_shell**(-0.15)*N*rho_w*v_shell**2
#     v_noz1 = mdot1/(rho_w*A_noz)
#     del_p_noz1 = 2*0.5*rho_w*v_noz1**2
#     del_p1 = del_p_shell+del_p_noz1

#     return del_p1, del_p2


def get_Re_shell(mdot1,Y,Nb):
    B = L/(Nb+1)
    A_shell = ds/Y * (Y-do)*B
    v_shell = mdot1/(rho_w*A_shell)
    A_pipe = pi*ds**2/4
    dprime_shell = ds*(A_shell/A_pipe)
    Re_shell = rho_w*v_shell*dprime_shell/mu

    return Re_shell

def get_del_p1(mdot1, Y, N, Nb, square):
    B = L/(Nb+1)
    A_shell = ds/Y * (Y-do)*B
    v_shell = mdot1/(rho_w*A_shell)
    Re_shell = get_Re_shell(mdot1, Y, Nb)
    pitch_const = 0.34 if square else 0.2
    del_p_shell = 4*pitch_const*Re_shell**(-0.15)*N*rho_w*v_shell**2
    v_noz1 = mdot1/(rho_w*A_noz)
    del_p_noz1 = 2*0.5*rho_w*v_noz1**2
    del_p1 = del_p_shell+del_p_noz1

    return del_p1


def get_Re_tube(mdot2, N, Nb):
    B = L/(Nb+1)
    mdot_tube = mdot2 / N
    v_tube = mdot_tube / (rho_w*Ai_water)
    Re_tube = rho_w*v_tube*di/mu

    return Re_tube

def get_del_p2(mdot2, N, Nb):
    Re_tube = get_Re_tube(mdot2,N,Nb)
    mdot_tube = mdot2 / N
    v_tube = mdot_tube / (rho_w*Ai_water)
    v_noz2 = mdot2/(rho_w*A_noz)

    # determine pressure drops in hot stream
    del_p_tube = 0.5*rho_w*v_tube**2*(f(Re_tube)*L/di)

    sigma = N*(di/ds)**2
    if Re_tube > 3000:
        kc = kc_turb(Re_tube, sigma)
        ke = ke_turb(Re_tube, sigma)
    else:
        kc = kc_lam(L/di, Re_tube, sigma)
        ke = ke_lam(L/di, Re_tube, sigma)

    del_p_end = 0.5*rho_w*v_tube**2*(kc+ke)
    del_p_noz2 = 2*0.5*rho_w*v_noz2**2
    del_p2 = del_p_tube+del_p_end+del_p_noz2

    return del_p2

def cold_residual(m_dot1, Y, N, Nb, square):
    Q = m_dot1 / rho_w
    system_dp = get_del_p1(m_dot1, Y, N, Nb, square)
    curve_dp = cold_curve(Q)
    return system_dp - curve_dp

def hot_residual(m_dot2, N, Nb):
    Q = m_dot2 / rho_w
    system_dp = get_del_p2(m_dot2, N, Nb)
    curve_dp = hot_curve(Q)
    return system_dp - curve_dp


def solve_mass_flows(Y, N, Nb, square):
    cold_solution = root_scalar(
        cold_residual,
        args=(Y,N,Nb,square),
        bracket=[1e-2, 1],
        method='brentq'
    )

    hot_solution = root_scalar(
        hot_residual,
        args=(N,Nb),
        bracket=[1e-2, 1],
        method='brentq'
    )

    return cold_solution.root, hot_solution.root

    print("Cold mass flow rate:", cold_solution.root)
    print("Hot mass flow rate:", hot_solution.root)

if __name__ == "__main__":
    mdot1 = 0.5
    mdot2 = 0.45
    N = 13
    L=0.35
    Nb = 9
    Y = 0.014
    square = True

    cold_mdot, hot_mdot = solve_mass_flows(Y, N, Nb, square)
    print(cold_mdot, hot_mdot)
    print(get_Re_shell(cold_mdot,Y,Nb))
    print(get_Re_tube(hot_mdot,N,Nb))
