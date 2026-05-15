from parameters import *
from math import log10
from loss_coefficients import kc_turb, ke_turb, kc_lam

f = lambda Re: (1.82*log10(Re) - 1.64)**(-2)

def get_mass_and_pressure(mdot1, mdot2, N, L, Nb, Y, square):
    # determine mass flow rates for energy balance
    B = L/(Nb+1)
    mdot_tube = mdot2 / N
    v_tube = mdot_tube / (rho_w*A_tube)
    Re_tube = rho_w*v_tube*di/mu
    v_noz2 = mdot2/(rho_w*A_noz)

    # determine pressure drops in hot stream
    del_p_tube = 0.5*rho_w*v_tube**2*(f(Re_tube)*L/di)


    sigma = N*(di/ds)**2
    if Re_tube > 3000:
        kc = kc_turb(Re_tube, sigma)
        ke = ke_turb(Re_tube, sigma)
    else:
        kc = kc_lam(L/di, Re_tube, sigma)

    del_p_end = 0.5*rho_w*v_tube**2*(kc+ke)
    del_p_noz2 = 2*0.5*rho_w*v_noz2**2
    del_p2 = del_p_tube+del_p_end+del_p_noz2

    # determine the pressure drop in the cold stream

    A_shell = ds/Y * (Y-do)*B
    v_shell = mdot1/(rho_w*A_shell)
    A_pipe = pi*ds**2/4
    dprime_shell = ds*(A_shell/A_pipe)
    Re_shell = rho_w*v_shell*dprime_shell/mu
    pitch_const = 0.34 if square else 0.2
    del_p_shell = 4*pitch_const*Re_shell**(-0.15)*N*rho_w*v_shell**2
    v_noz1 = mdot1/(rho_w*A_noz)
    del_p_noz1 = 2*0.5*rho_w*v_noz1**2
    del_p1 = del_p_shell+del_p_noz1

    print(del_p1, del_p2)

mdot1 = 0.5
mdot2 = 0.45
N = 13
L=0.35
Nb = 9
Y = 0.014
square = True

get_mass_and_pressure(mdot1,mdot2,N,L,Nb,Y,square)
