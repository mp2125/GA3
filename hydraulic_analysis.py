import parameters
from compressor_characteristics import cold_curve, hot_curve, min_cold_mass_flow, max_cold_mass_flow, min_hot_mass_flow, max_hot_mass_flow
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
import numpy as np
from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX



# f = lambda Re: (1.82*log10(Re) - 1.64)**(-2)

# This code is all redundant now
"""def get_Re_shell(mdot1,L,Y,Nb,square):
    B = L/(Nb+1)
    A_shell = ds/Y * (Y-do)*B
    v_shell = mdot1/(rho_w*A_shell)
    A_pipe = pi*ds**2/4
    # dprime_shell = ds*(A_shell/A_pipe)
    # Re_shell = rho_w*v_shell*dprime_shell/mu
    De = 4*(Y**2 - np.pi*do**2/4) / (np.pi*do) if square else 4*(3**0.5 * Y**2/4 - np.pi*do**2/8) / (np.pi*do/2)
    Re_shell = rho_w*v_shell*De/mu

    return Re_shell

def get_del_p1(mdot1,L, Y, N, Nb, square):
    B = L/(Nb+1)
    A_shell = ds/Y * (Y-do)*B
    v_shell = mdot1/(rho_w*A_shell)
    Re_shell = get_Re_shell(mdot1,L, Y, Nb, square)
    pitch_const = 0.34 if square else 0.2
    # del_p_shell = 4*pitch_const*Re_shell**(-0.15)*N*rho_w*v_shell**2
    Gs = mdot1 / A_shell
    del_p_shell = pitch_const * Re_shell**(-0.15) * (Nb+1) * Gs**2 / (2*rho_w)
    v_noz1 = mdot1/(rho_w*A_noz)
    del_p_noz1 = 2*0.5*rho_w*v_noz1**2
    del_p1 = del_p_shell+del_p_noz1
    return del_p1


def get_Re_tube(mdot2,L, N, Nb):
    B = L/(Nb+1)
    mdot_tube = mdot2 / N
    v_tube = mdot_tube / (rho_w*Ai_water)
    Re_tube = rho_w*v_tube*di/mu

    return Re_tube

def get_del_p2(mdot2,L, N, Nb):
    Re_tube = get_Re_tube(mdot2,L,N,Nb)
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

    return del_p2"""

def cold_residual(cold_mass_flow, hx: HX):
    model_dp = hx.cold_side_pressure_drop(cold_mass_flow)
    compressor_dp = cold_curve(cold_mass_flow)

    return model_dp - compressor_dp

def hot_residual(hot_mass_flow, hx: HX):
    model_dp = hx.hot_side_pressure_drop(hot_mass_flow)
    compressor_dp = hot_curve(hot_mass_flow)
    return model_dp - compressor_dp


def solve_mass_flows(hx):
    # returns the mass flows for the operating point (intersection with compressor curve)
    cold_solution = root_scalar(
        cold_residual,
        args=(hx,),
        bracket=[min_cold_mass_flow,max_cold_mass_flow],
        method='brentq'
    )

    hot_solution = root_scalar(
        hot_residual,
        args=(hx,),
        bracket=[min_hot_mass_flow,max_hot_mass_flow],
        method='brentq'
    )

    return cold_solution.root, hot_solution.root


if __name__ == "__main__":
    cold_mass_flow_guess = 0.5
    hot_mass_flow_guess = 0.45
    number_of_tubes = 13
    tube_length = 0.35
    number_of_baffles = 9
    pitch = 0.014
    square = True

    test_hx = HX(
        number_of_tubes,
        number_of_baffles,
        tube_length,
        pitch,
        square
    )

    # Mass flow range
    cold_mass_flows = np.linspace(min_cold_mass_flow, max_cold_mass_flow, 500)
    hot_mass_flows = np.linspace(min_hot_mass_flow, max_hot_mass_flow, 500)


    # HX characteristics
    cold_hx_dp = []
    hot_hx_dp = []

    # Compressor characteristics
    cold_comp_dp = []
    hot_comp_dp = []

    for cold_mdot, hot_mdot in zip(cold_mass_flows,hot_mass_flows):
        # HX pressure drops
        cold_hx_dp.append(test_hx.cold_side_pressure_drop(cold_mdot))
        hot_hx_dp.append(test_hx.hot_side_pressure_drop(hot_mdot))

        # Compressor curves
        cold_comp_dp.append(cold_curve(cold_mdot))
        hot_comp_dp.append(hot_curve(hot_mdot))

    # Convert to arrays
    cold_hx_dp = np.array(cold_hx_dp)
    hot_hx_dp = np.array(hot_hx_dp)

    cold_comp_dp = np.array(cold_comp_dp)
    hot_comp_dp = np.array(hot_comp_dp)

    # Plot
    plt.figure(figsize=(10, 6))

    # Cold side
    plt.plot(
        cold_mass_flows,
        cold_hx_dp,
        label="Cold HX Characteristic",
        linewidth=2
    )

    plt.plot(
        cold_mass_flows,
        cold_comp_dp,
        label="Cold Compressor Characteristic",
        linestyle="--",
        linewidth=2
    )

    # Hot side
    plt.plot(
        hot_mass_flows,
        hot_hx_dp,
        label="Hot HX Characteristic",
        linewidth=2
    )

    plt.plot(
        hot_mass_flows,
        hot_comp_dp,
        label="Hot Compressor Characteristic",
        linestyle="--",
        linewidth=2
    )

    # Operating points
    cold_op_mass_flow, hot_op_mass_flow = solve_mass_flows(test_hx)

    plt.scatter(
        cold_op_mass_flow,
        test_hx.cold_side_pressure_drop(cold_op_mass_flow),
        s=80,
        label="Cold Operating Point"
    )

    plt.scatter(
        hot_op_mass_flow,
        test_hx.hot_side_pressure_drop(hot_op_mass_flow),
        s=80,
        label="Hot Operating Point"
    )

    plt.xlabel("Mass Flow Rate [kg/s]")
    plt.ylabel("Pressure Drop / Pressure Rise [Pa]")
    plt.title("HX and Compressor Characteristics")
    plt.grid(True)
    plt.legend()
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.show()

