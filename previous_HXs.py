from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX
from parameters import *
from scipy.optimize import curve_fit
import numpy as np

# Define parameters for previous Heat Exchangers
heat_exchanger_data = [

    # [8,12,0.299,0.0215,False,2,2],
    # num tubes, num baffles, tube length, tube pitch, is square, tube passes, shell passes
    [14,5,0.25,None,False, 2, 2], # 2025A
    [12,6,0.26,None,False, 4, 2], # 2025B
    [16,7,0.195,None,False, 2, 2], # 2025C
    [12,8,0.26,None,True, 2, 2],  # 2025D
    [12,6,0.262,None,True, 2, 2],  # 2025E

    [12,8,0.278,None,False, 2, 1], #2024A
    [12,8,0.260,None,False, 2, 1], #2024B
    [12,8,0.290,None,False, 2, 2], #2024C
    # [14,6,0.250,12e-3,False, 2, 1], #2024D
    # [15,7,0.233,14e-3,False, 3, 1], #2024E

    # [14,10,0.25,14e-3,False, 2, 1], #2023A
    # [16, 6,0.194,14e-3,False, 4, 2], #2023B 
    # [18,12,0.160,14e-3,False, 4, 2], #2023C - Very anomolous result (supposedly no leaks?!?)
    # [10,8,0.246,14e-3,False, 2, 1], #2023D
    # [12,11,0.261,14e-3,False, 2, 1], #2023E

    # [14,13,0.236,14e-3,False, 1, 1], #2022A
    # [16,6,0.219,14e-3,False, 1, 1], #2022B
    # # [20,4,0.174,12e-3,False, 2, 1], #2024C
    # [12,8,0.259,14e-3,False,2,1], # 2022D
    # [24,4,0.145,10e-3,False, 4, 1], #2022E
    # [20,6,0.172,14e-3,False, 4, 2], #2022F

]

experimental_data = [
    # Tcold_in, Tcold_out, dp_cold_measured,
    # Thot_in, Thot_out, dp_hot_measured,
    # Qdot, m_dot_cold, m_dot_hot
    # [20,0,0,60,0,0,0,0,0],
    # --- 2025 data ---
    [20.8, 25.2, 0.270e5, 57.7, 52.1, 0.099e5, 9.02, 0.487, 0.388],  # 25 Group-A
    [22.2, 25.9, 0.227e5, 54.0, 47.7, 0.226e5, 8.14, 0.520, 0.313],  # 25 Group-B
    [23.5, 28.4, 0.340e5, 52.8, 47.4, 0.105e5, 8.40, 0.396, 0.385],  # 25 Group-C
    [24.3, 28.2, 0.276e5, 51.9, 46.4, 0.149e5, 7.95, 0.462, 0.364],  # 25 Group-D
    [25.1, 28.6, 0.293e5, 50.8, 46.1, 0.134e5, 6.88, 0.438, 0.374],  # 25 Group-E

    # --- 2024 data ---
    [19.9, 24.0, 0.161e5, 54.7, 48.8, 0.188e5, 10.29, 0.577, 0.434],  # 24 Group-A
    [21.0, 24.0, 0.201e5, 48.2, 43.3, 0.235e5, 7.29, 0.528, 0.388],   # 24 Group-B
    [20.3, 24.0, 0.195e5, 49.7, 44.1, 0.236e5, 8.92, 0.561, 0.392],   # 24 Group-C
    # [20.5, 24.6, 0.150e5, 51.5, 43.7, 0.138e5, 11.45, 0.577, 0.399],  # 24 Group-D
    # [21.1, 24.6, 0.158e5, 48.2, 42.3, 0.272e5, 8.74, 0.569, 0.371],   # 24 Group-E

    # # --- 2023 data ---
    # [22.1, 25.7, 0.287e5, 48.2, 42.9, 0.186e5, 9.00, 0.594, 0.409], # A
    # [23.1, 27.1, 0.381e5, 49.3, 43.2, 0.264e5, 8.69, 0.495, 0.357], # B
    # [21.6, 24.0, 0.303e5, 49.3, 45.0, 0.252e5, 6.00, 0.561, 0.354], # C
    # [24.2, 28.0, 0.270e5, 55.3, 48.9, 0.219e5, 9.84, 0.602, 0.378], # D
    # [21.0, 25.4, 0.307e5, 49.8, 43.5, 0.180e5, 10.48, 0.569, 0.399], # E

    # # --- 2022 data ---
    # [19.4, 25.3, 0.305e5, 55.4, 50.0, 0.122e5, 10.38, 0.405, 0.479],
    # [20.8, 24.9, 0.160e5, 53.8, 48.8, 0.118e5, 9.62, 0.540, 0.483],
    # # [22.0, 26.3, 0.179e5, 51.5, 46.1, 0.190e5, 9.67, 0.539, 0.430],  hot side leak 
    # [23.0, 27.0, 0.187e5, 51.0, 45.1, 0.221e5, 9.41, 0.528, 0.405],
    # [24.1, 28.0, 0.216e5, 48.8, 43.9, 0.206e5, 8.26, 0.495, 0.413],
    # [25.0, 28.9, 0.276e5, 48.1, 42.4, 0.299e5, 7.82, 0.467, 0.344],
]

hxs = [
    HX(num_tubes, num_baffles, tube_length, tube_pitch=tube_pitch, is_square_layout=is_square, tube_passes=tube_passes, shell_passes=shell_passes)
    for num_tubes, num_baffles, tube_length, tube_pitch, is_square, tube_passes, shell_passes in heat_exchanger_data
]

mass_flows = [data[-2:] for data in experimental_data] # cold, hot
pressure_changes = [[data[2],data[5]] for data in experimental_data] # cold, hot
temperatures = [data[0:2] + data[3:5] for data in experimental_data] # tcoldin tcoldout thotin thotout
heat_transfers = [data[-3] for data in experimental_data]

# this function is for calculating optimal c for Kern model from results. Note that D and E from 2024 are anomolous
# and should be exluded from this. Results from this are then used to manually tune H function in thermal analysis,
# as calling it every time would be unnecessary compute.
def back_calculate_c(hx_list, experimental_data, heat_exchanger_data):
    """
    Back-calculates the shell-side Nusselt constant c from experimental data.
    """
    results = []
    
    for i, (hx, exp, hx_params) in enumerate(zip(hx_list, experimental_data, heat_exchanger_data)):
        # Unpack experimental data
        Tcold_in, Tcold_out, dp_cold, Thot_in, Thot_out, dp_hot, Q_real, mdot1, mdot2 = exp
        n_tubes, n_baffles, length, pitch, is_square, tube_passes, shell_passes = hx_params
        
        Q_real  *= 1000  # kW -> W if necessary
        shape    = 'square' if is_square else 'triangle'

        # Geometry
        Aheat    = np.pi * di * length * n_tubes

        # LMTD (counterflow assumed)
        dT1  = Thot_in  - Tcold_out
        dT2  = Thot_out - Tcold_in
        if abs(dT1 - dT2) < 1e-6:
            LMTD = dT1
        else:
            LMTD = (dT1 - dT2) / np.log(dT1 / dT2)

        # Back-calculate H from real Q and LMTD
        H_real = Q_real / (Aheat * LMTD)

        # Tube-side coefficient (Dittus-Boelter, unchanged)
        velocity_tube = hx.tube_side_velocity(mdot2)
        ReTu          = hx.tube_side_reynolds_number(velocity_tube)
        Nui           = 0.023 * ReTu**0.8 * Pr**0.3
        hi            = Nui * k_w / di

        # Back-calculate ho from H_real
        # 1/H = 1/hi + (di*ln(do/di))/(2*k_tube) + (di/do)/ho
        wall_resistance = (di * np.log(do / di)) / (2 * k_tube)
        ho_real = (di / do) / (1/H_real - 1/hi - wall_resistance)

        # Back-calculate c from ho
        velocity_shell = hx.shell_side_velocity(mdot1)
        ReSh           = hx.shell_side_reynolds_number(velocity_shell)
        Nuo_real       = ho_real * do / k_w
        c_real         = Nuo_real / (ReSh**0.6 * Pr**0.3 * (ds / (length/(n_baffles+1))))

        results.append({
            'case'   : i + 1,
            'shape'  : shape,
            'H_real' : H_real,
            'hi'     : hi,
            'ho_real': ho_real,
            'ReSh'   : ReSh,
            'ReTu'   : ReTu,
            'c_real' : c_real,
            'baffle_spacing': hx.baffle_spacing,
            'Nuo_real'      : ho_real * do / k_w,
        })

        # print(f"Case {i+1:2d} ({shape:8s}): "
        #       f"H_real={H_real:.1f}  hi={hi:.1f}  ho_real={ho_real:.1f}  "
        #       f"ReSh={ReSh:.0f}  c={c_real:.4f}")



    def shell_nusselt(X, c, n):
        ReSh, baffle_spacing = X
        return c * ReSh**0.6 * Pr**0.3 * (ds / baffle_spacing)**n

    valid_tri = [r for r in results if r['shape'] == 'triangle' and r['c_real'] > 0]
    valid_sq  = [r for r in results if r['shape'] == 'square'   and r['c_real'] > 0]

    def shell_nusselt_fixed_n(X, c):
        ReSh, baffle_spacing = X
        n = 1.0  # fixed
        return c * ReSh**0.6 * Pr**0.3 * (ds / baffle_spacing)**n

    def fit_group(cases):
        ReSh_data = np.array([r['ReSh']          for r in cases])
        B_data    = np.array([r['baffle_spacing'] for r in cases])
        Nuo_data  = np.array([r['Nuo_real']       for r in cases])
        popt, pcov = curve_fit(shell_nusselt_fixed_n, (ReSh_data, B_data), Nuo_data, p0=[0.2])
        perr = np.sqrt(np.diag(pcov))
        return popt[0], np.sqrt(pcov[0, 0])

    c_tri, err_tri = fit_group(valid_tri)
    c_sq,  err_sq  = fit_group(valid_sq)

    print(f"Triangular: c={c_tri:.4f}±{err_tri:.4f}")
    print(f"Square:     c={c_sq:.4f}±{err_sq:.4f}")

    return results

if __name__ == '__main__':
    results = back_calculate_c(hxs, experimental_data, heat_exchanger_data)