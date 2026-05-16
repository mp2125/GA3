from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX

# Define parameters for previous Heat Exchangers
heat_exchanger_data = [
    # num tubes, num baffles, tube length, tube pitch, is square, tube passes
    [14,5,0.33,14,False], # 2025A
    [12,6,0.34,14,False], # 2025B
    [16,7,0.27,14,False], # 2025C
    [12,8,0.34,14,True],  # 2025D
    [12,6,0.34,14,True],  # 2025E

    [12,8,0.278,14,False], #2024A
    [12,8,0.260,14,False], #2024B
    [12,8,0.290,14,False], #2024C
    [14,6,0.250,12,False], #2024D
    [15,7,0.233,14,False], #2024E
]

experimental_data = [
    # Tcold_in, Tcold_out, dp_cold_measured,
    # Thot_in, Thot_out, dp_hot_measured,
    # Qdot, m_dot_cold, m_dot_hot

    # --- 2025 data ---
    [20.8, 25.2, 0.260, 57.7, 52.1, 0.090, 9.02, 0.487, 0.388],  # 25 Group-A
    [22.2, 25.9, 0.219, 54.0, 47.7, 0.215, 8.14, 0.520, 0.313],  # 25 Group-B
    [23.5, 28.4, 0.325, 52.8, 47.4, 0.096, 8.40, 0.396, 0.385],  # 25 Group-C
    [24.3, 28.2, 0.265, 51.9, 46.4, 0.140, 7.95, 0.462, 0.364],  # 25 Group-D
    [25.1, 28.6, 0.280, 50.8, 46.1, 0.125, 6.88, 0.438, 0.374],  # 25 Group-E

    # --- 2024 data ---
    [19.9, 24.0, 0.120, 54.7, 48.8, 0.138, 10.29, 0.577, 0.434],  # 24 Group-A
    [21.0, 24.0, 0.156, 48.2, 43.3, 0.183, 7.29, 0.528, 0.388],   # 24 Group-B
    [20.3, 24.0, 0.152, 49.7, 44.1, 0.184, 8.92, 0.561, 0.392],   # 24 Group-C
    [20.5, 24.6, 0.109, 51.5, 43.7, 0.187, 11.45, 0.577, 0.399],  # 24 Group-D
    [21.1, 24.6, 0.116, 48.2, 42.3, 0.220, 8.74, 0.569, 0.371],   # 24 Group-E
]

hxs = [
    HX(num_tubes, num_baffles, tube_length, tube_pitch, is_square)
    for num_tubes, num_baffles, tube_length, tube_pitch, is_square in heat_exchanger_data
]
