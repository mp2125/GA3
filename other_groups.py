from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX
from optimiser import fullSolver
from weight_limits import getWeight
from pprint import pprint
from hydraulic_analysis import solve_mass_flows


group_d1 = HX(10, 11, 0.243, tube_pitch=0.02, is_square_layout=False, tube_passes=2, shell_passes=2)
group_d2 = HX(8, 12, 0.299, tube_pitch=0.020, is_square_layout=False, tube_passes=2, shell_passes=2)
group_b = HX(16, 7, 0.174, tube_pitch=None, is_square_layout=False, tube_passes=4, shell_passes=2)
group_a = HX(12, 14, 0.236, tube_pitch=0.012, is_square_layout=False, tube_passes=2, shell_passes=1)
group_c = HX(12, 7, 0.230, tube_pitch=0.012, is_square_layout=False, tube_passes=2, shell_passes=2)
group_e = HX(16, 8, 0.185, is_square_layout=False, tube_passes=2, shell_passes=1)

other_hxs = [
    group_a, group_b, group_c, group_d2, group_e
]

predicted_transfers = [fullSolver(group)[1] for group in other_hxs]
weights = [getWeight(group) for group in other_hxs]

if __name__ == '__main__':
    print("Group, Number of Tubes, Number of Baffles, Tube Length (mm), Tube Pitch, Hot Passes, Cold Passes, Heat Transfer (W)")
    for hx, prediction, group in zip(other_hxs, predicted_transfers, "ABCDE"):
        pprint(", ".join(map(str,[
            group,
            hx.number_of_tubes,
            hx.number_of_baffles,
            hx.tube_length*1e3,
            round(float(hx.tube_pitch*1e3),0),
            hx.tube_passes,
            hx.shell_passes,
            f'{prediction:.0f}',
            ])))

# mflow = solve_mass_flows(group_d2)
# print(mflow)
# print(group_d2.hot_side_pressure_drop(mflow[1]))