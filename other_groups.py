from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX
from optimiser import fullSolver
from weight_limits import getWeight


group_d1 = HX(10, 11, 0.243, tube_pitch=0.02, is_square_layout=False, tube_passes=2, shell_passes=2)
group_d2 = HX(8, 12, 0.299, tube_pitch=0.0215, is_square_layout=False, tube_passes=2, shell_passes=2)
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
    for prediction in predicted_transfers:
        print(f'{prediction:.0f} W')

    for weight in weights:
        print(weight)


