from optimiser import fullSolver
from previous_HXs import hxs, temperatures, heat_transfers

letters = ['A','B','C','D','E','A','B','C']
for n in range(0,8):
    print(letters[n], f'{fullSolver(hxs[n],temperatures[n][0],temperatures[n][2])[1]/1000:.2f}', heat_transfers[n])