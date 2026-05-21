from optimiser import fullSolver
from previous_HXs import hxs, temperatures

letters = ['A','B','C','D','E','A','B','C','D','E']
for n in range(0,10):
    print(letters[n], fullSolver(hxs[n],temperatures[n][0],temperatures[n][2])[1])