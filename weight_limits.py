from previous_HXs import hxs
from parameters import *
from numpy import pi, arange
from ShellAndTubeHeatExchanger import ShellAndTubeHeatExchanger as HX

def getWeight(hx):

    wholeLength = hx.tube_length + 0.075

    pipeWeight = hx.number_of_tubes * hx.tube_length * rho_Cu
    shellWeight = wholeLength * rho_Sh
    splitterWeight = (hx.shell_passes - 1) * ds * wholeLength * rho_splitter
    baffleWeight = hx.number_of_baffles * (ds**2)*pi/4 * 0.8 * rho_splitter

    print(splitterWeight)
    return pipeWeight + shellWeight + splitterWeight + baffleWeight + 0.1

# realWeights = [1.124, 1.095, 1.029, 1.081, 1.074,
#                1.110, 1.127, 1.364, 1.145, 1.162,
#                1.101, 1.096, 0.964, 0.965, 1.094,
#                ]

# weights = []
# for hx in hxs:
#     weights.append(getWeight(hx))

# x = arange(0,len(weights))
# import matplotlib.pyplot as plt

# plt.plot(x,weights,'o-',label='Calculated')
# plt.plot(x,realWeights,'s-',label='Real')
# plt.legend()
# plt.show()