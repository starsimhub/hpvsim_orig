import numpy as np
import matplotlib.pyplot as plt
import cell

def run():
    print('running')
    simulate()

    plt.plot()
    plt.legend(loc='upper right')
    plt.xlabel('Duration')
    plt.ylabel('Probability')
    plt.semilogy()
    plt.show()



def simulate():


    sim = cell.MarkovSim.__init__(time = t, lambda, theta)
