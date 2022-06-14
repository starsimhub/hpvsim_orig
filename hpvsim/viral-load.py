'''
Define the viral load distributions here
'''

import numpy as np
import matplotlib.pyplot as plt


def vl_functions(start, peak_vl):

    # This is going to be defined by the amount of shedding cells at each time point.


    x = np.linspace(0, 2*start)
    y = -(x - start)**2 + peak_vl

    plt.plot(x, y)
    plt.xlabel("Duration")
    plt.ylabel("Viral load (copies/cell)")
    plt.show()

# def assign_VL_curve(person, time):
#     if infected:
#
#
#     return vl



class Lesion:

    '''Inititalizing the lesion to have 1 infected cell'''
    def __init__(self, num_cells):

        self.num_cells = 1000

        self.basal_cells = np.zeros((0.6* self.num_cells)) # 6% of cells are stem cells
        self.parabasal_cells = np.zeros((2, self.num_cells - self.basal_cells)) # The rest are cells that cannot replicate
        # first row is for infected, second row is for shedding When they die??

        # Infect one basal cell
        self.basal_cells[0] = 1

        self. normal_cells =(self.num_cells - ( np.sum(self.basal_cells) + np.sum(self.parabasal_cells[0,:]) ) )/self.num_cells
        self.infected_cells = (1 - self.normal_cells)/self.num_cells # proportion of cells infected
        self.shedding_cells = np.sum(self.parabasal_cells[1,:])/self.num_cells # subset of infected cells and parabasal cells
        self.transformed_cells = self.normal_cells/self.num_cells

    def getNumCells(self):
        return self.num_cells

    def getBasalCells(self):
        return len(self.basal_cells)/self.num_cells

    def getParabasalCells(self):
        return len(self.parabasal_cells)/self.num_cells

    def getInfectedCells(self):
        return (1 - self.normal_cells)/self.num_cells

    def getSheddingCells(self):
        return np.sum(self.parabasal_cells[1,:])/self.num_cells

    def getTransformedCells(self):
        return self.normal_cells/self.num_cells

    def infect(self): ## This is the important aspect
        i_cells = np.sum(self.basal_cells) #number of infected basal cells

        temp = np.random.randint(0,3) #not necessarily random PLACEHOLDER
        if temp == 0:
            #asymmetric spread

        elif temp == 1:
        # symmetric spread

        else:
            self.parabasal_cells




    def progression(self):








##

