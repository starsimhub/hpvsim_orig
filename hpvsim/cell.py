'''
Defines the stochastic cell model for HPV (Markov chain)
'''

import numpy as np
import random
import matplotlib.pyplot as plt


'''
TO DO:
Allow for random viral load to be floating around and be in the system
Normal cell division should happen, then allow there to be an infected section 
'''

#Event-driven

def draw_tau():
    tau = random.expovariate(1.0)
    return tau

def draw_event_class(V_B, V_P):
    basal_split_bb_rate = 0
    basal_split_pp_rate = 0
    basal_split_bp_rate = 0
    pbasal_split_pp_rate = 0
    infect_rate = 0

    for cell in V_B: # need to check into if this is right
        basal_split_bp_rate += cell.event_rate
        basal_split_bb_rate += cell.event_rate
        basal_split_pp_rate += cell.event_rate
        infect_rate += cell.event_rate

    for cell in V_P:
        pbasal_split_pp_rate += cell.event_rate


    basal_bp_start = 0
    basal_bp_end = basal_split_bp_rate
    basal_bb_end = basal_bp_end + basal_split_bb_rate
    basal_pp_end = basal_bb_end + basal_split_pp_rate
    pbasal_pp_end = basal_pp_end + pbasal_split_pp_rate
    infect_end = pbasal_pp_end + infect_rate


    random_draw = random.uniform(basal_bp_start, infect_end)

    if random_draw < basal_bp_end:
        return 5 # asymmetric split (BP) from basal
    elif (random_draw >= basal_bp_end) & (random_draw < basal_bb_end):
        return 4 # symmetric split (BB) from basal
    elif (random_draw >= basal_bb_end) & (random_draw < basal_pp_end):
        return 3 # symmetric split (PP) from basal
    elif (random_draw >= basal_pp_end) & (random_draw < pbasal_pp_end):
        return 2 # symmetric split (PP) from parabasal
    elif (random_draw >= pbasal_pp_end) & (random_draw < infect_end):
        return 1 # infection event


def draw_event(max_rate, event_list):
    accepted = False
    random_event = None
    while not accepted:
        random_event = random.choice(event_list)
        accept_rate = random_event.event_rate / max_rate
        random_draw = random.uniform(0,1)
        if random_draw < accept_rate:
            accepted = True

    return random_event


class MarkovSim:
    def __init__(self, time, Lambda, theta, size):
        self.indices = list(range(1, 100)) # number of cells in the system
        self.time = time
        self.current_time = 0
        self.Lambda = Lambda #place holders for rates
        self.theta = theta
        self.V_B = [] #vec of basal cells (objects)
        self.V_B_ind = [] #vec of basal cell indices
        self.V_P = []
        self.V_P_ind = []
        self.vl_time = np.zeros(time) # measured in weeks?
        self.infected_cells = []

    def initialize(self):
        cell_zero_idx = random.randint(0, len(size)-1)
        self.V_B.append(Basal_Cell(cell_zero_idx, 1, 5)) #change split rate eventually
        self.V_B_ind.append(cell_zero_idx)
        self.indices.remove(cell_zero_idx)

    def run_sim(self): #Function to run simulation of events
        self.initialize()
        while self.current_time < self.time:

            tau = draw_tau()
            event_class = draw_event_class(self.V_B, self.V_P)
            index_1 = np.random.choice(self.indices)
            index_2 = np.random.choice(self.indices)
            while index_2 == index_1:
                index_2 = np.random.choice(self.indices)

            if event_class == 1:
                # Infect a basal cell
                infect_event = draw_event(np.max, self.V_B)
                self.infected_cells.append()
            if event_class == 2:
                # Make two new parabasal cells from parabasal cell
                pbasal_pp_event = draw_event(np.max(self.Lambda), self.V_B)
                new_p1, new_p2 = pbasal_pp_event.split(1, index_1, index_2, self.current_time)

            if event_class == 3:
                # Make two new parabasal cells from basal cell
                basal_pp_event = draw_event(np.max(self.Lambda), self.V_B)
                new_p1, new_p2 = basal_pp_event.split(2, index_1, index_2, self.current_time)

                #Bookkeeping and updating the lists of cells
                self.V_P.extend((new_p1, new_p2))
                self.V_P_ind.extend((new_p1.index, new_p2.index))
                self.indices.remove(new_p1.index)
                self.indices.remove(new_p2.index)
                self.V_B.remove(basal_pp_event)
                self.indices.append(basal_pp_event.index)

            if event_class == 4:
                # Make two new basal cells from basal cell
                basal_bb_event = draw_event(np.max(self.Lambda), self.V_B)
                new_b1, current = basal_bb_event.split(2, index_1, basal_bb_event.index, self.current_time)

                # Bookkeeping and updating the lists of infected cells
                self.V_B.append(new_b1)
                self.V_B_ind.append(new_b1.index)
                self.indices.remove(new_b1.index)

            if event_class == 5:
                # Make a basal and a parabasal cell from basal cell
                basal_pb_event = draw_event(np.max(self.Lambda), self.V_B)
                new_p, current = basal_pb_event.split(3, index_1, basal_pb_event.index, self.current_time)

                # Bookkeeping and updating the lists of  cells
                self.V_P.append(new_p)
                self.V_P_ind.append(new_p.index)
                self.indices.remove(new_p.index)

            # Need to account for cell death and shedding here NEED TO FIGURE OUT WANING VIRAL LOAD
            # NEED TO CHECK OTHER PARABASAL CELLS SHEDDING HERE AS WELL. PUT IN A FOR LOOP FOR LEN OF V_P
            if V_P[0].shed_time <= self.current_time:
                self.vl_time[self.current_time] = V_P[0].shed()
                self.vl_time[self.current_time + 1] = V_P[0].shed()
                self.vl_time[self.current_time + 2] = V_P[0].shed()/2
                self.V_P.remove(V_P[0])
                self.V_P_ind.remove(V_P[0].index)
                self.indices.append(V_P[0].index)

            self.current_time += tau


    #Set up a plotting tool for the VL
    #def plot_vl(self):


class Basal_Cell:
    # index: index of the basal cell
    # cell: cell that made it (and all its attributes)
    # event_rate: rate of event, the type of split that occurs


    def __init__(self, index, cell, split_rate):
        self.index = index
        self.cell = cell
        self.event_rate = split_rate

    def split(self, event_type, index_1, index_2, t):
        # Make the extra cell, default parabasal cell which is infected
        if event_type == 1: # two parabasal cells
            split_1 = Parabasal_Cell(index_1, 0, self, 1000, 3+t)
            split_2 = Parabasal_Cell(index_2, 0, self, 1000, 3+t)

        if event_type == 2: # two new basal cells
            split_1 = Basal_Cell(index_1, 1, 5)
            split_2 = self

        if event_type == 3: # one basal, one parabasal
            split_1 = Parabasal_Cell(index_1, 0, self, 1000, 3+t)
            split_2 = self

        return split_1, split_2




class Parabasal_Cell:
    # index: index of the parabasal cell
    # state: 0-dividing, 1-differentiating
    # cell: cell that made it (and all it attributes)
    # event_rate: rate of and event, the a split occurs



    def __init__(self, index, state, cell):
        self.index = index
        self.state = 0
        self.cell = cell

    def differentiate(self):
        self.state = 1


    def split(self, index_1, index_2, t):
        # Make the extra cell, default parabasal cell which is infected
        split_1 = Parabasal_Cell(index_1, 0, self)
        split_2 = Parabasal_Cell(index_2, 0, self)

        return split_1, split_2


class Cell:
    # Start cell for the process

class Virus:
    # total: number of cells that are infected
    # cells: associated cells that the virus is inhabiting
    # genotype: the genotype of the virus in the system
    # vl: viral load of the cell

    def __init__(self, genotype, vl):
        self.cell = []
        self.genotype = genotype
        self.vl = vl

    def update_vl(self):
        self.vl = 1000

    def shed(self):
        return self.vl