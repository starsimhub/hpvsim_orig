'''
Defines the stochastic cell model for HPV (Markov chain)
'''

import numpy as np
import random
import matplotlib.pyplot as plt


#Event-driven

def draw_tau():
    tau = random.expovariate(1.0)
    return tau

def draw_event_class(infected_list):
    symmetric_split_bb_rate = 0
    symmetric_split_pp_rate = 0
    asymmetric_split_bp_rate = 0

    for cell in infected_list: # need to check into if this is right
        asymmetric_split_bp_rate += cell.event_rate
        symmetric_split_bb_rate += cell.event_rate
        symmetric_split_pp_rate += cell.event_rate

    asym_start = 0
    asym_end = asymmetric_split_bp_rate
    sym_bb_end = asym_end + symmetric_split_bb_rate
    sym_pp_end = sym_bb_end + symmetric_split_pp_rate

    random_draw = random.uniform(asym_start, sym_pp_end)

    if random_draw < asym_end:
        return 3 # asymmetric split
    elif (random_draw >= asym_end) & (random_draw < sym_bb_end):
        return 2 # symmetric split (BB)
    elif (random_draw >= sym_bb_end) & (random_draw < sym_pp_end):
        return 1 # symmetric split (PP)


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

    def initialize(self):
        cell_zero_idx = random.randint(0,len(size)-1)
        self.V_B.append(B_Cell(cell_zero_idx,1,5)) #change split rate eventually
        self.V_B_ind.append(cell_zero_idx)
        self.indices.remove(cell_zero_idx)

    def run_sim(self): #Function to run simulation of events
        self.initialize()
        while self.current_time < self.time:

            tau = draw_tau()
            event_class = draw_event_class(self.V_B)
            index_1 = np.random.choice(self.indices)
            index_2 = np.random.choice(self.indices)
            while index_2 == index_1:
                index_2 = np.random.choice(self.indices)

            if event_class == 1:
                # Make two new parabasal cells
                symm_pp_event = draw_event(np.max(self.Lambda), self.V_B)
                new_p1, new_p2 = symm_pp_event.split(1, index_1, index_2, self.current_time)

                #Bookkeeping and updating the lists of infected cells
                self.V_P.extend((new_p1, new_p2))
                self.V_P_ind.extend((new_p1.index, new_p2.index))
                self.indices.remove(new_p1.index)
                self.indices.remove(new_p2.index)
                self.V_B.remove(symm_pp_event)
                self.indices.append(symm_pp_event.index)

            if event_class == 2:
                # Make two new basal cells
                symm_bb_event = draw_event(np.max(self.Lambda), self.V_B)
                new_b1, current = symm_bb_event.split(2, index_1, symm_bb_event.index, self.current_time)

                # Bookkeeping and updating the lists of infected cells
                self.V_B.append(new_b1)
                self.V_B_ind.append(new_b1.index)
                self.indices.remove(new_b1.index)

            if event_class == 3:
                # Make a basal and a parabasal cell
                asym_pb_event = draw_event(np.max(self.Lambda), self.V_B)
                new_p, current = asym_pb_event.split(3, index_1, asym_pb_event.index, self.current_time)

                # Bookkeeping and updating the lists of infected cells
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

class B_Cell: # CHANGE THIS TO BASAL_CELL
    # index: index of the basal cell
    # state: 0-normal, 1-infected, 2-transformed, 3-integrated
    # event_rate: rate of event, the type of split that occurs
    # vl: viral load measurement


    def __init__(self, index, state, split_rate):
        self.index = index
        self.state = state
        self.event_rate = split_rate

    def split(self, event_type, index_1, index_2, t):
        # Make the extra cell, default parabasal cell which is infected
        if event_type == 1: # two parabasal cells
            split_1 = P_Cell(index_1, 1000, 3+t)
            split_2 = P_Cell(index_2, 1000, 3+t)

        if event_type == 2: # two new basal cells
            split_1 = B_Cell(index_1, 1, 5)
            split_2 = self

        if event_type == 3: # one basal, one parabasal
            split_1 = P_Cell(index_1, 1000, 3+t)
            split_2 = self

        return split_1, split_2

    def transform(self):
        self.state = 2

    def integrate(self):
        self.state = 3

    def display_info(self):
        print('Basal cell index: ', self.index,' state: ', self.state, ' split_rate: ', self.event_rate)


class P_Cell:
    # index: index of the parabasal cell
    # vl: viral load let out
    # shed_time: time when the p-cell dies and sheds its viral load
    # TAKE THE BASAL CELL OBJECT INTO THE INFORMATION


    def __init__(self, index, vl, shed_time):
        self.index = index
        self.vl = vl
        self.shed_time = shed_time

    def shed(self):
        return self.vl

    def display_info(self):
        print('Parabasal cell index: ', self.index, ' viral load: ', self.vl, ' shed time: ', self.shed_time)

