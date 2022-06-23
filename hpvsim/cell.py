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
    tau = random.expovariate(1.0)               # draws a random time interval
    return tau

def draw_event_class(V_B, V_P):
    basal_split_bb_rate = 0                     # draws the event class possibility
    basal_split_pp_rate = 0
    basal_split_bp_rate = 0
    pbasal_split_pp_rate = 0
    infect_rate = 0
    diff_rate = 0

    for cell in V_B: # TODO: Make sure that this is correct and do not need to have different event rates
        basal_split_bp_rate += cell.event_rate
        basal_split_bb_rate += cell.event_rate
        basal_split_pp_rate += cell.event_rate
        infect_rate += cell.event_rate

    for cell in V_P:
        pbasal_split_pp_rate += cell.event_rate
        diff_rate += cell.event_rate


    basal_bp_start = 0
    basal_bp_end = basal_split_bp_rate
    basal_bb_end = basal_bp_end + basal_split_bb_rate
    basal_pp_end = basal_bb_end + basal_split_pp_rate
    pbasal_pp_end = basal_pp_end + pbasal_split_pp_rate
    infect_end = pbasal_pp_end + infect_rate
    diff_end = infect_end + diff_rate


    random_draw = random.uniform(basal_bp_start, diff_end)

    if random_draw < basal_bp_end:
        return 6                                                        # asymmetric split (BP) from basal
    elif (random_draw >= basal_bp_end) & (random_draw < basal_bb_end):
        return 5                                                        # symmetric split (BB) from basal
    elif (random_draw >= basal_bb_end) & (random_draw < basal_pp_end):
        return 4                                                        # symmetric split (PP) from basal
    elif (random_draw >= basal_pp_end) & (random_draw < pbasal_pp_end):
        return 3                                                        # symmetric split (PP) from parabasal
    elif (random_draw >= pbasal_pp_end) & (random_draw < infect_end):
        return 2                                                        # infection event
    elif (random_draw >= infect_end) & (random_draw < diff_end):
        return 1                                                        # differentiation event


def draw_event(max_rate, event_list):                       # out of the event list, random draw occurs of if that
    accepted = False                                        # event will occur.
    random_event = None
    while not accepted:
        random_event = random.choice(event_list)
        accept_rate = random_event.event_rate / max_rate
        random_draw = random.uniform(0,1)
        if random_draw < accept_rate:
            accepted = True

    return random_event


###### Class declarations



class MarkovSim:                            # Simulation to let the event-driven stochastic process occur
    def __init__(self, time, Lambda, theta, size):
        self.indices = list(range(1, 100))  # number of cells in the system
        self.time = time                    # total time of simulation
        self.current_time = 0               # current time in simulation
        self.Lambda = Lambda                # placeholders for rates
        self.theta = theta
        self.V_B = []                       # vec of basal cells (objects)
        self.V_B_ind = []                   # vec of basal cell indices
        self.V_P_non_diff = []              # vec of parabasal cells that are not differentiating
        self.V_P_non_diff_ind = []          # vec of indices of above
        self.V_P_diff = []                  # vec of parabasal cells that have potential to divide
        self.V_P_diff_ind = []              # indices of potential dividing parabasal cells
        self.vl_time = np.zeros(time)       # measured in weeks?
        self.infected_cells = []            # vector for infected cells
        self.infected_cells_indices = []    # indices of infected cells
        self.num_infected_t = []            # counter for the cells infected at each time point

    def initialize(self):
        cell_zero_idx = random.randint(0, len(size)-1)
        self.V_B.append(Basal_Cell(cell_zero_idx, 1, 5) # TODO: change the rates, currently placeholders
        self.V_B_ind.append(cell_zero_idx)
        self.indices.remove(cell_zero_idx)

    def run_sim(self): #Function to run simulation of events
        self.initialize()
        while self.current_time < self.time:

            tau = draw_tau()
            event_class = draw_event_class(self.V_B, self.V_P_non_diff)
            index_1 = np.random.choice(self.indices)
            index_2 = np.random.choice(self.indices)
            while index_2 == index_1:
                index_2 = np.random.choice(self.indices)

            if event_class == 1:
                # Differentiating a parabasal cell
                diff_event = draw_event(np.max(self.Lambda), self.V_P_non_diff)
                diff_event.differentiate()
                self.V_P_diff.append(diff_event)
                self.V_P_diff_ind.append(diff_event.index)
                self.V_P_non_diff.remove(diff_event)
                self.V_P_non_diff_ind.remove(diff_event.index)


            if event_class == 2:
                # Infect a basal cell
                infect_event = draw_event(np.max(self.Lambda), self.V_B)
                infection = Infected_Cell('16', 500, infect_event)
                self.infected_cells.append(infection)
                self.infected_cells_indices.append(infection.cell.index)

            if event_class == 3:
                # Make two new parabasal cells from parabasal cell
                pbasal_pp_event = draw_event(np.max(self.Lambda), self.V_P_non_diff)
                new_p1, current = pbasal_pp_event.split(1, index_1, pbasal_pp_event.index, self.current_time)

                #Bookkeeping and updating the lists of cells
                self.V_P_non_diff.append(new_p1)
                self.V_P_non_diff_ind.append(new_p1.index)
                self.indices.remove(new_p1.index)

                #Bookkeeping the infection
                if isinstance(pbasal_pp_event, Infected_Cell):
                    self.infected_cells.append((Infected_Cell(pbasal_pp_event.genotype, pbasal_pp_event.vl, new_p1),Infected_Cell(pbasal_pp_event.genotype, pbasal_pp_event.vl, new_p2)))
                    self.infected_cells_indices((new_p1.index,new_p2.index))

            if event_class == 4:
                # Make two new parabasal cells from basal cell
                basal_pp_event = draw_event(np.max(self.Lambda), self.V_B)
                new_p1, new_p2 = basal_pp_event.split(2, index_1, index_2, self.current_time)

                #Bookkeeping and updating the lists of cells
                self.V_P_non_diff.extend((new_p1, new_p2))
                self.V_P_non_diff_ind.extend((new_p1.index, new_p2.index))
                self.indices.remove(new_p1.index)
                self.indices.remove(new_p2.index)
                self.V_B.remove(basal_pp_event)
                self.indices.append(basal_pp_event.index)

                # Bookkeeping the infection
                if isinstance(new_p1, Infected_Cell):
                    self.infected_cells.append(new_b1)
                    self.infected_cells_indices(new_b1.index)

            if event_class == 5:
                # Make two new basal cells from basal cell
                basal_bb_event = draw_event(np.max(self.Lambda), self.V_B)
                new_b1, current = basal_bb_event.split(2, index_1, basal_bb_event.index, self.current_time)

                # Bookkeeping and updating the lists of cells
                self.V_B.append(new_b1)
                self.V_B_ind.append(new_b1.index)
                self.indices.remove(new_b1.index)

                # Bookkeeping the infection
                if isinstance(new_b1, Infected_Cell):
                    self.infected_cells.append(new_b1)
                    self.infected_cells_indices(new_b1.index)

            if event_class == 6:
                # Make a basal and a parabasal cell from basal cell
                basal_pb_event = draw_event(np.max(self.Lambda), self.V_B)
                new_p, current = basal_pb_event.split(3, index_1, basal_pb_event.index, self.current_time)

                # Bookkeeping and updating the lists of  cells
                self.V_P.append(new_p)
                self.V_P_ind.append(new_p.index)
                self.indices.remove(new_p.index)

                # Bookkeeping the infection
                if isinstance(new_p, Infected_Cell):
                    self.infected_cells.append(new_p)
                    self.infected_cells_indices(new_p.index)
            # Need to account for cell death and shedding here, only pull from V_P_diff vector


            self.current_time += tau


    #Set up a plotting tool for the VL
    #def plot_vl(self):


class Basal_Cell(Cell):
    # index: index of the basal cell
    # cell: cell that made it (and all its attributes)

    def __init__(self, index, split_rate, death_rate, cell):
        super().__init__(index, split_rate, death_rate)
        self.cell = cell

    def split(self, event_type, index_1, index_2, t, infect):
        # Make the extra cell, default parabasal cell which is infected
        if event_type == 1: # two parabasal cells
            split_1 = Parabasal_Cell(index_1, self.split_rate, self.death_rate, self)
            split_2 = Parabasal_Cell(index_2, self.split_rate, self.death_rate, self)



        if event_type == 2: # two new basal cells
            split_1 = Basal_Cell(index_1, self.split_rate, self.death_rate, self)
            split_2 = self

        if event_type == 3: # one basal, one parabasal
            split_1 = Parabasal_Cell(index_1, self.split_rate, self.death_rate, self)
            split_2 = self

        return split_1, split_2




class Parabasal_Cell(Cell):
    # index: index of the parabasal cell
    # state: 0-dividing, 1-differentiating
    # cell: cell that made it (and all it attributes)


    def __init__(self, index, split_rate, death_rate, cell):
        super().__init__(index, split_rate, death_rate)
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
    # index: index of the cell
    # split_rate: rate of division
    # death_rate: rate of death

    def __init__(self, index, split_rate, death_rate):
        self.index = index
        self.split_rate = split_rate
        self.death_rate = death_rate

    def update_split(self, new_rate):
        self.split_rate = new_rate

    def update_death(self, new_death):
        self.death_rate = new_death

    def die(self):
        self.index = 99

class Infected_Cell(Basal_Cell):
    # cell: the cell that has just become infected
    # genotype: the genotype of the virus in the system
    # vl: viral load of the cell

    def __init__(self, genotype, vl, cell):
        super().__init__(cell.index, cell.split_rate, cell.death_rate, cell)
        self.genotype = genotype
        self.vl = vl

    def update_vl(self):
        self.vl = 1000

    def shed(self):
        return self.vl

class Infected_Cell(Parabasal_Cell):
    # cells: associated cells that the virus is inhabiting
    # cell: the cell that has just become infected
    # genotype: the genotype of the virus in the system
    # vl: viral load of the cell

    def __init__(self, genotype, vl, cell):
        super().__init__(cell.index, cell.split_rate, cell.death_rate, cell.state, cell)
        self.genotype = genotype
        self.vl = vl


    def update_vl(self):
        self.vl = 1000

    def shed(self):
        return self.vl

