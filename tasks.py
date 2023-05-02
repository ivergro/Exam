import numpy as np
from Amino import Amino
import plots
from scipy.constants import Boltzmann
import random as r


def initialize(N : int, dim = 2):
    chain = []
    chain.append(Amino(type = r.randint(1,20), pos=(0,0), index= 0, dim=dim))
    neighbour = np.array([(1,0), (0,1), (-1,0), (0,-1)])
    for i in range(N - 1):
        prev_amino_pos = chain[i].get_pos()
        new_pos = tuple(prev_amino_pos + neighbour[r.randint(0,3)])
        chain.append(Amino(type = r.randint(1,20), pos=new_pos, index= i+1,dim=dim ))
        tries = 0
        while chain[i + 1].does_collide(chain):
            new_pos = tuple(prev_amino_pos + neighbour[r.randint(0,3)])
            chain[i + 1].set_pos(new_pos)
            assert tries < 20, "Exceeded 20 tries, something is wrong."
    #Finding nearest neigbhours
    for amino in chain:
        amino.find_nearest_neighbours(chain)
    return np.array(chain)
    


#Points mark interaction energy for given two indexes. E.g point (3,2) marks energy between amino of type 3 and two
def interaction_energy_matrix() -> np.ndarray:
    #Giving lower triangle values
    low_triangle_indices = np.tril_indices(20)
    random_vals = np.random.uniform(-2*Boltzmann, -4*Boltzmann, len(low_triangle_indices[0]))
    random_matrix = np.zeros((20,20))
    random_matrix[low_triangle_indices] = random_vals

    #copy lower triangle to upper triangle, except diagonal
    random_matrix += random_matrix.T - np.diag(random_matrix.diagonal()) 

    return random_matrix

def calculate_interaction_energy(chain, energy_matrix):
    total_J = 0
    for amino in chain:
        current_type = amino.get_type()
        for NN in amino.get_NN():
            neighbour_type = NN.get_type()
            total_J += energy_matrix[current_type-1][neighbour_type-1]
    total_J /= 2 #Removing double counting
    return total_J

def Monte_Carlo(chain : np.ndarray):
    N = len(chain)
    for i in range(N):
        index = r.randint(0, N -1)
        Monte_Carlo_step(chain[index], chain)

def Monte_Carlo_step(amino : Amino, chain : np.ndarray) :
    available_spots = amino.can_move_to(chain)
    if len(available_spots) > 0:
        new_spot = available_spots[r.randint(0,len(available_spots) - 1)]
        index = amino.get_index()
        old_energy = calculate_interaction_energy(chain, energy_matrix)
        new_chain = chain.copy()
        new_chain[index].set_pos(new_spot)
        new_chain[index].update_NN(new_chain)
        new_energy = calculate_interaction_energy(new_chain, energy_matrix)
        print(new_energy - old_energy)

    return


def task_21(): 
    chain= initialize(15, 2)
    total_int_energy = calculate_interaction_energy(chain, energy_matrix)
    print("Energy: ", total_int_energy)
    plots.plot_positions(chain, total_int_energy)
    return chain

#Defined here so it remains constant during sim
energy_matrix = interaction_energy_matrix()

Monte_Carlo(task_21())



# chain, gridsize = initialize(25,20)
# plots.plot_positions(chain, gridsize)