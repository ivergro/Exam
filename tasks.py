import numpy as np
from Amino import Amino
import plots
import random as r
import copy
from time import time
import os

def initialize(N : int, dim = 2, random = True):
    chain = []
    chain.append(Amino(type = r.randint(1,20), pos=(0,0), index= 0, dim=dim))
    neighbour = np.array([(1,0), (0,1), (-1,0), (0,-1)])
    for i in range(N - 1):
        prev_amino_pos = chain[i].get_pos()
        if random:
            new_pos = tuple(prev_amino_pos + neighbour[r.randint(0,3)])
            chain.append(Amino(type = r.randint(1,20), pos=new_pos, index= i+1,dim=dim ))
            tries = 0
            while chain[i + 1].does_collide(chain):
                new_pos = tuple(prev_amino_pos + neighbour[r.randint(0,3)])
                chain[i + 1].set_pos(new_pos)
                assert tries < 20, "Exceeded 20 tries, something is wrong."
        else:
            new_pos = (i+1,0)
            chain.append(Amino(type = r.randint(1,20), pos=new_pos, index= i+1, dim=dim))
    #Finding nearest neigbhours
    for amino in chain:
        amino.find_nearest_neighbours(chain)
    return np.array(chain)
    


#Points mark interaction energy for given two indexes. E.g point (3,2) marks energy between amino of type 3 and two
def interaction_energy_matrix() -> np.ndarray:
    #Giving lower triangle values
    low_triangle_indices = np.tril_indices(20)
    random_vals = np.random.uniform(-2, -4, len(low_triangle_indices[0]))
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

def center_of_mass(chain : np.ndarray) -> tuple:
    #Assuming mass of amino to be 1 for all types
    N = len(chain)
    M = N*1
    pos = [amino.get_pos() for amino in chain]
    sum_pos = np.sum(pos, axis=0)
    # R = tuple((1/M)*np.add([amino.get_pos() for amino in chain]))
    R = 1/M*sum_pos
    return tuple(R)

def calculate_ROG(chain) -> float:
    #Assuming mass of amino to be 1 for all types
    N = len(chain)
    M   = N*1
    CoM = center_of_mass(chain)
    #Calculate distance from center of mass, r
    r2_sum = sum([np.linalg.norm(np.subtract(amino.get_pos(), CoM))**2 for amino in chain])
    #sum up and then take square root of all squared r's
    #r2 = sum([r**2 for amino in chain])
    return np.sqrt(r2_sum/M), CoM

def end_to_end(chain):
    if np.linalg.norm(np.subtract(chain[-1].get_pos(),chain[0].get_pos())) > (len(chain)-1):
        print("----")
        print(chain[-1].get_pos())
        print(chain[0].get_pos())
        print("....")
    return np.linalg.norm(np.subtract(chain[-1].get_pos(),chain[0].get_pos()))

def Monte_Carlo(chain : np.ndarray, T):
    N = len(chain)
    old_energy = calculate_interaction_energy(chain, energy_matrix)
    flips = 0
    #One sweep
    for i in range(N):
        #Picking randomly
        index = r.randint(0, N -1)
        move_is_possible, new_step, E_flip = Monte_Carlo_step(chain[index], chain)
        flip = False
        if move_is_possible:
            #Favourable to flip if it decrease total energy when flipping
            if E_flip < 0:
                flip = True
            #Not favourable to flip, but still a possibility that it might happen
            elif E_flip >= 0:
                random = r.random()
                boltz_fac = np.exp(-E_flip/(T))
                if random < boltz_fac:
                    flip = True
                else:
                    flip = False
            if flip:
                chain[index].set_pos(new_step)
                chain[index].update_NN(chain)
                flips += 1
    
    
def Monte_Carlo_step(amino : Amino, chain : np.ndarray) :
    available_spots = amino.can_move_to(chain)
    if len(available_spots) > 0:
        new_spot = available_spots[r.randint(0,len(available_spots) - 1)]
        index = amino.get_index()
        old_energy = calculate_interaction_energy(chain, energy_matrix)
        #Performs a deep copy, since the Amino objects where mutable, and got transferred during a regular copy
        new_chain = copy.deepcopy(chain)
        new_chain[index].set_pos(new_spot)
        new_chain[index].update_NN(new_chain)
        new_energy = calculate_interaction_energy(new_chain, energy_matrix)
        delta_E = new_energy - old_energy
        return True, new_spot, delta_E
    return False, None, None

def logger(chain : np.ndarray, old_energy : float, flips : int):
    new_energy = calculate_interaction_energy(chain, energy_matrix)
    delta_E = (new_energy - old_energy)
    end_to_end = np.linalg.norm((chain[-1].get_pos(),chain[0].get_pos()))
    RoG, CoM = calculate_ROG(chain)
    #info_dict = {"delta_E": round(delta_E,3), "flips_during_sweep": flips, "total_E": round(new_energy,3), "RoG":round(RoG,2),"end_to_end": round(end_to_end,2)}
    #logging.info({"delta_E": round(delta_E,3), "flips_during_sweep": flips, "total_E": round(new_energy,3), "RoG":round(RoG,2),"end_to_end": round(end_to_end,2)})



def task_21(): 
    chain= initialize(15, dim = 2)
    total_int_energy = calculate_interaction_energy(chain, energy_matrix)
    print("Energy: ", total_int_energy)
    plots.plot_positions(chain, total_int_energy)
    return chain

def task_25_6(save = True):
    N = 100
    T = 3
    sweeps = 1000
    chain = initialize(N, dim = 2, random=False)
    start_energy = calculate_interaction_energy(chain, energy_matrix)
    plots.plot_positions(chain, start_energy)
    filename =f'loggers/N={N}_T={T}_sweeps={sweeps}/'
    if not os.path.exists(filename):
    # If it doesn't exist, create it
        os.makedirs(filename)
    
    #Saving important parameters to a list
    energies = np.zeros(sweeps + 1)
    end_to_ends = np.zeros(sweeps + 1)
    RoGs    = np.zeros(sweeps+1)

    #Adding original values
    energies[0] = start_energy
    end_to_ends[0] = end_to_end(chain)
    RoGs[0] = calculate_ROG(chain)[0]

    for x in range(1, sweeps + 1):
        Monte_Carlo(chain, T)
        new_E = calculate_interaction_energy(chain, energy_matrix)
        #Adding new values after a sweep
        energies[x] = new_E
        end_to_ends[x] = end_to_end(chain)
        RoGs[x] = calculate_ROG(chain)[0]
        if (x%(sweeps//3) == 0)and save:
            #plots.plot_positions(chain, new_E)
            np.save(filename + f"step={x}_positions", chain)
    #Saving for later use
    if save:
        np.save(filename + "energies", energies)
        np.save(filename + "end_to_ends", end_to_ends)
        np.save(filename + "RoGs", RoGs)
    plots.plot_RoG_ete_energy(RoGs, end_to_ends, energies)
    plots.plot_positions(chain, calculate_interaction_energy(chain, energy_matrix))
    
def task_27(): 
    #Runs for wanted parameters
    Ns = [15, 50, 100] #Kan endres
    Ts = [0.1, 0.5, 1, 2, 3, 4, 5, 6, 7, 8 , 9, 10]
    sweeps = 2000
    #a) plot E og RoG mot antall sweeps, for ulike T. Hvor lang tid tar det å nå likevekt for ulik T
    #b) 
    #c) Hva er kritisk temp, altså T hvor faseoverganger skjer? Er det likt for de ulike N'ene
    for N in Ns:
        chain = initialize(N, dim = 2, random=False)
        start_energy = calculate_interaction_energy(chain, energy_matrix)

        for T in Ts:
            filename =f'loggers/Task_27/N={N}_T={T}_sweeps={sweeps}/'
            if not os.path.exists(filename):
            # If it doesn't exist, create it
                os.makedirs(filename)
            
            #Saving important parameters to a list
            energies = np.zeros(sweeps + 1)
            end_to_ends = np.zeros(sweeps + 1)
            RoGs    = np.zeros(sweeps+1)

            #Adding original values
            energies[0] = start_energy
            end_to_ends[0] = end_to_end(chain)
            RoGs[0] = calculate_ROG(chain)[0]

            for x in range(1, sweeps + 1):
                Monte_Carlo(chain, T)
                new_E = calculate_interaction_energy(chain, energy_matrix)
                #Adding new values after a sweep
                energies[x] = new_E
                end_to_ends[x] = end_to_end(chain)
                RoGs[x] = calculate_ROG(chain)[0]
                if (x == sweeps):
                    #plots.plot_positions(chain, new_E)
                    np.save(filename + f"step={x}_positions", chain)
            #Saving for later use
            np.save(filename + "energies", energies)
            np.save(filename + "end_to_ends", end_to_ends)
            np.save(filename + "RoGs", RoGs)
        
#Defined here so it remains constant during sim
energy_matrix = interaction_energy_matrix()

t_0 = time()

task_27()

print("Time of sim: ", time() - t_0)



        