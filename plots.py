import matplotlib.pyplot as plt
from Amino import Amino
import matplotlib.cm as cm
import logging
import numpy as np

#Endre colorbar verdier til -2kb til -4kb
def plot_energy_matrix(energy_matrix):
    fig, ax = plt.subplots()

    im = ax.imshow(energy_matrix, origin="lower",extent=(0,len(energy_matrix), 0, len(energy_matrix)))
    cmap = plt.cm.plasma
    fig.colorbar(im,cmap=cmap)
    
    plt.show()

def plot_positions(chain, total_int_energy):
    fig, ax = plt.subplots()
    positions = []
    colors    = []
    cmap = cm.get_cmap('tab20')
    for amino in chain:
        positions.append(amino.get_pos())
        colors.append(cmap(amino.get_type())) #Asigning each type a specific color to be used in plot, to destinguish between different amino acid types
    
    #Extracting each position value
    x_values, y_values = zip(*positions)
    ax.plot(x_values, y_values, color="k")
    ax.scatter(x_values, y_values, s=200, c = colors)
    # ax.set_ylim(0, gridsize)
    # ax.set_xlim(0,gridsize)
    ax.set_aspect("equal")
    ax.set_adjustable("datalim")
    ax.annotate("Interaction energy: {%.3f}Kb" % (total_int_energy),xy = (0.0,1.05),fontsize = 10, xycoords="axes fraction") #{%.3d}" % total_int_energy,xy

    plt.show()

def plot_multiple_positions(chain_list, total_int_energy_list):
    fig, axes = plt.subplots(len(chain_list))
    for i in range(len(chain_list)):
        positions = []
        colors    = []
        cmap = cm.get_cmap('tab20')
        for amino in chain_list[i]:
            positions.append(amino.get_pos())
            colors.append(cmap(amino.get_type())) #Asigning each type a specific color to be used in plot, to destinguish between different amino acid types
        
        #Extracting each position value
        x_values, y_values = zip(*positions)
        axes[i].plot(x_values, y_values, color="k")
        axes[i].scatter(x_values, y_values, s=100, c = colors)
        # ax.set_ylim(0, gridsize)
        # ax.set_xlim(0,gridsize)
        axes[i].set_aspect("equal")
        axes[i].set_adjustable("datalim")
        axes[i].annotate("Interaction E: %.3f$k_B$" % (total_int_energy_list[i]),xy = (0.0,1.05),fontsize = 10, xycoords="axes fraction") #{%.3d}" % total_int_energy,xy
        axes[i].set_xticks([])
        axes[i].set_title(f"sweeps={10**i}")
        #axes[i].set_yticks([])
    
    #fig.savefig("Results/task216/3_config")
    plt.show()


def plot_E(energies):
    fig, ax = plt.subplots()
    x = np.linspace(0,len(energies), num = len(energies))
    ax.plot(x, energies)
    ax.set_xlabel("MC sweeps")
    ax.set_ylabel("Energy [$k_B$]")
    
    plt.show()

def plot_end_to_end(filename):
    fig, ax = plt.subplots()
    end_to_ends = np.load(filename)
    x = np.linspace(0,len(end_to_ends), num = len(end_to_ends))
    ax.plot(x, end_to_ends)
    ax.set_xlabel("MC sweeps")
    ax.set_ylabel("End to end distance")
    
    plt.show()

def plot_RoG(filename):
    fig, ax = plt.subplots()
    RoGs = np.load(filename)
    x = np.linspace(0,len(RoGs), num = len(RoGs))
    ax.plot(x, RoGs)
    ax.set_xlabel("MC sweeps")
    ax.set_ylabel("Unit distance")
    ax.set_title("Radius of gyration")
    
    plt.show()

def plot_RoG_ete(RoGs, end_to_ends):
    fig, ax = plt.subplots()
    x = np.linspace(0,len(RoGs), num = len(RoGs))

    ax.plot(x, RoGs, label = "Radius of gyration")
    ax.plot(x, end_to_ends, label = "End to end")
    ax.set_xlabel("MC sweeps")
    ax.set_ylabel("Unit distance")
    ax.set_title(f"N={N}, T={T}")
    ax.legend()

    plt.show()

def plot_RoG_ete_energy(RoGs, end_to_ends,energies):
    fig, ax1 = plt.subplots()
    x = np.linspace(0,len(RoGs), num = len(RoGs))

    f_size = 13
    ax1.plot(x, RoGs, label = "Radius of gyration", color = "tab:green")
    ax1.plot(x, end_to_ends, label = "End to end")
    ax1.set_xlabel("MC sweeps", fontsize=f_size)
    ax1.set_ylabel("Unit distance", fontsize=f_size)
    ax1.set_title(f"N={N}, T={T}")

    c = "tab:red"
    ax2 = ax1.twinx()
    ax2.plot(x, energies, label = "Energy", color=c)
    ax2.set_ylabel("Interaction energy [$k_B$]", color=c, fontsize=f_size)
    ax2.tick_params(axis='y', labelcolor=c)
    fig.legend()

    #fig.savefig("Results/task216/energy_ete_RoG_plot")
    plt.show()

def phase_diagrams():
    #Plot average E as a func of selected temps [0.1,0.5,1,2,...,10] and average RoG
    #Average should go from chosen starting point. starting when...
    Ts = [0.1,0.5,1,2,3,4,5,6,7,8,9,10]
    sweeps = 1000
    N = 15
    E_T = []
    RoG_T = []
    for T in Ts:
        filename = f"loggers/Task_27/sweeps=1000/N=15_T={T}_sweeps=1000/"
        energies = np.load(filename + "energies.npy")
        RoGs     = np.load(filename + "RoGs.npy")
        #Choosing starting point
        starting_point = 0
        N              = len(energies[starting_point:])
        E_T.append(sum(energies[starting_point:])/N)
        RoG_T.append(sum(RoGs[starting_point:])/N)
    
    fig, ax1 = plt.subplots()

    ax1.plot(Ts, E_T, label="E(T)", color=E_color)
    ax2 = ax1.twinx()
    ax2.plot(Ts, RoG_T, label="RoG(T)", color=RoG_color)
    ax1.set_xlabel("T")
    ax1.set_ylabel("Unit distance", color=RoG_color)
    ax2.set_ylabel("Energy [$k_B$]", color = E_color)
    ax2.tick_params(axis='y', labelcolor=E_color)
    ax1.tick_params(axis='y', labelcolor=RoG_color)
    fig.legend()

    plt.show()



#----Colors----
RoG_color = "tab:green"
E_color   = "tab:red"
ete_color = "tab:blue"

N = 15
T = 1
# sweeps = 1000
# filename = f"loggers/N={N}_T={T}_sweeps={sweeps}/"

# RoGs = np.load(filename + "RoGs.npy")
# end_to_ends = np.load(filename + "end_to_ends.npy")
# energies = np.load(filename + "energies.npy")
# positions = [np.load(filename + "step=1_positions.npy", allow_pickle=True), np.load(filename + "step=10_positions.npy", allow_pickle=True), np.load(filename + "step=100_positions.npy", allow_pickle=True)]
# plot_positions(position_list[2], energies[100])
# plot_RoG_ete(RoGs, end_to_ends)
# plot_E(energies)
# plot_RoG_ete_energy(RoGs, end_to_ends, energies)
# plot_multiple_positions(positions, [energies[1], energies[10], energies[100]])

phase_diagrams()