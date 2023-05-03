import matplotlib.pyplot as plt
from Amino import Amino
from scipy.constants import Boltzmann
import matplotlib.cm as cm

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


