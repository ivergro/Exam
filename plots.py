import matplotlib.pyplot as plt
from Amino import Amino

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
    for amino in chain:
        positions.append(amino.get_pos())
    
    #Extracting each position value
    x_values, y_values = zip(*positions)
    ax.plot(x_values, y_values, color="k")
    ax.scatter(x_values, y_values, s=200)
    # ax.set_ylim(0, gridsize)
    # ax.set_xlim(0,gridsize)
    ax.set_aspect("equal")
    ax.annotate(f"Interaction energy: {total_int_energy}",xy = (0.0,1.05),fontsize = 10, xycoords="axes fraction") #{%.3d}" % total_int_energy,xy

    plt.show()


