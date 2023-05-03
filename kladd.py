import numpy as np



mot_tup = np.subtract((1,1,1), (1,1,0))

print(np.linalg.norm(mot_tup))

import matplotlib.cm as mcolors

# get the default colormap in matplotlib
cmap = mcolors.get_cmap('tab20')

# create a list of 20 distinct colors from the colormap
colors = [cmap(i) for i in range(20)]

# print the list of colors
print(colors)