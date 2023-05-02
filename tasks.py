import numpy as np
from Amino import Amino

def task_21_1():
    N = 15 #number of monomers
    chain = np.zeros(15)
    for i in range(N):
        chain[i] = Amino(i+1, (i,0), i, dim=2)
    
    #Finding nearest neigbhours
    for i in range(N):
        chain[i].find_nearest_neighbour(chain)
    
def task_21_2():
    pass
    
