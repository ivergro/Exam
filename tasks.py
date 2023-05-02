import numpy as np
from Amino import Amino

def task_21_1():
    N = 15 #number of monomers
    chain = np.zeros(15)
    for i in range(N):
        chain[i] = Amino(i+1, (i,0), dim=2)
    
