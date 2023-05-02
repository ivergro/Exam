import numpy as np

#Kan evt sjekke om avstanden til punktene er lik en lattice konstanten, heller enn om de er 1 unna hverandre
#Kan gjøre sånn at når en amino blir en annens NN, får den andre også den sin NN
#Må kanskje ikke ha NN avhening av pos, kan også bare ha en liste
#NN må skille på chain og bare nearest neigbhour
class Amino:
    """
    Amino class, or monomer class. 
    Represents a monomer in a protein structure

    Variables
    ---------

    type : int
        Number ranging from 1-20, depending on which amino acid type they are
    
    pos : tuple
        Tuple indicating grid position. Could take 2D or 3D format.
    
    index : int
        Chain index of monomer

    dim : int
        Dimensions of grid. Used to determine size of Nearest neigbour array
    
    NN : np.ndarray
        Array of all nearest neigbours. 0 marks unoccupied space. Ranges from 4-6 in size depending og dim
    """
    def __init__(self, type : int, pos : tuple, index : int, dim : int = 2) -> None:
        self._type = type
        self._NN = np.zeros(2*dim)
        self._pos = pos
        self._index = index

    
    #-----------Getters------------
    def get_type(self) -> int:
        return self._type
    
    def get_NN(self) -> np.ndarray:
        return self._NN
    
    def get_pos(self) -> tuple:
        return self._pos
    
    def get_index(self) -> int:
        return self._index
    
    #-----------Setters------------
    def set_pos(self, pos) -> None:
        self._pos = pos

    #-0-
    #3.1
    #-2-     . marks amino position, 3D gives 4 and 5 upwards then downwards
    def set_NN(self, amino, relative_pos : tuple):
        if len(relative_pos) == 2:
            match relative_pos:
                case (0,1):
                    self._NN[0] = amino
                case (1,0):
                    self._NN[1] = amino
                case (0,-1):
                    self._NN[2] = amino
                case (-1,0):
                    self._NN[3] = amino
        
        elif len(relative_pos) == 3:
             match relative_pos:
                case (0,1,0):
                    self._NN[0] = amino
                case (1,0,0):
                    self._NN[1] = amino
                case (0,-1,0):
                    self._NN[2] = amino
                case (-1,0,0):
                    self._NN[3] = amino
                case (0,0,1):
                    self._NN[4] = amino
                case (0,0,-1):
                    self._NN[5] = amino


    #-----------Help functions----------
    def is_nearest_neighbour(self, amino) -> bool:
        #Checks if they are covalently bonded
        if abs(self.get_index() - amino.get_index()) is not 1:
            #Substracting each x-y-z component to check distance vector
            difference = tuple(map(lambda i, j: i - j, self.get_pos, amino.get_pos()))

            #Sjekker om avstanden mellom de er 1, er den mer så er de ikke naboer, er den 0 så er de den samme monomeren
            if np.linalg.norm(difference) == 1:
                return True
        return False
    
    #def check_NN(self, amino) -> bool:

    