import numpy as np

#Kan evt sjekke om avstanden til punktene er lik en lattice konstanten, heller enn om de er 1 unna hverandre
#Kan gjøre sånn at når en amino blir en annens NN, får den andre også den sin NN
#Må kanskje ikke ha NN avhening av pos, kan også bare ha en liste
#NN må skille på chain og bare nearest neigbhour
#Må kunne fjærne nærmeste aminoer
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
        assert (type >= 1 and type <= 20), f"Type got value: {type}, only valid between 1-20."
        self._type = type
        #self._NN = np.zeros(2*dim)
        self._NN = []
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
    
    def set_NN(self, amino):
        (self._NN).append(amino)

    def remove_NN(self, amino):
        (self._NN).remove(amino)
        return amino

    #-0-
    #3.1
    #-2-     . marks amino position, 3D gives 4 and 5 upwards then downwards
    #Blir nok tungvindt, må bytte pos om en flytter på seg
    # def set_NN(self, amino, relative_pos : tuple):
    #     if len(relative_pos) == 2:
    #         match relative_pos:
    #             case (0,1):
    #                 self._NN[0] = amino
    #             case (1,0):
    #                 self._NN[1] = amino
    #             case (0,-1):
    #                 self._NN[2] = amino
    #             case (-1,0):
    #                 self._NN[3] = amino
        
    #     elif len(relative_pos) == 3:
    #          match relative_pos:
    #             case (0,1,0):
    #                 self._NN[0] = amino
    #             case (1,0,0):
    #                 self._NN[1] = amino
    #             case (0,-1,0):
    #                 self._NN[2] = amino
    #             case (-1,0,0):
    #                 self._NN[3] = amino
    #             case (0,0,1):
    #                 self._NN[4] = amino
    #             case (0,0,-1):
    #                 self._NN[5] = amino


    #-----------Help functions----------
    def is_nearest_neighbour(self, amino) -> bool:
        #Checks if they are covalently bonded
        if (abs(self.get_index() - amino.get_index()) != 1):
            #Substracting each x-y-z component to check distance vector
            difference = tuple(np.subtract(self.get_pos(), amino.get_pos()))

            #Sjekker om avstanden mellom de er 1, er den mer så er de ikke naboer, er den 0 så er de den samme monomeren
            if np.linalg.norm(difference) == 1:
                return True                                 
        return False
    
    def find_nearest_neighbours(self, chain) -> None:
        for amino in chain:
            #Checks if theyre neigbhours, and if theyre already in the list
            if (self.is_nearest_neighbour(amino)) and (amino not in self.get_NN()):
                #Sets both as each others NN
                self.set_NN(amino)
                amino.set_NN(self)

    #Could make a tempchain, pop NN values, and loop over it the second time to minimize runtime
    def update_NN(self, chain) -> None:
        #Removing old ones first, looping backwards to avoid removing the first object, and never reach the second, and so fort
        for i in range(len(self.get_NN())-1,-1, -1):
            print(i)
            old_NN = self.get_NN()[i]
            if not self.is_nearest_neighbour(old_NN):
                self.remove_NN(old_NN)
                old_NN.remove_NN(self)
        
        #Updating
        #for amino in chain:
        #     if self.is_nearest_neighbour(amino):
        #         self.set_NN(amino)
        self.find_nearest_neighbours(chain)

    def does_collide(self, chain) -> bool:
        """
        Returns true if monomer collides with chain
        """
        for amino in chain:
            if self.get_pos() == amino.get_pos() and self.get_index() != amino.get_index():
                return True
        return False
    
    def can_move_to(self, chain):
        index = self.get_index()
        # if index != 0 and index != len(chain) -1:
        #     pass
        old_pos = self.get_pos()
        occupied_spots = [amino.get_pos() for amino in chain]
        available_spots = []

        #Tail or head
        if index == 0 or index == (len(chain) - 1):
            if index == 0:
                neighbour_index = index + 1
            else:
                neighbour_index = index - 1
            #Neighbour spots of previous monomer
            neighbour_spots = [(1,0),(0,1),(-1,0),(0,-1)]
            for ns in neighbour_spots:
                temp_p = tuple(np.add(ns, chain[neighbour_index].get_pos()))
                if occupied_spots.count(temp_p) == 0:
                    #print(np.add(ns, chain[neighbour_index].get_pos()))
                    available_spots.append(temp_p)

        #Cornered amino
        elif np.linalg.norm(np.subtract(chain[index + 1].get_pos(), chain[index - 1].get_pos())) < 2:
            #New pos takes both self pos and new pos just to make it general
            new_pos = [(chain[index - 1].get_pos()[0],chain[index + 1].get_pos()[1]), (chain[index + 1].get_pos()[0],chain[index - 1].get_pos()[1])]
            for new_p in new_pos:
                if occupied_spots.count(new_p) == 0:
                    available_spots.append(new_p)
  
        return available_spots

    