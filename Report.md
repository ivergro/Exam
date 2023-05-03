# General
All energy is in kb scale
Varying energy matrix and chain?
Increase in T leads to a increase in the Boltzmann factor, which makes it easier for unfavourable flips to happen(T gives energy to overcome energy barrier of bondings)

# Tasks
## Task 2.1.1
    

## Task 2.1.2

## Task 2.1.3
Negative J implies that the monomers favors attraction to minimize the total free energy in the system

## Task 2.1.4
Assuming tail and head moves like a free joint, so it can have two possible spots to move to.
If it has two possitions, i randomly chooses one (can take the one with smallest energy too)'

## Task 2.1.5
Plot 3 configs of protein after x sweeps. Use T=10
Plot energy, end-to-end distance and RoG as func of MC
Comment om tertiary structure that appears, link the two plots

From the 10 sweeps figure, one can se that the ends start moving. This is because the monomers that lay in a straight line is rigid. After some time, more of the monomers will start moving, creating bendt structures and forming non covalent bondings between neighbouring amino acids. One can see from the plot that the energy fluctuates between -16k_B and 0k_B, where 0 k_B marks no neighbouring interactions, apart from the already existing covalent bondings. The average energy seem to decrease after 50-60 sweeps, when the end-to-end distance decreases. This could mean that more neighbouring bonds can take form when the length of the protein decreases, since the spatial distance between each monomer decreases, and therefore increases the chance of an energy efficient flip.


## Task 2.1.6
Do same plots as in 2.1.5, but now use T=1
Comment on temp effect. 
How many steps before steady state? - just above 60 steps

A decrease in temp leads to an decrease in MC sweeps before reaching a steady state. Looking at the plot, the protein seems to have reached a steady state just after 60 sweeps, with only minor fluctuations occuring afterwards, because of the Metropolis MC method flipping some of the monomers. These monomers quickly flips back though, or some others may take their place so that the energy is minimized.

## Task 2.1.7
Use mean energy for each temp to plot phase diagram, E(T)



# Kilder
https://en.wikipedia.org/wiki/Center_of_mass center of mass
https://physicscatalyst.com/mech/radius-of-gyration.php RoG