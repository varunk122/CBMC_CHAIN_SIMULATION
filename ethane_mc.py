import init
import utility
import visualize
import numpy as np
import random

#parameters

from config import *

def MC_step(positions, Npart):
    
    idx = random.randint(0,Npart-1)
    Uo = utility.energy_of_a_chain(idx, positions, box_length)
    prev_position = np.copy(positions[idx])
    positions[idx] = init.generate_random_molecule(bond_length)
    Un = utility.energy_of_a_chain(idx,positions,box_length)
    print(Un, Uo)
    if random.random() < np.exp(-beta *(Un-Uo)):
        print(f"New state accepted")
        return positions, Un-Uo
    else :
        positions[idx] = prev_position
        print(f"Old state accepted")
        return positions, 0

    

positions = init.init_system(box_length, Npart)
energy = utility.total_energy(positions,box_length)

for i in range(nsteps):
    positions, energy_change = MC_step(positions,Npart)
    print(energy_change)
    energy+= energy_change
    print("energies" ,energy, utility.total_energy(positions, box_length))


