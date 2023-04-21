import init
import utility
import visualize
import numpy as np
import random

#parameters
from config import *

accepted_states = 0

def MC_step(positions, Npart):

    global accepted_states
    
    idx = random.randint(0,Npart-1)
    Uo = utility.energy_of_a_chain(idx, positions, box_length)
    Po = utility.virial_pressure_molecular(idx, positions)
    prev_position = np.copy(positions[idx])
    if molecule_type == 'methane':
        positions[idx] = init.generate_random_atom()
    elif molecule_type == 'ethane':
        positions[idx] = init.generate_random_molecule(bond_length)
    elif molecule_type == 'propane':
        positions[idx] = init.generate_propane_molecule(bond_length)
    elif molecule_type == 'butane':
        positions[idx] = init.generate_butane_molecule(bond_length)

    Un = utility.energy_of_a_chain(idx,positions,box_length)
    
    if random.random() < np.exp(-beta *(Un-Uo)):
        print(f"New state accepted")
        accepted_states += 1
        Pn = utility.virial_pressure_molecular(idx, positions)

        return positions, Un-Uo, Pn - Po
    else :
        positions[idx] = prev_position
        print(f"Old state accepted")
        return positions, 0, 0

positions = init.init_system(box_length, Npart)
energy = utility.total_energy(positions,box_length)
pressure = utility.calculate_pressure(positions) *(100/6.022) #converting into bar
energy_list = [energy]
accepted_states_list = []
pressure_list = [pressure]
for i in range(nsteps):
    positions, energy_change, pressure_change = MC_step(positions,Npart)
    # print(energy_change)
    energy+= energy_change
    pressure += pressure_change * (100/6.022) #converting into bar
    energy_list.append(energy)
    accepted_states_list.append(accepted_states*100/(i+1))
    pressure_list.append(pressure)
    print(f"Energy and Pressure of the system at {i} step is {energy} and {pressure} respectively")

np.save(f'{project_name}_{molecule_type}_energy.npy', np.array(energy_list, dtype=object), allow_pickle=True)
np.save(f'{project_name}_{molecule_type}_ar_ratio.npy', np.array(accepted_states_list, dtype=object), allow_pickle=True)
np.save(f'{project_name}_{molecule_type}_pressure.npy', np.array(pressure_list, dtype=object), allow_pickle=True)



