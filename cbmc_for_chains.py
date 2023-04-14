import ethane_cbmc
import propane_cbmc
import init
import utility
import numpy as np

#import parameters
from config import *

positions  = init.init_system(box_length, Npart)
energy = utility.total_energy(positions, box_length)

energy_list = [energy]
ar_list = []

for i in range(nsteps):
    if molecule_type == 'ethane':
        positions, accepted_states = ethane_cbmc.CBMC_step(positions, Npart)
    elif  molecule_type == 'propane':
        positions, accepted_states = propane_cbmc.CBMC_step(positions, Npart)

    energy = utility.total_energy(positions, box_length)
    energy_list.append(energy)
    ar_list.append(accepted_states*100/(i+1))
    print(f"Energy of the system at {i} step is {energy}")

np.save(f'{project_name}_{molecule_type}_energy.npy', np.array(energy_list, dtype=object), allow_pickle=True)
np.save(f'{project_name}_{molecule_type}_ar_ratio.npy', np.array(ar_list, dtype=object), allow_pickle=True)


