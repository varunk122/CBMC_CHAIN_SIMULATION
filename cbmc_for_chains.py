import ethane_cbmc
import propane_cbmc
import init
import utility
import numpy as np

#import parameters
from config import *

positions  = init.init_system(box_length, Npart)
energy = utility.total_energy(positions, box_length)
pressure = utility.calculate_pressure(positions)

energy_list = [energy]
ar_list = []
pressure_list = [pressure]

for i in range(nsteps):
    if molecule_type == 'ethane':
        positions, accepted_states = ethane_cbmc.CBMC_step(positions, Npart)
    elif  molecule_type == 'propane':
        positions, accepted_states = propane_cbmc.CBMC_step(positions, Npart)

    energy = utility.total_energy(positions, box_length)
    pressure = utility.calculate_pressure(positions)
    energy_list.append(energy)
    pressure_list.append(pressure)
    ar_list.append(accepted_states*100/(i+1))
    print(f"Energy and Pressure of the system at {i} step is {energy:.4f} and {pressure:.4f} respectively")

np.save(f'{project_name}_{molecule_type}_energy.npy', np.array(energy_list, dtype=object), allow_pickle=True)
np.save(f'{project_name}_{molecule_type}_ar_ratio.npy', np.array(ar_list, dtype=object), allow_pickle=True)
np.save(f'{project_name}_{molecule_type}_pressure.npy', np.array(pressure_list, dtype=object), allow_pickle=True)

