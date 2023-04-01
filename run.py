from ethane_cbmc import *

positions = init.init_system(box_length, Npart)

for i in range(nsteps):
    # visualize.visualize(positions)
    energy = utility.total_energy(positions,box_length)
    print(f"Energy of the system at step {i} is: {energy}")
    
    positions = CBMC_step(positions,Npart)
