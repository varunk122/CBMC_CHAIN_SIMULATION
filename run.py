from ethane_cbmc import *

positions = init.init_system(box_length, Npart)
print(f"Number of molecules inserted {len(positions)}")
# positions = utility.read_positions_from_file("positions.npy")

energy_file = open("energy.txt",'a' )

for i in range(nsteps):
    # visualize.visualize(positions)
    energy = utility.total_energy(positions,box_length)
    print(f"Energy of the system at step {i} is: {energy}")
    if i%write_interval == 0:
        energy_file.write(str(energy)+"\n")
    positions = CBMC_step(positions,Npart)

np.save('positions.npy', np.array(positions, dtype=object), allow_pickle=True)
