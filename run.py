from ethane_cbmc import *
from utility import *
file="output.out"
file_xyz="output_xyz.gro"
out_gro=open(file_xyz,'w')

out=open(file,'w')
out.write("# step  PE")
out.write("\n")
# positions = init.init_system(box_length, Npart)
positions = read_positions_from_gromacs_file('nvt.gro')
# visualize.visualize(positions)
print(len(positions))
energy = utility.total_energy(positions,box_length)
print(energy)
for i in range(nsteps):
    # visualize.visualize(positions)
    energy = utility.total_energy(positions,box_length)
    print(energy)
    positions = CBMC_step(positions,Npart)
    if (i%100==0):
        output_energy(i, energy,file)
        output_xyz(Npart,box_length,positions,i,file_xyz)
out.close()

