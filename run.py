from ethane_cbmc import *
from utility import *
file="output4.out"
file_xyz="output_xyz4.gro"

out_gro=open(file_xyz,'w')
out=open(file,'w')
out.write("# step  PE\n")
init_pos =np.array(read_positions_from_gromacs_file('init.gro',1))
#print(init_pos)
positions = init_2.init_system(init_pos[0],box_length,Npart)
output_xyz(Npart,box_length,positions,0,'initial_strc.gro')
print(len(positions))
#positions = read_positions_from_gromacs_file('step9900.gro')
# visualize.visualize(positions)
#print(len(positions))
energy = utility.total_energy(positions,box_length)
print(energy)
for i in range(nsteps):
    # visualize.visualize(positions)
    energy = utility.total_energy(positions,box_length)
    print(energy)
    positions = CBMC_step(positions,Npart)
    if (i%write_interval==0):
        output_energy(i, energy,file)
        output_xyz(Npart,box_length,positions,i,file_xyz)
out.close()

