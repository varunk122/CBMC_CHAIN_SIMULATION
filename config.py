#parameters
rho = 0.01      #density of the system
Npart = 108     #number of molecules
box_length = (Npart / rho) ** (1/3)  #box length of the system
bond_length = 2.353  #angstrom       #bond length of the system

k = 30         #number of trial orientations
temp = 300 #kelvin     #temperature of the system
beta = 1               #value of beta

nsteps = 5000         #total number of simulation steps
