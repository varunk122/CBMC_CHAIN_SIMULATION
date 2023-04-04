#parameters
rho = 0.01      #density of the system
Npart = 108     #number of molecules
box_length =  2.2#(Npart / rho) ** (1/3)  #box length of the system
bond_length = 0.15  #angstrom       #bond length of the system
write_interval = 100            #write interval

k = 30         #number of trial orientations
temp = 300 #kelvin     #temperature of the system
beta = 1               #value of beta / KbT

sigma = 0.375   #for LJ potential
eps = 0.46   # e/KbT #for LJ Potential

nsteps = 10000  #total number of simulation steps
