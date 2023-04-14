#parameters
rho = 0.001      #density of the system
Npart = 108     #number of molecules
box_length = (Npart / rho) ** (1/3)  #box length of the system
bond_length = 1.5  #angstrom       #bond length of the system
write_interval = 100            #write interval

k = 30         #number of trial orientations
temp = 300 #kelvin     #temperature of the system
beta = 1               #value of beta / KbT

sigma = 3.75   #for LJ potential
eps = 0.46   # e/KbT #for LJ Potential

theta_mean = 1.936 #radian
k_prop = 212 # KT/rad^2 

nsteps = 100000  #total number of simulation steps

project_name = "run11" #name of the project
