#parameters
rho =   0.602  #density of the system
Npart = 108     #number of molecules
box_length =  (Npart / rho) ** (1/3)  #box length of the system
bond_length = 0.15  #nanometer       #bond length of the system
write_interval = 100            #write interval
k = 10         #number of trial orientations
temp = 180     #kelvin     #temperature of the system
beta = 1/(8.314*0.001 * temp)              #value of beta / KbT
sigma = [[0.375,0.385,0.385,0.375],[0.385,0.395,0.395,0.385],[0.385,0.395,0.395,0.385],[0.375,0.385,0.385,0.375]]   #for LJ potential
eps =  [[0.814,0.558,0.558,0.814],[0.558,0.382,0.382,0.558],[0.558,0.382,0.382,0.558],[0.814,0.558,0.558,0.814]]  # e/KbT #for LJ Potential
rcut = 3*0.375 #only for calculating pressure
theta_mean = 1.91 #radian
k_prop = 212 # KT/rad^2 

nsteps = 100000  #total number of simulation steps

project_name = "simulation_for_butane" #name of the project
molecule_type = "butane" #name of molecule to simulate

