from shutil import copytree, copyfile
import subprocess, datetime
import os, statistics
import pymbar
import matplotlib.pyplot as pyplot
import numpy as np
import mdtraj as md

polymer_name = "tetramer"
date = str(datetime.datetime.now()).split()[0]
build_input_files = True
run_minimization = False
equilibrate = False
run_simulation = True
run_directory = str(str(polymer_name)+'/fresh_run_'+str(date))

if not os.path.exists(run_directory): 
  os.mkdir(run_directory)

if build_input_files:
  if not os.path.exists(str(str(run_directory)+'/input_files')):
    copytree('tetramer/input_files',str(str(run_directory)+'/input_files'))
  if not os.path.exists(str(str(run_directory)+'/input_files/gaff')):
    copytree('gaff',str(str(run_directory)+'/input_files/gaff'))
  os.chdir(str(str(run_directory)+'/input_files'))
# Make the topology file
  subprocess.run(["gmx","pdb2gmx","-f","tetramer.pdb","-o","em"])
  exit()

if run_minimization:
# Setup an energy minimization of the initial structure guess
  subprocess.run(["gmx","grompp","-f","em.mdp","-p","topol.top","-c","solvated.gro","-o","em"]) 
# Run the energy minimization
  subprocess.run(["gmx","mdrun","-v","-deffnm","em"])
  subprocess.run(["gmx","trjconv","-f","em.trr","-o","em.pdb"])
if equilibrate:
# Setup an equilibration run with Berendsen barostat
  subprocess.run(["gmx","grompp","-f","berendsen.mdp","-p","topol.top","-c","em.gro","-o","berendsen"])
# Run the equilibration
  subprocess.run(["gmx","mdrun","-v","-deffnm","berendsen"])
  subprocess.run(["gmx","trjconv","-f","berendsen.trr","-o","berendsen.pdb"])
if run_simulation:
# Setup a Parrinello-Rahman pressure control simulation run
  subprocess.run(["gmx","grompp","-f","npt.mdp","-p","topol.top","-c","berendsen.gro","-o","npt"])
  


exit()

# Read the energies
trajectory_file = "test.xvg"
file_obj = open(trajectory_file,"r")
lines = file_obj.readlines()
file_obj.close()
index = 0
# Determine energies array size
for line in lines:
  if not any([str(line[0]) == i for i in ["#","@"]]):
    index = index + 1
time = np.zeros(index)
potential_energy = np.zeros(index)

index = 0
for line in lines:
  if not any([str(line[0]) == i for i in ["#","@"]]):
    time[index] = line.split()[0]
    potential_energy[index] = line.split()[7]
    index = index + 1

figure_index = 1
figure = pyplot.figure(figure_index)
pyplot.plot(time,potential_energy,figure=figure)
pyplot.xlabel("Simulation Time (ps)")
pyplot.ylabel("Potential Energy (kJ/mol)")
pyplot.savefig(str("potential_energy.png"))
#pyplot.show()
pyplot.close()

# Get the minimum energy simulation frame
minimum_potential_index = np.argmin(potential_energy)
# Extract the coordinates for the minimum energy frame from the PDB file
#pdb_file = "test.pdb"
file = "minimum_energy_pose.pdb"
print("Reading the PDB file.")
lowest_energy_trajectory = md.load(file)
#lowest_energy_trajectory = md.load_frame(pdb_file,minimum_potential_index)
#file = "minimum_energy_pose.pdb"
#lowest_energy_trajectory.save_pdb(file)

# Define data structures needed to calculate the pi-stacking distance

# Define a set of atom labels/indices for each phenyl group in our model
particle_name_list = [str(atom).split('TMR1-')[1] for atom in lowest_energy_trajectory.topology.atoms]
polymer_length = 4
number_phenyl_groups_per_monomer = 3
index = 0
phenyl_list = []
# The center-of-mass coordinates for each phenyl group are categorized according
# to their index in the monomer (1, 2, or 3), as stacking would suggest there could be a relationship between the COM distances in specific phenyl group types.
com_coordinates_list = {'1': [], '2': [], '3': []}
for monomer_index in range(polymer_length):
  for phenyl_index in range(number_phenyl_groups_per_monomer):
    phenyl = []
    for carbon_index in range(6):
      if index == 0:
        particle_name = "C"
      else:
        particle_name = str("C"+str(index))
      for particle_index in range(len(particle_name_list)):
        if particle_name == str(particle_name_list[particle_index]):
          phenyl.append(particle_index)
      index = index + 1
    phenyl_list.append(phenyl)
    phenyl_coordinates_list = [lowest_energy_trajectory.xyz[0][i] for i in phenyl]
    com_coordinates = [float(sum([float(phenyl_coordinates_list[i][j]) for i in range(6)])/6.0) for j in range(3)]
    com_coordinates_list[str(phenyl_index+1)].append(com_coordinates)

# Calculate the hydrogen bond distance for 

exit()
