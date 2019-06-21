# In order to run as expected, this script requires that the system
# contain an installed copy of GROMACS, available for reference
# with the standard "gmx ..." syntax
from terphenyl_folding.src.simulation import *
from terphenyl_folding.src.analysis import *
import datetime
import os, statistics
import pymbar
import matplotlib.pyplot as pyplot
import numpy as np
import mdtraj as md

# Begin user input
polymer_name = "o-terphenyl"
polymer_length = "octamer"
polymer_abbreviation = ['O','C','T']
polymer_code = ''.join(polymer_abbreviation)
make_parameter_files = True
add_solvent = True
run_minimization = True
run_equilibration = True
run_simulation = True
analyze_simulation_data = True
# End user input

date = str(datetime.datetime.now()).split()[0]
input_directory = str(str(os.path.abspath(os.path.dirname(__file__)))+"/"+str(polymer_name)+'/'+str(polymer_length)+'/input_files')
run_directory = str(str(str(os.path.abspath(os.path.dirname(__file__)))+"/"+str(polymer_name)+'/'+str(polymer_length)+'/run_'+str(date)))
pdb_file = str(str(input_directory)+"/"+str(polymer_length)+".pdb")

build_directories(polymer_name,polymer_length,run_directory,fresh_run=False)
os.chdir(str(str(run_directory)+'/input_files'))

if make_parameter_files:
# Parameterize our polymer using 'antechamber', from AmberTools.
  parameterize(polymer_length,polymer_code,pdb_file,run_directory)

if add_solvent:
# Add solvent to a simulation box containing the system
  solvate(input_pdb="em2.gro",solvent_density=0.6)
exit()

if run_minimization:
# Minimize our initial structure
  minimize()

if run_equilibration:
# Run NPT equilibration of our minimized structure
  equilibrate()

if run_simulation:
# Run a Parrinello-Rahman pressure control simulation run
  simulate() 
  compress_large_files(run_directory) 

if analyze_simulation_data:
# Define the paths for our simulation output files
  gmx_simulation_energies = str(str(run_directory)+"/"+str(simulation)+".edr")
  gmx_trajectory = str(str(run_directory)+"/"+str(simulation)+".xtc")
# Read in simulation data
  traj = read_trajectory(pdb_file,gmx_trajectory)
  energies = read_energies(gmx_simulation_energies)
# Get equilibrium frames, and data
  equilibrium_frames = get_equilibrium_frames(energies)
  energies = energies[equilibrium_frames]
  traj = traj[equilibrium_frames]
# Analyze trajectory to get internal coordinates of interest
  torsion_definitions = get_torsion_definitions(pdb_file)

if archive:
  compress_large_files(run_directory)

exit()

# Read the energies
trajectory_file = ".xvg"
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

# Compress large files
filesList = []
for path, subdirs, files in os.walk(str(str(run_directory)+'/input_files')):
    for name in files:
        filesList.append(os.path.join(path, name))

for file in filesList:
    # Getting the size in a variable
    fileSize = os.path.getsize(str(file))
    print(fileSize)

    # Print the files that meet the condition
    if int(fileSize) >= int(1.0e8):
      print(file)
      shutil.make_archive(str(file),"zip",file)

exit()
