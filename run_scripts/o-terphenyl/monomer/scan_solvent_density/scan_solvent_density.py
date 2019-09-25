from terphenyl_folding.simulation import *
from terphenyl_folding.analysis import *
import os, statistics
import matplotlib.pyplot as pyplot
import numpy as np

polymer_name = "o-terphenyl"
polymer_length = "monomer"
polymer_abbreviation = ['M','O','N']
polymer_code = ''.join(polymer_abbreviation)
solvent_density_list = [float(density*0.1) for density in range(1,100,10)]
fresh_run = True
make_parameter_files = True
add_solvent = True
run_minimization = True
run_equilibration = True
run_simulation = True
analyze_simulation_data = True

if fresh_run:
  run_directory,pdb_file,solvent_file,topology_file = build_directories(polymer_name,polymer_length,fresh_run=fresh_run)

if make_parameter_files:
# Parameterize our polymer using 'antechamber', from AmberTools.
  param_directory = str(str(run_directory)+"/parameterization")
  solute_gro_file = parameterize(param_directory,pdb_file,topology_file,polymer_code,polymer_length)

for solvent_density in solvent_density_list:
 if add_solvent:
# Add solvent to a simulation box containing the system
  solvation_directory = str(str(run_directory)+"/solvation")
  solvated_gro_file = solvate(solvation_directory,solute_gro_file,topology_file,solvent_file,solvent_density=solvent_density)

 if run_minimization:
# Minimize our initial structure
  minimized_pdb_file = minimize(solvation_directory,topology_file,solvated_gro_file)

 atom_pair = get_terminal_atoms(minimized_pdb_file)
 end_to_end_distance = get_end_to_end_distance(minimized_pdb_file,atom_pair)
 end_to_end_distance_list.append(end_to_end_distance)

solvent_density_list = np.array([float(solvent_density) for solvent_density in solvent_density_list])
end_to_end_distance_list = np.array([float(distance) for distance in end_to_end_distance_list])

figure_index = 1
figure = pyplot.figure(figure_index)
pyplot.plot(solvent_density_list,end_to_end_distance_list,figure=figure)
pyplot.xlabel("Solvent Density (mol/L)")
pyplot.ylabel("End-to-end distance (nm)")
pyplot.savefig(str("end-to-end-distance.png"))
pyplot.close()
