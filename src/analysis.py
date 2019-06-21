
import mdtraj as md
import numpy as np
import panedr
import os
from pymbar import timeseries

def get_phenyl_carbon_indices(polymer_name,polymer_length):
    """
    Parameters
    ----------

    polymer_name: string (ie: "o-terphenyl")

    polymer_length: string (ie: "monomer")

    Returns
    -------

    phenyl_indices_list: List( dict('1': [phenyl-1 carbon indices], '2': [phenyl-2 carbon indices], '3': [phenyl-3 carbon indices]))

    """
    if polymer_name == "o-terphenyl":
      if polymer_length == "monomer":
        phenyl_indices_list = [{'1': [14,15,16,17,18,19], '2': [8,9,10,11,12,13], '3': [2,3,4,5,6,7]}]
      if polymer_length == "dimer":
        phenyl_indices_list = [{'1': [3,4,5,6,7,8], '2': [9,10,11,12,13,14], '3': [15,16,17,18,19,20]},{'1': [2,22,23,24,25,26], '2': [27,28,29,30,31,32], '3': [33,34,35,36,37,38]}]
      if polymer_length == "trimer":
        phenyl_indices_list = [{'1': [16,17,18,19,20,21], '2': [10,11,12,13,14,15], '3': [4,5,6,7,8,9]},{'1': [23,24,25,26,27,28], '2': [29,30,31,32,33,34], '3': [35,36,37,38,39,40]},{'1': [22,43,44,45,46,47], '2': [48,49,50,51,52,53], '3': [54,55,56,57,58,59]}]
      if polymer_length == "tetramer":
        phenyl_indices_list = [{'1': [11,21,22,17,8,16], '2': [7,20,12,15,19,9], '3': [4,6,18,13,14,10]},{'1': [37,38,39,40,41,42], '2': [31,32,33,34,35,36], '3': [25,26,27,28,29,30]},{'1': [44,45,46,47,48,49], '2': [50,51,52,53,54,55], '3': [56,57,58,59,60,61]},{'1': [43,64,65,66,67,68], '2': [69,70,71,72,73,74], '3': [75,76,77,78,79,80]}]
      if polymer_length == "octamer":
        phenyl_indices_list = [{'1': [72,73,83,78,86,77], '2': [9,74,82,79,75,85], '3': [15,20,84,80,81,76]},{'1': [52,53,54,55,56,57], '2': [58,59,60,61,62,63], '3': [64,65,66,67,68,69]},{'1': [45,46,47,48,49,50], '2': [39,40,41,42,43,44], '3': [33,34,35,36,37,38]},{'1': [13,14,31,28,184,181], '2': [19,22,29,26,27,30], '3': [10,11,23,24,25,18]},{'1': [87,95,90,101,102,96], '2': [91,94,98,88,165,99], '3': [89,93,92,97,159,154]},{'1': [117,118,119,120,121,122], '2': [111,112,113,114,115,116], '3': [105,106,107,108,109,110]},{'1': [124,125,126,127,128,129], '2': [130,131,132,133,134,135], '3': [136,137,138,139,140,141]},{'1': [123,146,143,160,161,157], '2': [144,147,148,145,152,155], '3': [149,156,164,163,151,150]}]
    return(phenyl_indices_list)

def get_phenyl_centers_of_mass(structure,polymer_name,polymer_length):
    """
    Parameters
    ----------

    structure: pdb file path

    polymer_name: string (ie: "o-terphenyl")

    polymer_length: string (ie: "monomer")

    """
    phenyl_indices_list = get_phenyl_carbon_indices(polymer_name,polymer_length)
    phenyl_com_coordinates = [{'1': [],'2': [],'3': []} for monomer in range(len(phenyl_indices_list))]
    # Iterate over monomers
    for monomer_index in range(len(phenyl_com_coordinates)):
      # Iterate over phenyl groups (3)
      for carbon_indices in phenyl_indices_list[monomer_index]:
        print(carbon_indices)
        exit()
        phenyl_coordinates = []
        # Iterate over carbon atoms
        for carbon in carbon_indices:
          phenyl_coordinates.append(structure[carbon])
        

    phenyl_list.append(phenyl)
    phenyl_coordinates_list = [lowest_energy_trajectory.xyz[0][i] for i in phenyl]
    com_coordinates = [float(sum([float(phenyl_coordinates_list[i][j]) for i in range(6)])/6.0) for j in range(3)]
    com_coordinates_list[str(phenyl_index+1)].append(com_coordinates)


    return(center_of_mass_coordinates)

def get_terphenyl_centers_of_mass(center_of_mass_coordinates):
    return()


def get_internal_coordinate_definitions(structure,polymer_name,polymer_length):
    """
    Parameters
    ----------

    structure: pdb file path

    polymer_name: string (ie: "o-terphenyl")

    polymer_length: string (ie: "monomer")

    """
    if polymer_name == "o-terphenyl":
      internal_coordinate_definitions = {'1': [], '2': [], '3': [], 'phi': [], 'psi': [], '4': [], '5': [], 'alpha_1': [], 'theta_1': [], 'alpha_2': [], 'theta_2': []}
      if polymer_length == "monomer":
        internal_coordinate_definitions['1'].append([15,14,13,8])
        internal_coordinate_definitions['2'].append([14,13,8,7])
        internal_coordinate_definitions['3'].append([1,7,8,13])
        
      if polymer_length == "dimer":
        internal_coordinate_definitions['1'].append([6,7,9,10],[24,25,27,28])
        internal_coordinate_definitions['2'].append([7,9,10,15],[25,27,28,33])
        internal_coordinate_definitions['3'].append([16,15,10,9],[34,33,28,27])
        internal_coordinate_definitions['phi'].append([18,21,0,39])
        internal_coordinate_definitions['psi'].append([23,22,39,0])
        
      if polymer_length == "trimer":
        internal_coordinate_definitions['1'].append([18,17,15,14],[26,27,29,30],[45,46,48,49])
        internal_coordinate_definitions['2'].append([17,15,14,9],[27,29,30,35],[46,48,49,54])
        internal_coordinate_definitions['3'].append([8,9,14,15],[27,29,30,25],[55,54,49,48])
        internal_coordinate_definitions['phi'].append([6,3,0,41],[38,42,1,60])
        internal_coordinate_definitions['psi'].append([0,41,24,25],[1,60,43,44])

      if polymer_length == "tetramer":
        internal_coordinate_definitions['1'].append([17,22,7,9],[39,38,36,35],[47,48,50,51],[66,67,69,70])
        internal_coordinate_definitions['2'].append([22,7,9,14],[38,36,35,30],[48,50,51,56],[67,69,70,75])
        internal_coordinate_definitions['3'].append([10,14,9,7],[29,30,35,36],[57,56,51,50],[76,75,70,69])
        internal_coordinate_definitions['phi'].append([6,5,1,24],[27,23,0,62],[59,63,2,81])
        internal_coordinate_definitions['psi'].append([40,41,24,1],[46,45,62,0],[65,64,81,2])
        internal_coordinate_definitions['alpha_2'].append([82,23,5],[23,5,63])
        internal_coordinate_definitions['theta_2'].append([82,23,5,63])

      if polymer_length == "octamer":
        internal_coordinate_definitions['1'].append([77,72,9,85],[55,56,58,59],[47,46,44,43],[14,31,30,27],[96,102,165,88],[119,118,116,115],[127,128,130,131],[160,143,144,147])
        internal_coordinate_definitions['2'].append([72,9,85,80],[56,58,59,64],[46,44,43,38],[31,30,27,25],[102,165,88,93],[118,116,115,110],[128,130,131,136],[143,144,147,149])
        internal_coordinate_definitions['3'].append([84,80,85,9],[65,64,59,58],[37,38,43,44],[24,25,27,30],[89,93,88,165],[109,110,115,116],[137,136,131,130],[150,149,147,144])
        internal_coordinate_definitions['phi'].append([15,12,2,70],[67,71,3,32],[35,16,1,8],[11,21,0,100],[159,162,5,104],[107,103,4,142],[139,158,6,166])
        internal_coordinate_definitions['psi'].append([2,70,53,54],[3,32,49,48],[1,8,17,13],[0,100,95,87],[5,104,121,120],[4,142,125,126],[6,166,157,161])
        internal_coordinate_definitions['alpha_2'].append([12,71,16],[71,16,21],[16,21,162],[21,162,103],[162,103,158],[103,158,153])
        internal_coordinate_definitions['theta_2'].append([12,71,16,21],[71,16,21,162],[16,21,162,103],[21,162,103,158],[162,103,158,153])


    return(internal_coordinate_definitions)

def read_trajectory(structure, trajectory):
    """
    Parameters
    ----------

    structure: pdb file path

    trajectory: xtc file path

    """
    traj = md.load(trajectory, top=structure)
    return(traj)

def read_energies(edr_file):
    energies = panedr.edr_to_df(edr_file)
    return(energies)

def get_equilibrium_frames(energy):
    t0, g, Neff_max = timeseries.detectEquilibration(energy['Potential'].values)
    return(t0)

def construct_selector(atom_list, res_dict, n_residues):
    """
    Since our molecule is all in one residue, we need some interesting ways to extract residues, instead of a simple
    select string

    atom_list : list
        String name of base atoms in the first residue
    res_dict : dict
        Dictionary telling how many atoms of a given element are in each residue
    """
    atom_info = [[atom[0], int(atom[1:])] for atom in atom_list]

    selector_list = []
    for i in range(n_residues):
        atomnos = [atom[0]+str(atom[1] + i*res_dict[atom[0]]) for atom in atom_info]
        # print(atomnos)
        for j in range(len(atomnos)):
            selector_list.append(atomnos[j])
    
    selector = ' '.join(selector_list)
    return(selector)

def get_dihedrals(selector, traj):
    top = traj.topology
    i_dihe = top.select(selector)
    i_dihe = i_dihe.reshape(int(len(i_dihe)/4), 4)
    dihes = md.compute_dihedrals(traj, i_dihe, periodic=True).flatten()
    return(dihes)


def main():
    import matplotlib.pyplot as plt

    traj_files = ['../octamer/initial_configs/1/PR.trr', 
                  '../octamer/initial_configs/3/PR.trr',
                  '../octamer/initial_configs/4/PR.trr',
                  '../tetramer/PR.trr']
    struct_files = ['../heteropolymer_inputs/octamer/initial_configs/frame_433.gro',
                    '../heteropolymer_inputs/octamer/initial_configs/frame_618.gro',
                    '../heteropolymer_inputs/octamer/initial_configs/solvated.gro',
                    '../heteropolymer_inputs/tetramer/solvated.gro']
    energy_files = ['../octamer/initial_configs/1/PR.edr',
                    '../octamer/initial_configs/3/PR.edr',
                    '../octamer/initial_configs/4/PR.edr',
                    '../tetramer/PR.edr']
    fig_names = ['outputs/octamer_config_1_dihe_',
                'outputs/octamer_config_3_dihe_',
                'outputs/octamer_config_4_dihe_',
                'outputs/tetramer_dihe_']
    selectors = [['C3', 'C4','C6','C7'], ['C4', 'C6', 'C7', 'C12'], ['C6', 'C7', 'C12', 'C13']]

    chain_sizes = [8, 8, 8, 4]

    dihe_nums = ['1','2','3']

    for traj_file, struct_file, energy_file, fig_name, chain_size in zip(traj_files, struct_files, energy_files, fig_names, chain_sizes):
        
        if not np.all([os.path.exists(f) for f in [fig_name+dihe_num+'.npy' for dihe_num in dihe_nums]]):
            traj = read_in_trajectory(struct_file, traj_file)  # expensive
            energy = read_in_energy(energy_file)
            t0 = get_equilibrium_frames(energy)

        # Dihe of interest : C3, C4, C6, C7
        # Dihe of interest : C4, C6, C7, C12
        # Dihe of interest : C6, C7, C12, C13
        for sel, dihe_num in zip(selectors, dihe_nums):
            selector = construct_selector(sel, {'C':20, 'N':1, 'O':1}, chain_size)
            if os.path.exists(fig_name+dihe_num+'.npy'):
                dihes = np.load(fig_name+dihe_num+'.npy')
            else:
                dihes = get_dihedrals('name '+selector, traj[t0:])*180/np.pi # expensive
                np.save(fig_name+dihe_num+'.npy', dihes)
            plt.figure()
            plt.hist(dihes, bins = 30, density=True)
            plt.ylabel('Probability Density')
            plt.xlabel('Dihedral Angle')
            plt.savefig(fig_name+dihe_num+'.pdf')
            plt.savefig(fig_name+dihe_num+'.png')

            # plt.show()

if __name__ == '__main__':
    main()
