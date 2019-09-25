import shutil, os
from shutil import copytree, copyfile
from zipfile import ZipFile
import subprocess, datetime
import os, statistics, random
import pymbar
import matplotlib.pyplot as pyplot
import numpy as np
import mdtraj as md
import terphenyl_folding

def get_terphenyl_top_directory():
        """
        """
        terphenyl_top = os.environ['TERPHENYL_TOP']
        return(terphenyl_top)

def replace(file,original_text,replacement_text):
        """
        Given a file, a target search string, and a replacement string, this function replaces the text in 'file'.

        :param file: A file containing the text that will be replaced.
        :type file: file

        :param original_text: Text that will be replaced.
        :type original_text: str

        :param replacement_text: Text that will be used to replace the original text.
        :type replacement_text: str
        """
        with open(file, "rt") as fin:
          with open(str(file+"temp"), "wt") as fout:
            for line in fin:
               fout.write(line.replace(original_text,replacement_text))
        os.rename(str(file+"temp"),file)
        return

def make_topology(polymer_code,num_solvent_molecules=0):
        """
        """
        topology_file = "topol.top"
        file = open(topology_file,"w")
        file.write(";Forcefield\n")
        file.write('#include "gaff.itp"\n')
        file.write('\n')
        file.write(str(';'+str(polymer_code)+' Topology\n'))
        file.write('#include "'+str(polymer_code)+'.itp"\n')
        file.write('\n')
        if num_solvent_molecules != 0:
          file.write('TCM Topology\n')
          file.write('#include "TCM.itp"\n')
          file.write('\n')
        file.write('[ system ]\n')
        file.write('Simulation Box\n')
        file.write('\n')
        file.write('[ molecules ]\n')
        file.write(";Compounds     nmols\n")
        file.write(str(polymer_code)+"                1\n")
        if num_solvent_molecules != 0:
          file.write('TCM            '+str(str("{:>5}".format(num_solvent_molecules)))+'\n')
        file.close()
        return(topology_file)

def build_directories(polymer_name,polymer_length,fresh_run=False):
        """
        Given a set of input strings, this function builds the directories that are needed to perform GROMACS simulations with terphenyl oligomers.

        :param polymer_name: The name of the polymer
        :type polymer_name: str

        :param polymer_length: The length of the polymer we are modeling (in monomer units)
        :type polymer_length: str

        :param run_directory: The directory where simulations will be run
        :type run_directory: str

        :param fresh_run: A logical variable determining whether old run files should be removed.
        :type fresh_run: Logical

        :returns:
            - run_directory ( str ) - The path to a directory where simulations will be run.
            - pdb_file ( str ) - The path to pdb file that will be used for simulations.
            - solvent_file ( str ) - The path to a file containing a box of solvent molecules
            - topology_file ( str ) - The path to a file containing the topology for the pdb_file.
        """
        polymer_code=polymer_name[0:2].upper()
        date = str(datetime.datetime.now()).split()[0]
        terphenyl_top = get_terphenyl_top_directory()
        run_directory = str(str(terphenyl_top)+'/data/'+str(polymer_name)+'/'+str(polymer_length)+'/run_'+str(date))

        if not os.path.exists(run_directory):
          os.makedirs(run_directory)
        else:
          if fresh_run:
            shutil.rmtree(run_directory)
            os.makedirs(run_directory)

        input_files = str(str(terphenyl_top)+'/input_files/'+str(polymer_name)+'/'+str(polymer_length))

        pdb_file = str(str(input_files)+'/'+str(polymer_length)+'.pdb')
        solvent_file = str(str(input_files)+"/solvent.pdb")
        topology_file = make_topology(polymer_code)
        run_pdb_file = str(str(run_directory)+"/"+str(polymer_length)+".pdb")
        run_topology_file = str(str(run_directory)+"/topol.top")
        copyfile(pdb_file,run_pdb_file)
        copyfile(topology_file,run_topology_file)
        return(run_directory,pdb_file,solvent_file,topology_file)

def parameterize(param_directory,pdb_file,topology_file,polymer_code,polymer_length):
        """
        Given a directory path, PDB file, and a topology file, this function parameterizes the structure with GAFF.

        :param param_directory: The path to a directory where intermediate and output parameterization files will be written.
        :type param_directory: str

        :param pdb_file: The path to a PDB file containing data for the structure that will be parameterized.
        :type pdb_file: str

        :param topology_file: The path to a file containing the topology for the structure that will be parameterized.
        :type topology_file: str

        :param polymer_code: A three-letter code for the polymer (required for GAFF parameterization, in place of residue codes)
        :type polymer_code: str (length=3)

        :param polymer_length: The length of the input structure (number of monomers).
        :type polymer_length: str

        :returns:
          - solute_gro_file ( str ) - Path to the GAFF-parameterized .gro file for the input structure.
        """

        terphenyl_top = get_terphenyl_top_directory()
        if not os.path.exists(param_directory):
          os.mkdir(param_directory)
        param_topology = str(str(param_directory)+"/topol.top")
        copyfile(topology_file,param_topology)
        cwd = os.getcwd()
        if cwd != param_directory:
          os.chdir(param_directory)
        param_pdb = str(str(param_directory)+"/"+str(polymer_length)+".pdb")
        copyfile(pdb_file,param_pdb)

        # Parameterize our polymer using 'antechamber', from AmberTools.
#
        # We parameterize the PDB structure using the param.sh BASH script written by Ben Coscia as a template: "https://github.com/shirtsgroup/useful-scripts/blob/master/Paramaterization/GAFF/param.sh"
        gaff_directory = str(str(terphenyl_top)+"/setup_files/gaff")
        replace(param_topology,'$TERPHENYL_TOP',terphenyl_top)
        replace(param_topology,'$RUN_DIRECTORY',param_directory)
        replace(param_topology,'$POLYMER_CODE  ',str("{:<15}".format(polymer_code)))
        replace(param_topology,'$POLYMER_CODE ',str("{:<3}".format(polymer_code)))
        copyfile(str(str(gaff_directory)+"/acpype.py"),str(str(param_directory)+"/acpype.py"))
        copyfile(str(str(gaff_directory)+"/insertmol2charges.py"),str(str(param_directory)+"/insertmol2charges.py"))
#        copyfile(str(str(gaff_directory)+"/anneal.mdp"),str(run_directory+"/anneal.mdp"))
        # Replace the variable keyword '$NAME' in param.sh with the name of the current polymer length
        copyfile(str(str(gaff_directory)+"/param.sh"),str(str(param_directory)+"/param.sh"))
        replace(str(param_directory+"/param.sh"),'$NAME',polymer_length)
        replace(str(param_directory+"/param.sh"),'$RES',polymer_code)
       # Place the residue name in the input PDB file residue name columns
        with open(pdb_file, "rt") as fin:

          new_pdb_file = param_pdb
          with open(new_pdb_file, "wt") as fout:
              for line in fin:
                  line_list = [char for char in line]
                  line_start = ''.join(line_list[0:6])
                  residue_code = ''.join(line_list[17:20])
                  if line_start == 'HETATM' or line_start == 'ATOM  ':
                   if residue_code == '   ':
                    line_list[17:20] = str("{:<3}".format(polymer_code)).split()
                   #del line_list[29]
                  line = ''.join(line_list)
                  fout.write(line)
        subprocess.run(["chmod","+x",str(str(param_directory)+"/param.sh")])
        os.chdir(param_directory)
        subprocess.run([str(str(param_directory)+"/param.sh")])
        solute_gro_file = str(str(param_directory)+"/"+str(polymer_length)+".gro")
        solute_topology_file = str(str(param_directory)+"/"+str(polymer_code)+".top")
        if cwd != param_directory:
          os.chdir(cwd)
        return(solute_gro_file,solute_topology_file)

def minimize(run_directory,topology,input_structure,polymer_code):
        """
        Given a GROMACS .mdp MD commands file, a topology, and an input structure, this function performs a mimimization.

        :param run_directory: The path to a directory where the minimization will be run.
        :type run_directory: str

        :param topology: The path to a topology file
        :type topology: str

        :param input_structure: The path to an input structure to minimize.
        :type input_structure: type

        :returns:
           - minimized_pdb_file ( str ) - The path to a PDB file for the minimized structure. 
        """
        # Setup an energy minimization of the initial structure guess
        cwd = os.getcwd()
        if cwd != run_directory:
          os.chdir(run_directory)
        mdrun_file = open("em.mdp","w")
        mdrun_file.write("title = Energy Minization\n")
        mdrun_file.write("integrator = steep\n")
        mdrun_file.write("nsteps = -1\n")
        mdrun_file.write("cutoff-scheme = verlet\n")
        mdrun_file.write("nstlist = 40")
        mdrun_file.close()

        terphenyl_top = get_terphenyl_top_directory()
        gaff_directory = str(str(terphenyl_top)+"/setup_files/gaff")
        copyfile(str(str(gaff_directory)+"/gaff.itp"),str(str(run_directory)+"/gaff.itp"))
        copyfile(str(str(gaff_directory)+"/ffbonded.itp"),str(str(run_directory)+"/ffbonded.itp"))
        copyfile(str(str(gaff_directory)+"/ffnonbonded.itp"),str(str(run_directory)+"/ffnonbonded.itp"))
        copyfile(topology,str(str(run_directory)+"/"+str(polymer_code)+".itp"))
        copyfile(str(str(terphenyl_top)+"/setup_files/TCM.itp"),str(str(run_directory)+"/TCM.itp"))
        subprocess.run(["gmx","grompp","-f",mdrun_file,"-p",topology,"-c",input_structure,"-o","em3"])
        # Run the energy minimization
        subprocess.run(["gmx","mdrun","-v","-deffnm","em3"])
        minimized_pdb_file = "em3.pdb"
        subprocess.run(["gmx","trjconv","-f","em3.trr","-s",input_file,"-o",minimized_pdb_file])
        if cwd != run_directory:
          os.chdir(cwd)
        return(minimized_pdb_file)

def get_box_volume(solvent_file):
        """
        Given an input file containing data for a simulation box of solvent molecules, this function calculates the box volume.
        """
        box_volume = None
        file = solvent_file
        with open(file,"rt") as fin:
          for line in fin:
            if line[0:6] == "CRYST1":
              x_length = float(line[9:14])
              y_length = float(line[18:23])
              z_length = float(line[27:33])
              box_volume = x_length * y_length * z_length
              return(box_volume)
        return(box_volume)

def get_num_solvent_molecules(solvent_file):
        """
        
        """
        num_solvent_molecules = 0
        last_molecule_number = None
        file = solvent_file
        with open(file,"rt") as fin:
          for line in fin:
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
              molecule_number = int(line[22:26].strip())
              if last_molecule_number != None:
                if molecule_number != last_molecule_number:
                  last_molecule_number = molecule_number
                  num_solvent_molecules = num_solvent_molecules + 1
              else:
                last_molecule_number = molecule_number
        return(num_solvent_molecules)

def calculate_solvent_density(solvent_file):
        """
        Given an input file containing solvent molecules, this function computes the solvent density, in 
        """
        num_solvent_molecules = get_num_solvent_molecules(solvent_file)
        moles_solvent_molecules = num_solvent_molecules/6.02e23
        box_volume_cubic_angstroms = get_box_volume(solvent_file)
        if box_volume_cubic_angstroms == None:
          print("ERROR: The box volume wasn't calculated correctly.")
          exit()
        liter_conversion = 1.0e-27
        box_volume_liters = liter_conversion*box_volume_cubic_angstroms
        solvent_density = moles_solvent_molecules/box_volume_liters
        return(solvent_density)

def remove_random_molecules(solvent_file,num_molecules_to_remove):
        """
        """
        num_solvent_molecules = get_num_solvent_molecules(solvent_file)
        molecules_to_remove = []
        while len(molecules_to_remove) < num_molecules_to_remove:
          to_remove = random.randint(1,num_solvent_molecules)
          while to_remove in molecules_to_remove:
            to_remove = random.randint(1,num_solvent_molecules)
          molecules_to_remove.append(to_remove)

        num_solvent_molecules = 1
        last_molecule_number = 1
        file = solvent_file
        #print("Removing solvent molecules with the following IDs: "+str(molecules_to_remove))
        new_file = "new.pdb"
        with open(file,"rt") as fin:
         with open(new_file,"wt") as fout:
          for line in fin:
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
              molecule_number = int(line[22:26].strip())
              #print(molecule_number)
              if molecule_number not in molecules_to_remove:
               if molecule_number != last_molecule_number:
                last_molecule_number = molecule_number
                num_solvent_molecules = num_solvent_molecules + 1
               line_list = [char for char in line]
               line_list[22:26] = str("{:>4}".format(num_solvent_molecules)).split()
               line = ''.join(line_list)
               fout.write(line)
            else:
              fout.write(line)
        return(new_file)


def adjust_solvent_density(solvent_file,target_solvent_density):
        """
        """
        new_solvent_file = None
        if str(solvent_file.split('.')[1]) == "gro":
          subprocess.run(["gmx","editconf","-f",str(solvent_file),"-o",str(str(solvent_file.split('.')[0])+".pdb")])
          solvent_file = str(str(solvent_file.split('.')[0])+".pdb")
        print("The solvent file is : "+str(solvent_file))
        solvent_density = calculate_solvent_density(solvent_file)
        print("The current solvent density is: "+str(solvent_density))
        print("The target solvent density is: "+str(target_solvent_density))
        if solvent_density != target_solvent_density:
          print("Attempting to adjust the solvent density")
          num_solvent_molecules = get_num_solvent_molecules(solvent_file)
          print("The current number of solvent molecules is: "+str(num_solvent_molecules))
          box_volume_cubic_angstroms = get_box_volume(solvent_file)
          print("The solvent box has a volume of: "+str(box_volume_cubic_angstroms)+" cubic angstroms.")
          liter_conversion = 1.0e-27
          moles_solvent_molecules = num_solvent_molecules/6.02e23
          box_volume_liters = liter_conversion*box_volume_cubic_angstroms
          target_num_moles = target_solvent_density*box_volume_liters
          print("The target number of moles is: "+str(target_num_moles))
          target_num_molecules = round(target_num_moles*6.02e23)
          print("The target number of molecules is: "+str(target_num_molecules))
          if int(target_num_molecules) < int(num_solvent_molecules):
            num_molecules_to_remove = int(num_solvent_molecules)-int(target_num_molecules)
            print("Removing "+str(num_molecules_to_remove)+" molecules from the solvent box")
            new_solvent_file = remove_random_molecules(solvent_file,num_molecules_to_remove)
          if int(target_num_molecules) > int(num_solvent_molecules):
            num_molecules_to_add = int(target_num_molecules)-int(num_solvent_molecules)
            print("Adding "+str(num_molecules_to_add)+" molecules to the solvent box")
            new_solvent_file = add_random_molecules(solvent_file,num_molecules_to_add)
        else:
          print("The input solvent file had a density that is equivalent to 'target_solvent_density'")

        if new_solvent_file == None and solvent_density != target_solvent_density:
          print("ERROR: A new solvent file was not created correctly in 'adjust_solvent_density'")
          exit()
        num_solvent_molecules = get_num_solvent_molecules(new_solvent_file)
        if int(num_solvent_molecules) not in range(target_num_molecules-10,target_num_molecules+10):
            print("ERROR: Something went wrong when building a new solvent box.")
            exit()
        return(new_solvent_file)

def solvate(solvation_directory,solute_gro_file,solute_topology_file,solvent_file,polymer_code,solvent_density=None):
        """
        Given solute and solvent .gro files, this function produces a combined, solvated .gro file.

        :param solvation_directory: The path to a directory where solvation will be performed.
        :type solvation_directory: str

        :param solute_gro_file: The path to a .gro file for the solute
        :type solute_gro_file: str

        :param solute_topology_file: The path to a topology file for the solute.
        :type solute_topology_file: str

        :param solvent_file: The path to a .gro file containing solvent
        :type solvent_file: str

        :param solvent_density: The density for the solvent box.
        :type solvent_density: float

        :returns:
           - solvated_gro_file ( str ) - The path to a .gro file containing the combined, solvated system.
        """
        # Build a simulation box and fill it with solvent
        if not os.path.exists(solvation_directory):
          os.mkdir(solvation_directory)
        cwd = os.getcwd()
        solute_top = str(str(solvation_directory)+"/topol.top")
        copyfile(solute_topology_file,solute_top)
        if cwd != solvation_directory:
          os.chdir(solvation_directory)
        new_gro_file = str(str(solvation_directory)+"/solute.gro")
        copyfile(solute_gro_file,new_gro_file)

        solvated_gro_file = str(str(solvation_directory)+"/solvated.gro")
        subprocess.run(["gmx","solvate","-cp",new_gro_file,"-cs",solvent_file,"-p",solute_top,"-o",solvated_gro_file])
        if not os.path.exists(solvated_gro_file):
          print("ERROR: Something went wrong while solvating the molecule.\n")
          exit()
        if solvent_density != None:
          solvated_gro_file = adjust_solvent_density(solvated_gro_file,solvent_density)
        num_solvent_molecules = get_num_solvent_molecules(solvated_gro_file)
        solvated_topology_file = make_topology(polymer_code,num_solvent_molecules=num_solvent_molecules)
        if cwd != solvation_directory:
          os.chdir(cwd)
        return(solvated_gro_file,solvated_topology_file)

def equilibrate(run_directory,input_gro_file):
        """
        Given an input system, this function performs pressure (Berendsen) equilibration.

        :param run_directory: The path to a directory where the equilibration will be run.
        :type run_directory: str

        :param input_gro_file: The path to a .gro file that will be used for the equilibration run.
        :type input_gro_file: str

        :returns:
          - output_trajectory ( str ) - The path to a PDB file containing the output for the equilibration run.        
        """
        # Setup an equilibration run with Berendsen barostat
        cwd = os.getcwd()
        if cwd != run_directory:
          os.chdir(run_directory)
        terphenyl_top = os.path.abspath('../../')
        ber_file = str(str(terphenyl_top)+"/input_files/berendsen.mdp")
        run_file = str(str(run_directory)+"/berendsen.mdp")
        copyfile(ber_file,run_file)
        subprocess.run(["gmx","grompp","-f",run_file,"-p",topology,"-c",input_gro_file,"-o","berendsen"])
        # Run the equilibration
        subprocess.run(["gmx","mdrun","-v","-deffnm","berendsen"])
        output_trajectory = "berendsen.pdb"
        subprocess.run(["gmx","trjconv","-f","berendsen.trr","-s",input_gro_file,"-o",output_trajectory])
        if cwd != run_directory:
          os.chdir(cwd)
        return(output_trajectory)

def simulate():
        """
        """
        subprocess.run(["gmx","grompp","-f","npt.mdp","-p","topol.top","-c","berendsen.gro","-o","npt"])
        subprocess.run(["gmx","mdrun","-v","-deffnm","npt"])
        subprocess.run(["gmx","trjconv","-f","npt.trr","-o","npt.pdb"])
        return

def compress_large_files(directory,size_threshold=1.0e8):
        """
        """
        filesList = []
        for path, subdirs, files in os.walk(directory):
          for name in files:
            filesList.append(os.path.join(path, name))

        for file in filesList:
        # Getting the size in a variable
            fileSize = os.path.getsize(str(file))

            # Print the files that meet the condition
            if int(fileSize) >= int(size_threshold):
               shutil.make_archive(str(file),"zip",file)

        return
