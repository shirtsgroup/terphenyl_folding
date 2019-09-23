import shutil
from shutil import copytree, copyfile
from zipfile import ZipFile
import subprocess, datetime
import os, statistics
import pymbar
import matplotlib.pyplot as pyplot
import numpy as np
import mdtraj as md
import terphenyl_folding

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
        date = str(datetime.datetime.now()).split()[0]
        run_directory = str(os.path.abspath('../../')+'/data/'+str(polymer_name)+'/'+str(polymer_length)+'/run_'+str(date))

        if not os.path.exists(run_directory):
          os.makedirs(run_directory)
        else:
          if fresh_run:
            shutil.rmtree(run_directory)
            os.makedirs(run_directory)

        terphenyl_top = os.path.abspath('../../')
        input_files = str(str(terphenyl_top)+'/input_files/'+str(polymer_name)+'/'+str(polymer_length))

        pdb_file = str(str(input_files)+'/'+str(polymer_length)+'.pdb')
        solvent_file = str(str(input_files)+"/solvent.pdb")
        topology_file = str(str(terphenyl_top)+"/setup_files/topol.top")
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
        if not os.path.exists(param_directory):
          os.mkdir(param_directory)
        param_pdb = str(str(param_directory)+"/"+str(polymer_length)+".pdb")
        copyfile(pdb_file,param_pdb)

        # Parameterize our polymer using 'antechamber', from AmberTools.
#
        # We parameterize the PDB structure using the param.sh BASH script written by Ben Coscia as a template: "https://github.com/shirtsgroup/useful-scripts/blob/master/Paramaterization/GAFF/param.sh"
        terphenyl_top = os.path.abspath('../../')
        gaff_directory = str(str(terphenyl_top)+"/setup_files/gaff")
        param_topology = str(str(param_directory)+"/topol.top")
        copyfile(topology_file,param_topology)
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
       
        return(solute_gro_file)

def minimize(run_directory,topology,input_structure):
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
        mdrun_file = open("em.mdp","w")
        mdrun_file.write("title = Energy Minization\n")
        mdrun_file.write("integrator = steep\n")
        mdrun_file.write("nsteps = -1\n")
        mdrun_file.write("cutoff-scheme = verlet\n")
        mdrun_file.write("nstlist = 40")
        mdrun_file.close()

        subprocess.run(["gmx","grompp","-f",mdrun_file,"-p",topology,"-c",input_structure,"-o","em3"])
        # Run the energy minimization
        subprocess.run(["gmx","mdrun","-v","-deffnm","em3"])
        minimized_pdb_file = "em3.pdb"
        subprocess.run(["gmx","trjconv","-f","em3.trr","-s",input_file,"-o",minimized_pdb_file])
        return(minimized_pdb_file)

def solvate(solvation_directory,solute_gro_file,solute_topology_file,solvent_file,solvent_density=0.5):
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
        os.chdir(solvation_directory)
        copyfile(solute_gro_file,str(str(solvation_directory)+"/solute.gro"))

        solvated_gro_file = str(str(solvation_directory)+"/solvated.gro")
        subprocess.run(["gmx","solvate","-cp",solute_gro_file,"-cs",solvent_file,"-p",solute_topology_file,"-o",solvated_gro_file])
        os.chdir(cwd)
        return(solvated_gro_file)

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
        terphenyl_top = os.path.abspath('../../')
        ber_file = str(str(terphenyl_top)+"/input_files/berendsen.mdp")
        run_file = str(str(run_directory)+"/berendsen.mdp")
        copyfile(ber_file,run_file)
        subprocess.run(["gmx","grompp","-f",run_file,"-p",topology,"-c",input_gro_file,"-o","berendsen"])
        # Run the equilibration
        subprocess.run(["gmx","mdrun","-v","-deffnm","berendsen"])
        output_trajectory = "berendsen.pdb"
        subprocess.run(["gmx","trjconv","-f","berendsen.trr","-s",input_gro_file,"-o",output_trajectory])

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

def get_simulation_energies():
        return
