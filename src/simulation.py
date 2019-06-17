# In order to run as expected, this script requires that the system
# contain an installed copy of GROMACS, available for reference
# with the standard "gmx ..." syntax
import shutil
from shutil import copytree, copyfile
from zipfile import ZipFile
import subprocess, datetime
import os, statistics
import pymbar
import matplotlib.pyplot as pyplot
import numpy as np
import mdtraj as md

def build_directories(polymer_name,polymer_length,run_directory,fresh_run=False):
        """

        Parameters
        ----------

        """
        if not os.path.exists(run_directory):
          os.mkdir(run_directory)

        current_file = os.path.abspath(os.path.dirname(__file__))
        terphenyl_top = os.path.join(current_file, '../')
        template_input_files = str(str(polymer_name)+'/'+str(polymer_length)+'/input_files')
        input_files = str(str(run_directory)+'/input_files')
        gaff_directory = str(str(terphenyl_top)+'/gaff')
        gaff_run_directory = str(str(run_directory)+'/input_files/gaff')
        if os.path.exists(input_files):
          if fresh_run:
            shutil.rmtree(input_files)
            copytree(template_input_files,input_files)
          else:
            subprocess.run(["cp","-r",template_input_files,input_files])
        if not os.path.exists(str(str(run_directory)+'/input_files/gaff')):
          copytree(gaff_directory,gaff_run_directory)
        pdb_file = str(str(polymer_length)+".pdb")
        solvent_pdb_file_path = str(str(terphenyl_top)+"/"+str(polymer_name)+"/"+str(polymer_length)+"/input_files/solvent.pdb")
        pdb_path = str(str(terphenyl_top)+"/"+str(polymer_name)+"/"+str(polymer_length)+"/input_files/"+pdb_file)
        run_pdb_path = pdb_file
        copyfile(pdb_path,run_pdb_path)
        copyfile(solvent_pdb_file_path,str(str(input_files)+"/solvent.pdb"))
        return

def parameterize(polymer_length,polymer_code,polymer_abbreviation,pdb_file):
        """

        Parameters
        ----------

        """
        # Parameterize our polymer using 'antechamber', from AmberTools.
#
        # We parameterize the PDB structure using a BASH script written by Ben Coscia: "https://github.com/shirtsgroup/useful-scripts/blob/master/Paramaterization/GAFF/param.sh"
        #  print(os.getcwd())
        current_file = os.path.abspath(os.path.dirname(__file__))
        terphenyl_top = os.path.join(current_file, '../')
        gaff_directory = str(str(terphenyl_top)+"/gaff")
        copyfile(str(str(gaff_directory)+"/param.sh"),"param.sh")
        copyfile(str(str(gaff_directory)+"/acpype.py"),"acpype.py")
        copyfile(str(str(gaff_directory)+"/insertmol2charges.py"),"insertmol2charges.py")
        copyfile(str(str(gaff_directory)+"/em.mdp"),"em.mdp")
        copyfile(str(str(gaff_directory)+"/anneal.mdp"),"anneal.mdp")
        # Replace the variable keyword '$NAME' in param.sh with the name of the current polymer length
        with open("param.sh", "rt") as fin:
          with open("new_param.sh", "wt") as fout:
            for line in fin:
               fout.write(line.replace('$NAME', polymer_length).replace('$RES', polymer_code))
        os.rename("new_param.sh","param.sh")

        # Place the residue name in the input PDB file residue name columns
        with open(pdb_file, "rt") as fin:
          with open(str("new_"+pdb_file), "wt") as fout:
              for line in fin:
                  line_list = [char for char in line]
                  line_start = ''.join(line_list[1:6])
                  if line_start == 'HETATM' or line_start == 'ATOM  ':
                    line_list[18:20] = polymer_abbreviation
                    line = ''.join(line_list)
                  fout.write(line)
        os.rename(str("new_"+pdb_file),pdb_file)
        subprocess.run(["chmod","+x","param.sh"])
        subprocess.run(["./param.sh",pdb_file])

        return

def minimize():
        """

        Parameters
        ----------

        """
        # Setup an energy minimization of the initial structure guess
        subprocess.run(["gmx","grompp","-f","em.mdp","-p","topol.top","-c","solvated.gro","-o","em"])
        # Run the energy minimization
        subprocess.run(["gmx","mdrun","-v","-deffnm","em"])
        subprocess.run(["gmx","trjconv","-f","em.trr","-o","em.pdb"])
        return

def solvate(input_pdb,solvent_density=0.5):
        """

        Parameters
        ----------

        """
        # Build a simulation box and fill it with solvent
        subprocess.run(["gmx","solvate","-cp",input_pdb,"-cs","solvent.pdb","-p","topol.top","-o","solvated.gro"])

        return

def equilibrate():
        """

        Parameters
        ----------

        """
        # Setup an equilibration run with Berendsen barostat
        subprocess.run(["gmx","grompp","-f","berendsen.mdp","-p","topol.top","-c","em.gro","-o","berendsen"])
        # Run the equilibration
        subprocess.run(["gmx","mdrun","-v","-deffnm","berendsen"])
        subprocess.run(["gmx","trjconv","-f","berendsen.trr","-o","berendsen.pdb"])

        return

def simulate():
        """

        Parameters
        ----------

        """
        subprocess.run(["gmx","grompp","-f","npt.mdp","-p","topol.top","-c","berendsen.gro","-o","npt"])
        subprocess.run(["gmx","mdrun","-v","-deffnm","npt"])
        subprocess.run(["gmx","trjconv","-f","npt.trr","-o","npt.pdb"])
        return

def compress_large_files(directory,size_threshold=1.0e8):
        """

        Parameters
        ----------

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
