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

def replace(file,original_text,replacement_text):
        with open(file, "rt") as fin:
          with open(str(file+"temp"), "wt") as fout:
            for line in fin:
               fout.write(line.replace(original_text,replacement_text))
        os.rename(str(file+"temp"),file)
        return

def build_directories(polymer_name,polymer_length,run_directory,fresh_run=False):
        """

        Parameters
        ----------

        """
        if not os.path.exists(run_directory):
          os.mkdir(run_directory)
        else:
          if fresh_run:
            shutil.rmtree(run_directory)
            os.mkdir(run_directory)

        terphenyl_top = os.path.abspath(os.path.dirname(__file__)).split('/src')[0]
        input_files = str(str(terphenyl_top)+'/'+str(polymer_name)+'/'+str(polymer_length)+'/input_files')

        pdb_file = str(str(input_files)+'/'+str(polymer_length)+".pdb")
        solvent_pdb_file = str(str(input_files)+"/solvent.pdb")
        run_pdb_file = str(str(terphenyl_top)+"/"+str(polymer_name)+"/"+str(polymer_length)+"/run/"+str(polymer_length)+".pdb")
        run_solvent_file = str(str(terphenyl_top)+"/"+str(polymer_name)+"/"+str(polymer_length)+"/run/solvent.pdb")
        copyfile(pdb_file,run_pdb_file)
        copyfile(solvent_pdb_file,run_solvent_file)
        return

def parameterize(polymer_length,polymer_code,pdb_file,run_directory):
        """

        Parameters
        ----------

        """
        param_directory = str(run_directory+"/parameterization")
        if not os.path.exists(param_directory):
          os.mkdir(param_directory)
        os.chdir(param_directory)
        param_pdb = str(str(param_directory)+"/"+str(pdb_file.split('/run/input_files/')[1]))
        copyfile(pdb_file,param_pdb)

        # Parameterize our polymer using 'antechamber', from AmberTools.
#
        # We parameterize the PDB structure using the param.sh BASH script written by Ben Coscia as a template: "https://github.com/shirtsgroup/useful-scripts/blob/master/Paramaterization/GAFF/param.sh"
        terphenyl_top = str(os.path.abspath(os.path.dirname(__file__)).split('/src')[0])
        gaff_directory = str(terphenyl_top+"/setup_files/gaff")
        input_files = str(run_directory.split('/run')[0]+"/input_files")
        # Construct topology file
        topol = str(terphenyl_top+"/setup_files/topol.top")
        topol_new = str(param_directory+"/topol.top")
        copyfile(topol,topol_new)
        replace(topol_new,'$TERPHENYL_TOP',terphenyl_top)
        replace(topol_new,'$RUN_DIRECTORY',param_directory)
        replace(topol_new,'$POLYMER_CODE  ',str("{:<15}".format(polymer_code)))
        replace(topol_new,'$POLYMER_CODE ',str("{:<3}".format(polymer_code)))
        copyfile(str(gaff_directory+"/acpype.py"),str(param_directory+"/acpype.py"))
        copyfile(str(gaff_directory+"/insertmol2charges.py"),str(param_directory+"/insertmol2charges.py"))
#        copyfile(str(str(gaff_directory)+"/anneal.mdp"),str(run_directory+"/anneal.mdp"))
        # Replace the variable keyword '$NAME' in param.sh with the name of the current polymer length
        copyfile(str(gaff_directory+"/param.sh"),str(param_directory+"/param.sh"))
        replace(str(param_directory+"/param.sh"),'$NAME',polymer_length)
        replace(str(param_directory+"/param.sh"),'$RES',polymer_code)
       # Place the residue name in the input PDB file residue name columns
        with open(pdb_file, "rt") as fin:

          new_pdb_file = str(pdb_file+"temp")
          with open(new_pdb_file, "wt") as fout:
              for line in fin:
                  line_list = [char for char in line]
                  line_start = ''.join(line_list[0:6])
                  residue_code = ''.join(line_list[22:25])
                  if residue_code == '   ':
                    line_list[22:25] = polymer_code
                  if line_start == 'HETATM' or line_start == 'ATOM  ':
                    del line_list[26]
                  line = ''.join(line_list)
                  fout.write(line)
        os.rename(new_pdb_file,param_pdb)
        subprocess.run(["chmod","+x","param.sh"])
        subprocess.run(["./param.sh"])
       
        return

def minimize(polymer_code):
        """

        Parameters
        ----------

        """
        # Setup an energy minimization of the initial structure guess
        subprocess.run(["gmx","grompp","-f","em.mdp","-p","topol.top","-c","solvated.gro","-o","em3"])
        exit()
        # Run the energy minimization
        subprocess.run(["gmx","mdrun","-v","-deffnm","em3"])
        subprocess.run(["gmx","trjconv","-f","em3.trr","-s","solvated.gro","-o","em3.pdb"])
        return("em3.pdb")

def solvate(polymer_code,input_pdb,solvent_density=0.5):
        """

        Parameters
        ----------

        """
        # Build a simulation box and fill it with solvent
        subprocess.run(["gmx","solvate","-cp",input_pdb,"-cs","solvent.pdb","-p",str(str(polymer_code)+".top"),"-o","solvated.gro"])

        return

def equilibrate():
        """

        Parameters
        ----------

        """
        # Setup an equilibration run with Berendsen barostat
        subprocess.run(["gmx","grompp","-f","berendsen.mdp","-p","topol.top","-c","em3.gro","-o","berendsen"])
        # Run the equilibration
        subprocess.run(["gmx","mdrun","-v","-deffnm","berendsen"])
        subprocess.run(["gmx","trjconv","-f","berendsen.trr","-s","solvated.gro","-o","berendsen.pdb"])

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
