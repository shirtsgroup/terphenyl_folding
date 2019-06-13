import os, sys
import numpy as np
sys.path.insert(0,'/home/gmeek/software/coordinateTransform')

def Lattice_parameters_to_Crystal_matrix(lattice_parameters):
   """
   This function takes the lattice parameters and returns the crystal matrix

   **Required Inputs
   lattice_parameters = crystal lattice parameters as an array ([a,b,c,alpha,beta,gamma])
   """
   # Computing pieces of the crystal lattice matrix
   Vxx = lattice_parameters[0]
   Vxy = lattice_parameters[1]*np.cos(np.radians(lattice_parameters[5]))
   Vxz = lattice_parameters[2]*np.cos(np.radians(lattice_parameters[4]))

   Vyy = lattice_parameters[1]*np.sin(np.radians(lattice_parameters[5]))
   Vyz = lattice_parameters[2]*(np.cos(np.radians(lattice_parameters[3])) - np.cos(np.radians(lattice_parameters[4]))*np.cos(np.radians(lattice_parameters[5])))/np.sin(np.radians(lattice_parameters[5]))

   Vzz = np.sqrt(lattice_parameters[2]**2 - Vxz**2 - Vyz**2)
   # Combining the pieces of the matrix together
   if np.absolute(Vxy) < 1e-10:
       Vxy = 0.
   if np.absolute(Vxz) < 1e-10:
       Vxz = 0.
   if np.absolute(Vyz) < 1e-10:
       Vyz = 0.
   crystal_matrix = [[Vxx, Vxy, Vxz], [0., Vyy, Vyz], [0., 0., Vzz]]
   return crystal_matrix

def frac2cart(cellParams, fracCoords):
  cartCoords = []
  cellParam = Lattice_parameters_to_Crystal_matrix(cellParams)
  for i in fracCoords:
    xPos = i[0]*cellParam[0][0] + i[1]*cellParam[1][0] + i[2]*cellParam[2][0]
    yPos = i[0]*cellParam[0][1] + i[1]*cellParam[1][1] + i[2]*cellParam[2][1]
    zPos = i[0]*cellParam[0][2] + i[1]*cellParam[1][2] + i[2]*cellParam[2][2]
    cartCoords.append([xPos, yPos, zPos])
  return cartCoords


fractional_coordinates_file = open("fractional_coordinates.frac","r")
fractional_coordinates_lines = fractional_coordinates_file.readlines()
fractional_coordinates_file.close()
atom_name_list = []
for line in fractional_coordinates_lines:
  atom_name = line.split()[0]
  atom_name_list.append(atom_name)
# Fix the syntax for the fractional coordinates
fractional_coordinates = []
for line in fractional_coordinates_lines:
  coordinates = line.split()[1:4]
  fixed_coordinates = []
  for coor in coordinates:
    fixed_coordinates.append(float(coor.split('(')[0]))
  fractional_coordinates.append(fixed_coordinates)

# a, b, c
cell_params = [17.997,6.1519,19.636,90.0,106.549,90.0]
cartesian_coordinates = frac2cart(cell_params,fractional_coordinates)
output_file = open("cartesian_coordinates.pdb","w")
for atom_index in range(len(atom_name_list)):
  atom_name = atom_name_list[atom_index]
  coordinates = [round(cartesian_coordinates[atom_index][i],3) for i in range(3)]
  output_file.write(str("HETATM"+str("{:>5}".format(atom_index))+str("{:>3}".format(atom_name[0]))+str("{:>6}".format(str("   ")))+"  "+str("{:>4}".format(1))+"     "+str("{:>7}".format(coordinates[0]))+" "+str("{:>7}".format(coordinates[1]))+" "+str("{:>7}".format(coordinates[2]))+"  1.00  0.00\n"))
output_file.write("END")
output_file.close()
exit()

