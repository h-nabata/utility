import math

# constants
num_atoms = 20  # the number of atoms
angle_between_atoms = 360 / num_atoms  # Angles between adjacent atoms
distance_between_atoms = 1.6  # Å, Distance between adjacent atoms

# Calculate radius (using equilateral triangle formula)
radius = distance_between_atoms / (2 * math.sin(math.radians(angle_between_atoms / 2)))

# Calculate XYZ coordinates of each atom
atoms_coordinates = []
for i in range(num_atoms):
    angle = math.radians(angle_between_atoms * i)  # Convert angles to radians
    x = radius * math.cos(angle)
    y = radius * math.sin(angle)
    z = 0  # Assumed to be on the XY plane
    atoms_coordinates.append((x, y, z))

# Display coordinates
for i in range(len(atoms_coordinates)):
    print('Ar  ','{:>14.10f}'.format(atoms_coordinates[i][0]), '{:>14.10f}'.format(atoms_coordinates[i][1]), '{:>14.10f}'.format(atoms_coordinates[i][2]))
    # print('Ar  ','{:>14.10f}'.format(atoms_coordinates[i][0]), '{:>14.10f}'.format(atoms_coordinates[i][1]), '{:>14.10f}'.format(atoms_coordinates[i][2]), str(i+1))

# for i in range(len(atoms_coordinates)-1):
#     print("AddKeepPotential = "+str(i+1)+","+str(i+2)+" = 1.0 ; "+str(i+1)+" "+str(i+2)+" 1.6")
# print("AddKeepPotential = 1,"+str(len(atoms_coordinates))+" = 1.0 ; 1 "+str(len(atoms_coordinates))+" 1.6")
