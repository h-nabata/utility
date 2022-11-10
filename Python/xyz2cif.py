import math
import os


print("Convert xyz cordinates into cif format. (by xyz2cif.py v1.0)")
filepath = 'C:/Users/'+os.getlogin()+'/Downloads/'
filename_list = ["sample_crystal.txt", "ZnO.xyz"]

# 3D coordinates (in angs.) like below can be converted to cif format (TVs need to be specified)
# Si    1.367182000000      4.101546000000      1.367182000000
# Si    0.000000000000      0.000000000000      2.734364000000
# Si    1.367182000000      1.367182000000      4.101546000000
# Si    0.000000000000      2.734364000000      0.000000000000
# Si    4.101546000000      4.101546000000      4.101546000000
# Si    2.734364000000      0.000000000000      0.000000000000
# Si    4.101546000000      1.367182000000      1.367182000000
# Si    2.734364000000      2.734364000000      2.734364000000
# TV    5.468730000000      0.000000000000      0.000000000000
# TV    0.000000000000      5.468730000000      0.000000000000
# TV    0.000000000000      0.000000000000      5.468730000000

for i in range(len(filename_list)):
    target_file = filepath + filename_list[i]
    with open(target_file, encoding='utf-8') as f:
        tmp_atomnum = 0;  tmp_nodenum = 0;  atomcount = 0
        atomnamelist=[]
        coordlist = []
        TVlist = []
        TVcount = 0
        while line := f.readline():
            elem = line.split()
            if len(elem) >= 4:
                if elem[0] == "TV":
                    TVlist.append([float(elem[1]), float(elem[2]), float(elem[3])])
                    TVcount += 1
                    if TVcount == 4:
                        print("Warning: only three TVs are allowed.")
                else:
                    atomnamelist.append(str(elem[0]))
                    coordlist.append([float(elem[1]), float(elem[2]), float(elem[3])])
            if not line:
                break

    # angström
    a_length = TVlist[0][0]
    b_length = math.sqrt(TVlist[1][0]**2.0 + TVlist[1][1]**2.0)
    c_length = math.sqrt(TVlist[2][0]**2.0 + TVlist[2][1]**2.0 + TVlist[2][2]**2.0)
    # radian
    alpha = math.acos(sum([i*j for (i, j) in zip(TVlist[1], TVlist[2])])/(b_length * c_length)) # b-c
    beta  = math.acos(sum([i*j for (i, j) in zip(TVlist[2], TVlist[0])])/(c_length * a_length)) # c-a
    gamma = math.acos(sum([i*j for (i, j) in zip(TVlist[0], TVlist[1])])/(a_length * b_length)) # a-b
    # degree
    degree_alpha = math.acos(sum([i*j for (i, j) in zip(TVlist[1], TVlist[2])])/(b_length * c_length)) * (180/math.pi) # b-c
    degree_beta  = math.acos(sum([i*j for (i, j) in zip(TVlist[2], TVlist[0])])/(c_length * a_length)) * (180/math.pi) # c-a
    degree_gamma = math.acos(sum([i*j for (i, j) in zip(TVlist[0], TVlist[1])])/(a_length * b_length)) * (180/math.pi) # a-b

    print("\n\nTV coordinates to cell parameters < " + filename_list[i])   # TVs as a lower triangular matrix
    print('{:>11.8f}'.format(a_length), '{:>11.8f}'.format(b_length), '{:>11.8f}'.format(c_length), sep=", ")
    print('{:>11.8f}'.format(degree_alpha), '{:>11.8f}'.format(degree_beta), '{:>11.8f}'.format(degree_gamma), sep=", ")
    print("Surface area (ang^2) =", '{:>14.10f}'.format(TVlist[0][0] * TVlist[1][1]))
    print("Cell volume  (ang^3) =", '{:>14.10f}'.format(TVlist[0][0] * TVlist[1][1] * TVlist[2][2]))

    # create a cif file
    outputall_file = filepath + filename_list[i] + ".cif"
    with open(outputall_file, mode='w', encoding='utf-8') as outf:
        outf.write("# generated using xyz2cif.py (by HN)\n" + filename_list[i] + "\n")
        outf.write("_symmetry_space_group_name_H-M   \'P 1\'\n_chemical_name_common ?\n_chemical_formula_sum ?\n")
        outf.write("_cell_length_a   "+'{:>11.8f}'.format(a_length) + "\n")
        outf.write("_cell_length_b   "+'{:>11.8f}'.format(b_length) + "\n")
        outf.write("_cell_length_c   "+'{:>11.8f}'.format(c_length) + "\n")
        outf.write("_cell_angle_alpha   "+'{:>11.8f}'.format(degree_alpha) + "\n")
        outf.write("_cell_angle_beta    "+'{:>11.8f}'.format(degree_beta) + "\n")
        outf.write("_cell_angle_gamma   "+'{:>11.8f}'.format(degree_gamma) + "\n")
        outf.write("_symmetry_Int_Tables_number   1\n_cell_volume   " + '{:>11.8f}'.format(TVlist[0][0] * TVlist[1][1] * TVlist[2][2]) + "\n_cell_formula_units_Z   " + str(len(coordlist)) + "\n")
        outf.write("loop_\n _symmetry_equiv_pos_site_id\n _symmetry_equiv_pos_as_xyz\n  1  'x, y, z'\n")
        outf.write("loop_\n _atom_site_type_symbol\n _atom_site_label\n _atom_site_symmetry_multiplicity\n _atom_site_fract_x\n _atom_site_fract_y\n _atom_site_fract_z\n _atom_site_occupancy\n")
        for j in range(len(atomnamelist)):
            v = math.sqrt(1 - math.cos(alpha)**2 - math.cos(beta)**2 - math.cos(gamma)**2 + 2*math.cos(alpha)*math.cos(beta)*math.cos(gamma))
            a_frac = coordlist[j][0] - coordlist[j][1] / math.tan(gamma) + coordlist[j][2] * ((math.cos(alpha)*math.cos(gamma)-math.cos(beta))/(v*math.sin(gamma)))
            b_frac = coordlist[j][1] / math.sin(gamma) + coordlist[j][2] * ((math.cos(beta)*math.cos(gamma)-math.cos(alpha))/(v*math.sin(gamma)))
            c_frac = coordlist[j][2] * math.sin(gamma) / v
            outf.write('  {:>2s}'.format(atomnamelist[j]))
            outf.write('  {:>2s}'.format(atomnamelist[j]) + str(j) + "  1")
            outf.write('  {:>9.8f}'.format(a_frac / a_length))  # x componet only
            outf.write('  {:>9.8f}'.format(b_frac / b_length))  # y componet
            outf.write('  {:>9.8f}'.format(c_frac / c_length))  # z componet
            outf.write("  1\n")
