import math

'''
dispersion correction for SIESTA input
made by H. Nabata (2022/04/08)
'''

# https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.20495
grimme_C6Parameter_list = [0.14, 0.08, 1.61, 1.61, 3.13, 1.75, 1.23, 0.7, 0.75, 0.63, 5.71, 5.71, 10.79, 9.23, 7.84, 5.57, 5.07, 4.61, 10.8, 10.8, 10.8, 10.8, 10.8, 10.8, 10.8, 10.8, 10.8, 10.8, 10.8, 10.8, 16.99, 17.1, 16.37, 12.64, 12.47, 12.01, 24.67, 24.67, 24.67, 24.67, 24.67, 24.67, 24.67, 24.67, 24.67, 24.67, 24.67, 24.67, 37.32, 38.71, 38.44, 31.74, 31.5, 29.99]
grimme_vdw_list = [1.001, 1.012, 0.825, 1.408, 1.485, 1.452, 1.397, 1.342, 1.287, 1.243, 1.144, 1.364, 1.639, 1.716, 1.705, 1.683, 1.639, 1.595, 1.485, 1.474, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.65, 1.727, 1.76, 1.771, 1.749, 1.727, 1.628, 1.606, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.672, 1.804, 1.881, 1.892, 1.892, 1.881]
element_list = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

## input (element names)
species_list = ["B", "N", "H", "O"]

species_index_list = []
for i in range(len(species_list)):
    species_index_list.append(element_list.index(species_list[i]))
# print(species_index_list)

for i in range(len(species_list)):
    for j in range(i+1):
        print(i+1, j+1, "Grimme", '{:7.5}'.format(math.sqrt(grimme_C6Parameter_list[species_index_list[i]]*grimme_C6Parameter_list[species_index_list[j]])*10.36416), '{:7.5}'.format(grimme_vdw_list[species_index_list[i]] + grimme_vdw_list[species_index_list[j]]) , "  # "+species_list[i]+"-"+species_list[j])   # C6Parameter
# 10.36416 is the scaling factor to convert from [J mol−1 nm6] to [eV Å6]

'''
### In this case, this code generates a set of parameters as below.
1 1 Grimme   32.44    2.97   # B-B
2 1 Grimme  20.336   2.882   # N-B
2 2 Grimme  12.748   2.794   # N-N
3 1 Grimme  6.8607   2.486   # H-B
3 2 Grimme  4.3008   2.398   # H-N
3 3 Grimme   1.451   2.002   # H-H
4 1 Grimme  15.341   2.827   # O-B
4 2 Grimme  9.6169   2.739   # O-N
4 3 Grimme  3.2445   2.343   # O-H
4 4 Grimme  7.2549   2.684   # O-O
'''
