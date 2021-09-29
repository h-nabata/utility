'''
"molecule_rotation.py"
@h-nabata  (2021/09/30 ver1.0)
A program that randomly rotates the coordinates of molecules around the center of gravity.
*input:  xyz-format
*output: xyz-format
'''

import os
import math
import random

# 回転させる分子のxyzファイルのパス（複数指定可）
xyz_files = ["C:/Users/"+str(os.getlogin())+"/Downloads/tryptophan.xyz",
             "C:/Users/"+str(os.getlogin())+"/Downloads/tyrosine.xyz"]
# 生成する構造数
generate_num = 20


def prod_dot(list1, list2):  # ベクトルの内積
    if len(list1) == len(list2):
        return sum([list1[i] * list2[i] for i in range(len(list1))])
    else:
        return print("ERROR: The number of elements in the list does not match.")
        sys.exit()

def prod_3Dcross(list1, list2):  # 3次元ベクトルの外積
    if len(list1) == len(list2):
        list1a = list1 + list1[0:2]  # コードの簡単化のために配列を複製して追加
        list2a = list2 + list2[0:2]
        return [list1a[i+1]*list2a[i+2]-list1a[i+2]*list2a[i+1] for i in range(len(list1))]
    else:
        return print("ERROR: The number of elements in the list does not match.")
        sys.exit()

def Elementname(a):  # 元素記号
    if a == 1:
        return "H"
    if a == 2:
        return "He"
    if a == 3:
        return "Li"
    if a == 4:
        return "Be"
    if a == 5:
        return "B"
    if a == 6:
        return "C"
    if a == 7:
        return "N"
    if a == 8:
        return "O"
    if a == 9:
        return "F"
    if a == 10:
        return "Ne"
    if a == 11:
        return "Na"
    if a == 12:
        return "Mg"
    if a == 13:
        return "Al"
    if a == 14:
        return "Si"
    if a == 15:
        return "P"
    if a == 16:
        return "S"
    if a == 17:
        return "Cl"
    if a == 18:
        return "Ar"
    if a == 19:
        return "K"
    if a == 20:
        return "Ca"
    if a == 21:
        return "Sc"
    if a == 22:
        return "Ti"
    if a == 23:
        return "V"
    if a == 24:
        return "Cr"
    if a == 25:
        return "Mn"
    if a == 26:
        return "Fe"
    if a == 27:
        return "Co"
    if a == 28:
        return "Ni"
    if a == 29:
        return "Cu"
    if a == 30:
        return "Zn"
    if a == 31:
        return "Ga"
    if a == 32:
        return "Ge"
    if a == 33:
        return "As"
    if a == 34:
        return "Se"
    if a == 35:
        return "Br"
    if a == 36:
        return "Kr"
    if a == 37:
        return "Rb"
    if a == 38:
        return "Sr"
    if a == 39:
        return "Y"
    if a == 40:
        return "Zr"
    if a == 41:
        return "Nb"
    if a == 42:
        return "Mo"
    if a == 43:
        return "Tc"
    if a == 44:
        return "Ru"
    if a == 45:
        return "Rh"
    if a == 46:
        return "Pd"
    if a == 47:
        return "Ag"
    if a == 48:
        return "Cd"
    if a == 49:
        return "In"
    if a == 50:
        return "Sn"
    if a == 51:
        return "Sb"
    if a == 52:
        return "Te"
    if a == 53:
        return "I"
    if a == 54:
        return "Xe"
    if a == 55:
        return "Cs"
    if a == 56:
        return "Ba"
    if a == 57:
        return "La"
    if a == 58:
        return "Ce"
    if a == 59:
        return "Pr"
    if a == 60:
        return "Nd"
    if a == 61:
        return "Pm"
    if a == 62:
        return "Sm"
    if a == 63:
        return "Eu"
    if a == 64:
        return "Gd"
    if a == 65:
        return "Tb"
    if a == 66:
        return "Dy"
    if a == 67:
        return "Ho"
    if a == 68:
        return "Er"
    if a == 69:
        return "Tm"
    if a == 70:
        return "Yb"
    if a == 71:
        return "Lu"
    if a == 72:
        return "Hf"
    if a == 73:
        return "Ta"
    if a == 74:
        return "W"
    if a == 75:
        return "Re"
    if a == 76:
        return "Os"
    if a == 77:
        return "Ir"
    if a == 78:
        return "Pt"
    if a == 79:
        return "Au"
    if a == 80:
        return "Hg"
    if a == 81:
        return "Tl"
    if a == 82:
        return "Pb"
    if a == 83:
        return "Bi"
    if a == 84:
        return "Po"
    if a == 85:
        return "At"
    if a == 86:
        return "Rn"
    if a == 87:
        return "Fr"
    if a == 88:
        return "Ra"
    if a == 89:
        return "Ac"
    if a == 90:
        return "Th"
    if a == 91:
        return "Pa"
    if a == 92:
        return "U"
    if a == 93:
        return "Np"
    if a == 94:
        return "Pu"
    if a == 95:
        return "Am"
    if a == 96:
        return "Cm"
    if a == 97:
        return "Bk"
    if a == 98:
        return "Cf"
    if a == 99:
        return "Es"
    if a == 100:
        return "Fm"
    if a == 101:
        return "Md"
    if a == 102:
        return "No"
    if a == 103:
        return "Lr"
    if a == 104:
        return "Rf"
    if a == 105:
        return "Db"
    if a == 106:
        return "Sg"
    if a == 107:
        return "Bh"
    if a == 108:
        return "Hs"
    if a == 109:
        return "Mt"
    if a == 110:
        return "Ds"
    if a == 111:
        return "Rg"
    if a == 112:
        return "Cn"
    if a == 113:
        return "Nh"
    if a == 114:
        return "Fl"
    if a == 115:
        return "Mc"
    if a == 116:
        return "Lv"
    if a == 117:
        return "Ts"
    if a == 118:
        return "Og"

def Atomicmass(a):  # 原子量
    if a == "H":
        return 1.008
    if a == "He":
        return 4.003
    if a == "Li":
        return 6.941
    if a == "Be":
        return 9.012
    if a == "B":
        return 10.81
    if a == "C":
        return 12.01
    if a == "N":
        return 14.01
    if a == "O":
        return 16.00
    if a == "F":
        return 19.00
    if a == "Ne":
        return 20.18
    if a == "Na":
        return 22.99
    if a == "Mg":
        return 24.31
    if a == "Al":
        return 26.98
    if a == "Si":
        return 28.09
    if a == "P":
        return 30.97
    if a == "S":
        return 32.07
    if a == "Cl":
        return 35.45
    if a == "Ar":
        return 39.95
    if a == "K":
        return 39.10
    if a == "Ca":
        return 40.08
    if a == "Sc":
        return 44.96
    if a == "Ti":
        return 47.88
    if a == "V":
        return 50.94
    if a == "Cr":
        return 52.00
    if a == "Mn":
        return 54.94
    if a == "Fe":
        return 55.85
    if a == "Co":
        return 58.93
    if a == "Ni":
        return 58.69
    if a == "Cu":
        return 63.55
    if a == "Zn":
        return 65.39
    if a == "Ga":
        return 69.72
    if a == "Ge":
        return 72.61
    if a == "As":
        return 74.92
    if a == "Se":
        return 78.96
    if a == "Br":
        return 79.90
    if a == "Kr":
        return 83.80
    if a == "Rb":
        return 85.47
    if a == "Sr":
        return 87.62
    if a == "Y":
        return 88.91
    if a == "Zr":
        return 91.22
    if a == "Nb":
        return 92.91
    if a == "Mo":
        return 95.94
    if a == "Tc":
        return 99.00
    if a == "Ru":
        return 101.10
    if a == "Rh":
        return 102.90
    if a == "Pd":
        return 106.40
    if a == "Ag":
        return 107.90
    if a == "Cd":
        return 112.40
    if a == "In":
        return 114.80
    if a == "Sn":
        return 118.70
    if a == "Sb":
        return 121.80
    if a == "Te":
        return 127.60
    if a == "I":
        return 126.90
    if a == "Xe":
        return 131.30
    if a == "Cs":
        return 132.90
    if a == "Ba":
        return 137.30
    if a == "La":
        return 138.90
    if a == "Ce":
        return 140.10
    if a == "Pr":
        return 140.90
    if a == "Nd":
        return 144.20
    if a == "Pm":
        return 145.00
    if a == "Sm":
        return 150.40
    if a == "Eu":
        return 152.00
    if a == "Gd":
        return 157.30
    if a == "Tb":
        return 158.90
    if a == "Dy":
        return 162.50
    if a == "Ho":
        return 164.90
    if a == "Er":
        return 167.30
    if a == "Tm":
        return 168.90
    if a == "Yb":
        return 173.00
    if a == "Lu":
        return 175.00
    if a == "Hf":
        return 178.50
    if a == "Ta":
        return 180.90
    if a == "W":
        return 183.80
    if a == "Re":
        return 186.20
    if a == "Os":
        return 190.20
    if a == "Ir":
        return 192.20
    if a == "Pt":
        return 195.10
    if a == "Au":
        return 197.00
    if a == "Hg":
        return 200.60
    if a == "Tl":
        return 204.40
    if a == "Pb":
        return 207.20
    if a == "Bi":
        return 209.00
    if a == "Po":
        return 210.00
    if a == "At":
        return 210.00
    if a == "Rn":
        return 222.00
    if a == "Fr":
        return 223.00
    if a == "Ra":
        return 226.00
    if a == "Ac":
        return 227.00
    if a == "Th":
        return 232.00
    if a == "Pa":
        return 231.00
    if a == "U":
        return 238.00
    if a == "Np":
        return 237.00
    if a == "Pu":
        return 239.00
    if a == "Am":
        return 243.00
    if a == "Cm":
        return 247.00
    if a == "Bk":
        return 247.00
    if a == "Cf":
        return 252.00
    if a == "Es":
        return 252.00
    if a == "Fm":
        return 257.00
    if a == "Md":
        return 256.00
    if a == "No":
        return 259.00
    if a == "Lr":
        return 260.00
    if a == "TV":
        return 0.000
    return 0.000


mol_atomnum = []; mol_atomname = []; mol_atomcoord = []  # 取得データを格納するリスト

# xyz座標の取得
for i in range(len(xyz_files)):
    atomnum_count = 0
    tmp_atomname = []
    tmp_atomcoord = []
    with open(xyz_files[i], encoding='utf-8') as f:
        while True:
            line = f.readline()
            elem = line.split()
            # print(elem)
            if len(elem) == 4:
                if str.isdigit(elem[0]):
                    tmp_atomname.append(Elementname(int(elem[0])))  # 1列目の要素が原子番号の場合は元素記号に変換
                else:
                    tmp_atomname.append(elem[0])  # 1列目の要素が元素記号の場合はそのまま格納
                tmp_atomcoord.append([float(elem[1]), float(elem[2]), float(elem[3])])
                atomnum_count += 1
            # print(line.split())
            if not line:
                mol_atomnum.append(atomnum_count)
                mol_atomname.append(tmp_atomname)
                mol_atomcoord.append(tmp_atomcoord)
                break

# 取得データの確認（オプショナル）
for i in range(len(xyz_files)):
    print(mol_atomnum[i])
    print(mol_atomname[i])
    print(mol_atomcoord[i])

# 重心の計算
GravityCenter_coord = []  # 平行移動前の重心の座標は GravityCenter_coord に格納
for i in range(len(xyz_files)):
    tmp_GravityCenter_coord = []
    for j in range(3):
        tmp_GravityCenter_coord.append(sum([mol_atomcoord[i][k][j] * Atomicmass(mol_atomname[i][k]) for k in range(mol_atomnum[i])]) / sum([Atomicmass(mol_atomname[i][k]) for k in range(mol_atomnum[i])]))  # jは原子の番号
    GravityCenter_coord.append(tmp_GravityCenter_coord)
# print("GravityCenter_coord", GravityCenter_coord)  # 重心の座標の確認

# 重心を原点に合わせるように分子を平行移動
mol_normatomcoord = []  # 平行移動後のxyz座標は mol_normatomcoord に格納
for i in range(len(xyz_files)):
    tmp_mol_normatomcoord = []
    for j in range(mol_atomnum[i]):
        tmp_mol_normatomcoord.append([mol_atomcoord[i][j][0] - GravityCenter_coord[i][0], mol_atomcoord[i][j][1] - GravityCenter_coord[i][1], mol_atomcoord[i][j][2] - GravityCenter_coord[i][2]])
    mol_normatomcoord.append(tmp_mol_normatomcoord)

# 平行移動後のxyz出力（オプショナル）
out_xyz_files = [str(xyz_files[i]).replace(".xyz", "")
                 + "_GCset.xyz" for i in range(len(xyz_files))]
# print(out_xyz_files)  # 出力パスの確認
for i in range(len(xyz_files)):
    with open(out_xyz_files[i], mode='w', encoding='utf-8') as outf:
        outf.write(str(mol_atomnum[i]) + "\n")
        outf.write(str(os.path.basename(out_xyz_files[i]).replace(".xyz", "")) + "\n")
        for j in range(mol_atomnum[i]):
            output_linestr = str(mol_atomname[i][j]) + " " + str(mol_normatomcoord[i][j][0]) + " " + str(mol_normatomcoord[i][j][1]) + " " + str(mol_normatomcoord[i][j][2]) + "\n"
            outf.write(output_linestr)

## 重心（今は原点）まわりの3D回転操作
# # 原点からの距離
# DistancefromOrigin = []  # 各原子の原点からの距離は DistancefromOrigin に格納
# for i in range(len(xyz_files)):
#     tmp_DistancefromOrigin = []
#     for j in range(mol_atomnum[i]):
#         tmp_DistancefromOrigin.append(math.sqrt(
#             mol_normatomcoord[i][j][0] ** 2 + mol_normatomcoord[i][j][1] ** 2 + mol_normatomcoord[i][j][2] ** 2))  # jは原子の番号
#     DistancefromOrigin.append(tmp_DistancefromOrigin)
# # print(DistancefromOrigin)

# 回転操作（x,y,z軸周りでそれぞれ回転操作を施す）
for tmp_generate_num in range(generate_num):
    theta = [2*math.pi * random.random(), 2*math.pi * random.random(), 2*math.pi * random.random()]
    rotate_x = [[1.0, 0.0, 0.0], [0.0, math.cos(theta[0]), -math.sin(theta[0])], [0.0, math.sin(theta[0]), math.cos(theta[0])]]
    rotate_y = [[math.cos(theta[1]), 0.0, math.sin(theta[1])], [0.0, 1.0, 0.0], [-math.sin(theta[1]), 0.0, math.cos(theta[1])]]
    rotate_z = [[math.cos(theta[2]), -math.sin(theta[2]), 0.0], [math.sin(theta[2]), math.cos(theta[2]), 0.0], [0.0, 0.0, 1.0]]
    mol_rotatedcoord = []
    for i in range(len(xyz_files)):
        tmp_mol_rotatedcoord = []
        for j in range(mol_atomnum[i]):
            tmp_mol_normatomcoord = [mol_normatomcoord[i][j][0], mol_normatomcoord[i][j][1], mol_normatomcoord[i][j][2]]
            tmp_tmp_mol_rotatedcoord = [tmp_mol_normatomcoord[0] * rotate_x[k][0] + tmp_mol_normatomcoord[1] * rotate_x[k][1] + tmp_mol_normatomcoord[2] * rotate_x[k][2] for k in range(3)]  # x軸周りに回転
            tmp_tmp_mol_rotatedcoord = [tmp_tmp_mol_rotatedcoord[0] * rotate_y[k][0] + tmp_tmp_mol_rotatedcoord[1] * rotate_y[k][1] + tmp_tmp_mol_rotatedcoord[2] * rotate_y[k][2] for k in range(3)]  # y軸周りに回転
            tmp_tmp_mol_rotatedcoord = [tmp_tmp_mol_rotatedcoord[0] * rotate_z[k][0] + tmp_tmp_mol_rotatedcoord[1] * rotate_z[k][1] + tmp_tmp_mol_rotatedcoord[2] * rotate_z[k][2] for k in range(3)]  # z軸周りに回転
            tmp_mol_rotatedcoord.append(tmp_tmp_mol_rotatedcoord)
        mol_rotatedcoord.append(tmp_mol_rotatedcoord)
    # print(mol_rotatedcoord)
    # print(mol_normatomcoord)

    # 回転後のxyz出力（オプショナル）
    out_xyz_files = [str(xyz_files[i]).replace(".xyz", "") + "_rotated.xyz" for i in range(len(xyz_files))]
    # print(out_xyz_files)  # 出力パスの確認
    for i in range(len(xyz_files)):
        with open(out_xyz_files[i], mode='a', encoding='utf-8') as outf:
            outf.write(str(mol_atomnum[i]) + "\n")
            outf.write(str(os.path.basename(out_xyz_files[i]).replace(".xyz", "")) + "_" + str(tmp_generate_num) + "\n")
            for j in range(mol_atomnum[i]):
                output_linestr = str(mol_atomname[i][j]) + " " + str(mol_rotatedcoord[i][j][0]) + " " + str(mol_rotatedcoord[i][j][1]) + " " + str(mol_rotatedcoord[i][j][2]) + "\n"
                outf.write(output_linestr)
