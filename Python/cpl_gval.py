#coding:utf-8
#2021-03-15 3:22 Yuuya Nagata
import sys
import os
import re
import numpy as np

#換算係数
electric_const = 254.35
magnetic_const = -0.97401

#ファイルの読み込み
args = sys.argv
with open(args[1]) as f:
    s = f.read()
P1 = s .rfind(' Ground to excited state transition electric dipole moments')
P2 = s .rfind(' Ground to excited state transition velocity dipole moments')
P3 = s .rfind(' Ground to excited state transition magnetic dipole moments')
P4 = s .rfind(' Ground to excited state transition velocity quadrupole moments')
P5 = s .rfind(' Excitation energies and oscillator strengths')
P6 = s .rfind(' SavETr')

electric = s[P1:P2]
magnetic = s[P3:P4]
excited = s[P5:P6]

#electric
electric = electric.replace(' Ground to excited state transition electric dipole moments (Au):\n','')
electric = electric.replace('       state          X           Y           Z        Dip. S.      Osc.\n         ','')
electric = electric.replace('\n         ','\n')
electric = electric.replace('\n        ','\n')
electric = electric.replace('\n       ','\n')
electric = re.sub(' +',',',electric)
electric_num = electric.count('\n')
electric = electric.rstrip('\n')
electric = electric.replace('\n',',')
electric_array = electric.split(',')
electric_array = np.array(electric_array)
electric_array = electric_array.reshape(electric_num, 6)
electric_array_state, electric_array_X, electric_array_Y, electric_array_Z, electric_arrayDipS, electric_arrayOsc = np.split(electric_array, [1,2,3,4,5], 1)
electric_array_X = electric_array_X.astype(np.longdouble) * electric_const
electric_array_Y = electric_array_Y.astype(np.longdouble) * electric_const
electric_array_Z = electric_array_Z.astype(np.longdouble) * electric_const
electric_array_XYZ = np.concatenate([electric_array_X, electric_array_Y, electric_array_Z], 1)

electric_L2 = ((electric_array_X ** 2) + (electric_array_Y ** 2) + (electric_array_Z ** 2)) ** 0.5

#magnetic
magnetic = magnetic.replace(' Ground to excited state transition magnetic dipole moments (Au):\n','')
magnetic = magnetic.replace('       state          X           Y           Z\n         ','')
magnetic = magnetic.replace('\n         ','\n')
magnetic = magnetic.replace('\n        ','\n')
magnetic = magnetic.replace('\n       ','\n')
magnetic = re.sub(' +',',',magnetic)
magnetic_num = magnetic.count('\n')
magnetic = magnetic.rstrip('\n')
magnetic = magnetic.replace('\n',',')
magnetic_array = magnetic.split(',')
magnetic_array = np.array(magnetic_array)
magnetic_array = magnetic_array.reshape(magnetic_num, 4)
magnetic_array_state, magnetic_array_X, magnetic_array_Y, magnetic_array_Z = np.split(magnetic_array, [1,2,3], 1)
magnetic_array_X = magnetic_array_X.astype(np.longdouble) * magnetic_const
magnetic_array_Y = magnetic_array_Y.astype(np.longdouble) * magnetic_const
magnetic_array_Z = magnetic_array_Z.astype(np.longdouble) * magnetic_const
magnetic_array_XYZ = np.concatenate([magnetic_array_X, magnetic_array_Y, magnetic_array_Z], 1)

magnetic_L2 = ((magnetic_array_X ** 2) + (magnetic_array_Y ** 2) + (magnetic_array_Z ** 2)) ** 0.5

#excited
excited = excited.split('\n')
excited = [line for line in excited if 'Excited State' in line]
excited = '\n'.join(excited)
excited = excited.replace(' Excited State   ','')
excited = excited.replace(':      ',',')
excited = excited.replace('      ',',')
excited = excited.replace(' eV  ',',')
excited = excited.replace(' nm  f=',',')
excited = excited.replace('  <S**2>=',',')
excited_num = excited.count('\n') + 1
excited = excited.rstrip('\n')
excited = excited.replace('\n',',')
excited_array = excited.split(',')
excited_array = np.array(excited_array)
excited_array = excited_array.reshape(excited_num, 6)
excited_array_state, excited_array_mode, excited_array_eV, excited_array_nm, excited_array_f, excited_array_S2 = np.split(excited_array, [1,2,3,4,5], 1)


#g値の計算用
cos_array = ((electric_array_X * magnetic_array_X) + (electric_array_Y * magnetic_array_Y) + (electric_array_Z * magnetic_array_Z)) / electric_L2 / magnetic_L2
g_array = 4 * electric_L2 * magnetic_L2 * cos_array / (electric_L2 * electric_L2 + magnetic_L2 * magnetic_L2)

#まとめ
all_array = np.concatenate([electric_array_state, electric_array_X, electric_array_Y, electric_array_Z, magnetic_array_X, magnetic_array_Y, magnetic_array_Z, electric_L2, magnetic_L2, cos_array, excited_array_eV, excited_array_nm, g_array, electric_arrayOsc], 1)
title_array = np.array(['state', 'electric X', 'electric Y', 'electric Z', 'magnetic X', 'magnetic Y', 'magnetic Z', 'mu', 'm', 'cos theta', 'eV', 'nm', 'g value', 'Oscillator Strength'])
all_array2 = np.vstack([title_array, all_array])

#デバッグ用
#print(all_array2)

#ファイルの書き出し
args[1] = args[1].replace('.log','')
filename = args[1] + '.csv'
np.savetxt(filename,all_array2,delimiter=',', fmt="%s")

