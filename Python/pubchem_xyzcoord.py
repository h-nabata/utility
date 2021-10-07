'''
"pubchem_xyzcoord.py"
@h-nabata  (2021/10/07 ver1.0)
A program that obtains the 3D structures of molecules from the compound database "PubChem" and converts them into xyz format.
PubChem:    https://pubchem.ncbi.nlm.nih.gov/
*input:     Compound CID
*output:    xyz-format
'''

import requests
import time
import os

# CIDを指定してsdfファイルの内容を取得（ここだけ入力すればOK）
CID = [10000 + i for i in range(10)]

for i in range(len(CID)):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/" + str(CID[i]) + "/record/SDF/?record_type=3d&response_type=display"
    r = requests.get(url)
    tmp_input_file = "tmp_3DCoordFile.txt"

    with open(tmp_input_file, encoding='utf-8', mode="w") as f:
        f.write(r.text)  # 取得データの書き出し

    # xyz座標の取得
    atomnum = 0; atomname = []; atomcoord = []
    with open(tmp_input_file, encoding='utf-8') as f:
        while True:  # 読み込める行が無くなるまで繰り返す
            line = f.readline()  # ファイル f を1行ずつ読み込んでいく
            elem = line.split()  # 読み込んだ行を要素に分割する
            if len(elem) > 10:  # 要素数が10を超える行のみ処理する
                atomname.append((elem[3]))
                atomcoord.append([float(elem[0]), float(elem[1]), float(elem[2])])
                atomnum += 1
            if not line:
                break

    # 取得データの確認（オプション）
    print(atomnum)
    print(atomname)
    print(atomcoord)

    # xyz形式でファイルに出力
    output_file = str(CID[i]) + ".xyz"
    with open(output_file, mode='w', encoding='utf-8') as outf:
        outf.write(str(atomnum) + "\n" + str(CID[i]) + "\n")
        for i in range(atomnum):
            output_linestr = str(atomname[i]) + " " + str(atomcoord[i][0]) + " " + str(atomcoord[i][1]) + " " + str(atomcoord[i][2]) + "\n"
            outf.write(output_linestr)

#     # 全構造をxyz形式でファイルに出力（オプション）
#     outputall_file = "C:/Users/" + str(os.getlogin()) + "/Downloads/mol_xyz/all.xyz"  # パスはご自由に指定して下さい
#     with open(outputall_file, mode='a', encoding='utf-8') as outf2:
#         outf2.write(str(atomnum) + "\n" + str(CID[i]) + "\n")
#         for j in range(atomnum):
#             output_linestr2 = str(atomname[j]) + " " + str(atomcoord[j][0]) + " " + str(atomcoord[j][1]) + " " + str(atomcoord[j][2]) + "\n"
#             outf2.write(output_linestr2)
#         outf2.write("\n")
            
    time.sleep(0.5)  # 0.5秒sleepする（オプション）
