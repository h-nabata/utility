import os

# sdfファイルのパス
filepath = 'C:/Users/'+os.getlogin()+'/Downloads/'
input_file = 'Ethanol.sdf'

# xyz座標の取得
atomnum = 0; atomname = []; atomcoord = []
with open(input_file, encoding='utf-8') as f:  # 読み込みモードで開く
    while True:  # 読み込める行が無くなるまで繰り返す
        line = f.readline()  # ファイル f を1行ずつ読み込んでいく
        elem = line.split()  # 読み込んだ行を要素に分割する
        if len(elem) > 10:  # 要素数が10を超える行のみ処理する
            atomname.append((elem[3]))  # 4列目の要素は元素記号
            atomcoord.append([float(elem[0]), float(elem[1]), float(elem[2])])  # 1～3列目の要素はx,y,z座標
            atomnum += 1  # 原子数をカウント
        if not line:  # 読み込める行が無くなったら break
            break

# 取得データの確認
print(atomnum)
print(atomname)
print(atomcoord)

output_file = filepath + input_file.replace(".sdf", "") + ".xyz"
with open(output_file, mode='w', encoding='utf-8') as outf:  # 書き込みモードで開く
    outf.write(str(atomnum) + "\n" + input_file.replace(".sdf", "") + "\n")
    for i in range(atomnum):
        xc = '{:>17.12f}'.format(float(atomcoord[i][0]))  # フォーマット整形
        yc = '{:>17.12f}'.format(float(atomcoord[i][1]))
        zc = '{:>17.12f}'.format(float(atomcoord[i][2]))
        output_linestr = str(atomname[i]) + " " + xc + " " + yc + " " + zc
        print(output_linestr)  # 整形データの確認
        outf.write(output_linestr + "\n")  # 改行が必要

print("Done.")
