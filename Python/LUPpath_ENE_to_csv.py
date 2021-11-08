'''
"LUPpath_ENE_to_csv.py"
@h-nabata  (2021/11/08 ver1.0)
A program to extract the values of LUP-path energy profile.
*input:  "LUPOUTt.log" files
*output: *csv file
usage: `py path_extract.py`
The JOBNAME list should be rewritten in need.
This code requires the pandas library.
'''

# coding: utf-8
import subprocess
import pandas as pd

JOBNAME = ["LUP1","LUP2","LUP3"]

for i in range(len(JOBNAME)):
    # ENERGY の値を grepコマンド で抽出して ENE.log に書き出す
    subprocess.call("grep ENE "+JOBNAME[i]+"_LUPOUTt.log > "+JOBNAME[i]+"_ENE.log", shell=True)
    tmpENE = []
    # ENE.log のファイルを開く
    with open(JOBNAME[i]+"_ENE.log", encoding='utf-8') as f:
        while True:  # 読み込める行が無くなるまで繰り返す
            line = f.readline()  # ファイル f を1行ずつ読み込んでいく
            elem = line.split()  # 読み込んだ行を要素に分割する
            if len(elem) > 3:  # 要素数が 3 を超える行のみ処理する
                tmpENE.append((elem[1]))
            if not line:
                break
    tmp_df = pd.DataFrame({JOBNAME[i]:tmpENE})
    if i == 0:
        df = tmp_df
    else:
        df = pd.concat([df, tmp_df], axis=1)

# 各JOBのLUPパスの ENERGY のデータを allENE.csv に列ごとに書き出す
df.to_csv('allENE.csv', index=None)
