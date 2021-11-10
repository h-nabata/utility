'''
"ris2ACSformat.py"
@h-nabata  (2021/11/10 ver1.0)
A program to format the citation to the ACS style.
*input:  .ris files, etc.
*output: ACS style citation (text)
usage example: `py ris2ACSformat.py`
'''

with open('achs_jpcafh115_2877.ris', encoding='utf-8') as f:
    while True:  # 読み込める行が無くなるまで繰り返す
        line = f.readline()  # ファイル f を1行ずつ読み込んでいく
        elem = line.split()  # 読み込んだ行を要素に分割する
        if "A1  -" in line or "AU  -" in line:  # 著者
            print(elem[3][0]+". "+elem[2], sep=",", end=" ")
        if "JO  -" in line:  # journal名
            journal = line[6:-1]
            # for i in range(len(elem)-2):
            #     print(elem[i+2], end=" ")
        if "Y1  -" in line:  # 出版年
            year = elem[2]
        if "VL  -" in line:  # volume
            VL = elem[2]
        if "SP  -" in line:  # 最初のページ
            SP = elem[2]
        if "EP  -" in line:  # 最後のページ
            EP = elem[2]
        if not line:
            break

print(journal, year, VL, SP+"–"+EP+".", sep=", ", end=" ")
