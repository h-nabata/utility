'''
"ris2ACSformat.py"
@h-nabata  (2021/11/10 ver1.0, 2022/09/30 ver2.0)
A program to format the citation to the ACS style.
*input:  .ris files
*output: ACS or Nature style citation format (text)
usage example: `py ris2ACSformat.py`
'''

import os

form_num = 2  # 0: ACS, 1: Nature, 2: Nature (TeX format)
addEP = 0

with open('C:/Users/'+os.getlogin()+'/Downloads/acs.jpcc.0c11107.ris', encoding='utf-8') as f:
    journal = ""
    year = ""
    volume = ""
    SP = ""
    EP = ""
    EPflag = 0
    while True:  # 読み込める行が無くなるまで繰り返す
        line = f.readline()  # ファイル f を1行ずつ読み込んでいく
        elem = line.split()  # 読み込んだ行を要素に分割する
        if "A1  -" in line or "AU  -" in line and len(elem) > 2:  # 著者
            if "," in line:
                line = line.replace(',',', ').replace('  ',' ')
                elem = line.split()
            if ". " in line:
                print(elem[2],elem[3], sep=" ", end=", ")
            else:
                print(elem[-1][0]+". "+elem[2], sep=",", end=" ")
        if "JO  -" in line or "JA  -" in line and len(elem) > 2:  # journal名
            journal = line[6:-1]
        if "Y1  -" in line and len(elem) > 2:  # 出版年
            year = elem[2]
        if "PY  -" in line and len(elem) > 2:  # 出版年
            year = elem[2]
        if "VL  -" in line and len(elem) > 2:  # volume
            volume = elem[2]
        if "SP  -" in line and len(elem) > 2:  # 最初のページ
            SP = elem[2]
        if "EP  -" in line and len(elem) > 2:  # 最後のページ
            EP = elem[2]
            EPflag = 1
        if not line:
            break

if form_num == 0:
    if EPflag == 1 and addEP == 1:  # 最後のページの記載があれば
        print(journal, year, volume, SP+"–"+EP+".", sep=", ", end=" ")
    else:
        print(journal, year, volume, SP+".", sep=", ", end=" ")
elif form_num == 1:
    if EPflag == 1 and addEP == 1:  # 最後のページの記載があれば
        print(journal, volume, SP+"–"+EP+" ("+year+").", sep=", ", end=" ")
    else:
        print(journal, volume, SP+" ("+year+").", sep=", ", end=" ")
elif form_num == 2:
    if EPflag == 1 and addEP == 1:  # 最後のページの記載があれば
        print("\\textit{"+journal+"}", "{\\bf"+volume+"}", SP+"–"+EP+" ("+year+").", sep=", ", end=" ")
    else:
        print("\\textit{"+journal+"}", "{\\bf"+volume+"}", SP+" ("+year+").", sep=", ", end=" ")
