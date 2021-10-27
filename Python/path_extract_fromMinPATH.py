'''
"path_extract_fromMinPATH.py"
@h-nabata  (2021/10/27 ver1.0)
A program to extract the numbers of all PTs and all TSs from MinPATH.rrm_* .
*input:  "MinPATH.rrm" files
*output: *list
usage: `py path_extract.py [num]`
([num] is an optional argument for the "head command")
If you have less than three MinPATH.rrm files, please rewrite the contents.
'''

# coding: utf-8
import sys
import subprocess

args = sys.argv
if len(args) > 2:
    print("ERROR: Too many args...")
    sys.exit()
elif len(args) == 2:
    print(args)
    # 文字列で取得してsplitで要素に分割
    cmd_MinPATH_0 = 'head *MinPATH.rrm_0 -n'+str(args[1])+'| grep "\->"'
    list_MinPATH_0 = subprocess.check_output(cmd_MinPATH_0, shell=True).split()
    cmd_MinPATH_1 = 'head *MinPATH.rrm_1 -n'+str(args[1])+'| grep "\->"'
    list_MinPATH_1 = subprocess.check_output(cmd_MinPATH_1, shell=True).split()
    cmd_MinPATH_2 = 'head *MinPATH.rrm_2 -n'+str(args[1])+'| grep "\->"'
    list_MinPATH_2 = subprocess.check_output(cmd_MinPATH_2, shell=True).split()
else:
    # 文字列で取得してsplitで要素に分割
    list_MinPATH_0 = subprocess.check_output('head *MinPATH.rrm_0 | grep "\->"', shell=True).split()
    list_MinPATH_1 = subprocess.check_output('head *MinPATH.rrm_1 | grep "\->"', shell=True).split()
    list_MinPATH_2 = subprocess.check_output('head *MinPATH.rrm_2 | grep "\->"', shell=True).split()


PT_list_for_refine_MinPATH_0 = []
TS_list_for_refine_MinPATH_0 = []

for i in range(len(list_MinPATH_0)):
    # print(list_MinPATH_0[i].decode('ascii')) # byte型を文字列にデコード
    if list_MinPATH_0[i].decode('ascii') == "PT":
        if not int(list_MinPATH_0[i+1].decode('ascii')) in PT_list_for_refine_MinPATH_0:
            PT_list_for_refine_MinPATH_0.append(int(list_MinPATH_0[i+1].decode('ascii')))
    if list_MinPATH_0[i].decode('ascii') == "TS":
        if not int(list_MinPATH_0[i+1].decode('ascii')) in TS_list_for_refine_MinPATH_0:
            TS_list_for_refine_MinPATH_0.append(int(list_MinPATH_0[i+1].decode('ascii')))

# # MinPATH_0 の結果
# PT_list_for_refine_MinPATH_0.sort()
# TS_list_for_refine_MinPATH_0.sort()
# print("Result : MinPATH.rrm_0")
# print(PT_list_for_refine_MinPATH_0)
# print(TS_list_for_refine_MinPATH_0)

# sys.exit()

PT_list_for_refine_MinPATH_1 = []
TS_list_for_refine_MinPATH_1 = []

for i in range(len(list_MinPATH_1)):
    # print(list_MinPATH_1[i].decode('ascii')) # byte型を文字列にデコード
    if list_MinPATH_1[i].decode('ascii') == "PT":
        if not int(list_MinPATH_1[i+1].decode('ascii')) in PT_list_for_refine_MinPATH_1:
            PT_list_for_refine_MinPATH_1.append(int(list_MinPATH_1[i+1].decode('ascii')))
    if list_MinPATH_1[i].decode('ascii') == "TS":
        if not int(list_MinPATH_1[i+1].decode('ascii')) in TS_list_for_refine_MinPATH_1:
            TS_list_for_refine_MinPATH_1.append(int(list_MinPATH_1[i+1].decode('ascii')))

# # MinPATH_1 の結果
# PT_list_for_refine_MinPATH_1.sort()
# TS_list_for_refine_MinPATH_1.sort()
# print("Result : MinPATH.rrm_1")
# print(PT_list_for_refine_MinPATH_1)
# print(TS_list_for_refine_MinPATH_1)


PT_list_for_refine_MinPATH_2 = []
TS_list_for_refine_MinPATH_2 = []

for i in range(len(list_MinPATH_2)):
    # print(list_MinPATH_2[i].decode('ascii')) # byte型を文字列にデコード
    if list_MinPATH_2[i].decode('ascii') == "PT":
        if not int(list_MinPATH_2[i+1].decode('ascii')) in PT_list_for_refine_MinPATH_2:
            PT_list_for_refine_MinPATH_2.append(int(list_MinPATH_2[i+1].decode('ascii')))
    if list_MinPATH_2[i].decode('ascii') == "TS":
        if not int(list_MinPATH_2[i+1].decode('ascii')) in TS_list_for_refine_MinPATH_2:
            TS_list_for_refine_MinPATH_2.append(int(list_MinPATH_2[i+1].decode('ascii')))

# # MinPATH_2 の結果
# PT_list_for_refine_MinPATH_2.sort()
# TS_list_for_refine_MinPATH_2.sort()
# print("Result : MinPATH.rrm_1")
# print(PT_list_for_refine_MinPATH_2)
# print(TS_list_for_refine_MinPATH_2)


# total の結果
for i in range(len(PT_list_for_refine_MinPATH_1)):
    if not PT_list_for_refine_MinPATH_1[i] in PT_list_for_refine_MinPATH_0:
        PT_list_for_refine_MinPATH_0.append(PT_list_for_refine_MinPATH_1[i])

for i in range(len(TS_list_for_refine_MinPATH_1)):
    if not TS_list_for_refine_MinPATH_1[i] in TS_list_for_refine_MinPATH_0:
        TS_list_for_refine_MinPATH_0.append(TS_list_for_refine_MinPATH_1[i])

for i in range(len(PT_list_for_refine_MinPATH_2)):
    if not PT_list_for_refine_MinPATH_2[i] in PT_list_for_refine_MinPATH_0:
        PT_list_for_refine_MinPATH_0.append(PT_list_for_refine_MinPATH_2[i])

for i in range(len(TS_list_for_refine_MinPATH_2)):
    if not TS_list_for_refine_MinPATH_2[i] in TS_list_for_refine_MinPATH_0:
        TS_list_for_refine_MinPATH_0.append(TS_list_for_refine_MinPATH_2[i])

PT_list_for_refine_MinPATH_0.sort()
TS_list_for_refine_MinPATH_0.sort()
print("Result : total")
print(PT_list_for_refine_MinPATH_0)
print(TS_list_for_refine_MinPATH_0)
