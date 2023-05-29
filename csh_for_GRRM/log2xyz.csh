#!/bin/csh

# usage: `./log2xyz.csh xxx.log`
# 反応経路自動探索プログラム「GRRM」を用いたMIN計算の全ITRをxyz形式にするCシェルスクリプトです。
# GRRMプログラムの .com ファイルと .log ファイルを参照します。
# 実行すると 0Geom_${JOBNAME}.xyz というxyzファイルが生成します。
# 適宜改変してご利用ください。（周期系に対応済み／ONIOM計算の結果では未テスト／EQ_list、PT_listに対応済み）


if ( `echo ${argv} | grep .` != "" ) then
  set tmpJOBNAME = `echo ${argv} | cut -d "." -f 1`
else
  set tmpJOBNAME = ${argv}
endif

set listmode = 0
if ( `echo ${argv} | grep _LUPOUTt.log` != "" ) then
  set itr_string = "NODE"
else if ( `echo ${argv} | grep _EQ_list` != "" ) then
  set itr_string = "Geometry of EQ"
  set listmode = 1
else if ( `echo ${argv} | grep _PT_list` != "" ) then
  set itr_string = "Geometry of TS"
  set listmode = 2
else
  set itr_string = "ITR."
endif


if ( $listmode > 0 ) then
  set JOBNAME = `echo $tmpJOBNAME | sed "s/_EQ_list//" | sed "s/_PT_list//" | sed "s/_TS_list//"`
  echo "(listmode = on)"
else
  set JOBNAME = $tmpJOBNAME
endif
echo "JOBNAME = $JOBNAME"

set total_atomnum = `grep "#" -A5000 ${JOBNAME}.com | grep -v TV | tail -n +4 | grep "Options" -B5000 | grep -v "Options" | grep -v rozen | grep -c " "`
set move_atomnum = `grep "#" -A5000 ${JOBNAME}.com | grep -v TV  | tail -n +4 | grep "rozen" -B5000 | grep -v rozen | grep -c " "`
set frozen_atomnum = `grep "rozen" -A5000 ${JOBNAME}.com | grep -v TV | grep "Options" -B5000 | grep -v "Options" | grep -v rozen | grep -c " "`

if ( $listmode == 0 ) then
  set itr_num = `grep "#" -c ${JOBNAME}.log`
  echo -n "" > 0Geom_${JOBNAME}.xyz
  sed "s/: EQ Converged\!\!//g" ${JOBNAME}.log > ${JOBNAME}.log2
  sed -i "s/: PP Converged\!\!//g" ${JOBNAME}.log2
else if ( $listmode == 1 ) then
  set itr_num = `grep "#" -c ${JOBNAME}_EQ_list.log`
  echo -n "" > 0Geom_${JOBNAME}_EQ_list.xyz
  cat ${JOBNAME}_EQ_list.log > ${JOBNAME}_EQ_list.log2
else if ( $listmode == 1 ) then
  set itr_num = `grep "#" -c ${JOBNAME}_PT_list.log`
  echo -n "" > 0Geom_${JOBNAME}_PT_list.xyz
  cat ${JOBNAME}_PT_list.log > ${JOBNAME}_PT_list.log2
endif
echo "$total_atomnum $move_atomnum $frozen_atomnum $itr_num"

@ i = 0
while ($i < $itr_num)
  if ($listmode == 1) then
    echo $total_atomnum >> 0Geom_${JOBNAME}_EQ_list.xyz
    if ($frozen_atomnum > 0) then
      grep "# ${itr_string} ${i}," -A$move_atomnum ${JOBNAME}_EQ_list.log2 >> 0Geom_${JOBNAME}_EQ_list.xyz
      grep "rozen" -A$frozen_atomnum ${JOBNAME}.com | grep -v rozen >> 0Geom_${JOBNAME}_EQ_list.xyz
    else
      grep "# ${itr_string} ${i}," -A$move_atomnum ${JOBNAME}_EQ_list.log2 >> 0Geom_${JOBNAME}_EQ_list.xyz
    endif
  else if ($listmode == 2) then
    echo $total_atomnum >> 0Geom_${JOBNAME}_PT_list.xyz
    if ($frozen_atomnum > 0) then
      grep "# ${itr_string} ${i}," -A$move_atomnum ${JOBNAME}_PT_list.log2 >> 0Geom_${JOBNAME}_PT_list.xyz
      grep "rozen" -A$frozen_atomnum ${JOBNAME}.com | grep -v rozen >> 0Geom_${JOBNAME}_PT_list.xyz
    else
      grep "# ${itr_string} ${i}," -A$move_atomnum ${JOBNAME}_PT_list.log2 >> 0Geom_${JOBNAME}_PT_list.xyz
    endif
  else if ($listmode == 0) then
    echo $total_atomnum >> 0Geom_${JOBNAME}.xyz
    if ($frozen_atomnum > 0) then
      grep -x "# ${itr_string} $i" -A$move_atomnum ${JOBNAME}.log2 >> 0Geom_${JOBNAME}.xyz
      grep "rozen" -A$frozen_atomnum ${JOBNAME}.com | grep -v rozen >> 0Geom_${JOBNAME}.xyz
    else
      grep -x "# ${itr_string} $i" -A$move_atomnum ${JOBNAME}.log2 >> 0Geom_${JOBNAME}.xyz
    endif
  endif
  echo "" >> 0Geom_${JOBNAME}.xyz
  @ i = $i + 1
end




