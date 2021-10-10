#!/bin/csh

# usage: `./log2xyz.csh xxx.log`
# 反応経路自動探索プログラム「GRRM」を用いたMIN計算の全ITRをxyz形式にするCシェルスクリプトです。
# GRRMプログラムの .com ファイルと .log ファイルを参照します。
# 実行すると 0Geom_${JOBNAME}.xyz というxyzファイルが生成します。
# 適宜改変してご利用ください。（周期系に対応済み／ONIOM計算の結果では未テスト）


if ( `grep . ${argv}` != "" ) then
  set JOBNAME = `echo ${argv} | cut -d "." -f 1`
else
  set JOBNAME = ${argv}
endif

set total_atomnum = `grep "#" -A1000 ${JOBNAME}.com | grep -v TV | tail -n +4 | grep "Options" -B1000 | grep -v "Options" | grep -v rozen | grep -c " "`
set move_atomnum = `grep "#" -A1000 ${JOBNAME}.com | grep -v TV  | tail -n +4 | grep "rozen" -B1000 | grep -v rozen | grep -c " "`
set frozen_atomnum = `grep "rozen" -A1000 ${JOBNAME}.com | grep -v TV | grep "Options" -B1000 | grep -v "Options" | grep -v rozen | grep -c " "`
set itr_num = `grep "#" -c ${JOBNAME}.log`
echo -n "" > 0Geom_${JOBNAME}.xyz

@ i = 0
while ($i < $itr_num)
  echo $total_atomnum >> 0Geom_${JOBNAME}.xyz
  if ($frozen_atomnum > 0) then
    grep -x "# ITR. $i" -A$move_atomnum ${JOBNAME}.log >> 0Geom_${JOBNAME}.xyz
    grep "rozen" -A$frozen_atomnum ${JOBNAME}.com | grep -v rozen >> 0Geom_${JOBNAME}.xyz
  else
    grep -x "# ITR. $i" -A$total_atomnum ${JOBNAME}.log >> 0Geom_${JOBNAME}.xyz
  endif
  echo "" >> 0Geom_${JOBNAME}.xyz
  @ i = $i + 1
end



