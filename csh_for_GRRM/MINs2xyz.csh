#!/bin/csh

# usage: `./MINs_to_xyz.csh`
# 反応経路自動探索プログラム「GRRM」を用いたMIN計算の結果をまとめるCシェルスクリプトの例です。
# すべてのGRRMプログラムの出力ファイルがカレントディレクトリ内に存在することを前提としており、存在するすべてのMIN計算の結果を参照します。
# "allopt.xyz"に最適化後の構造を、"allopt.log" にエネルギーを出力します。
# 適宜改変してご利用ください。（周期系に対応済み／ONIOM計算の結果では未テスト）

set name=`grep '# MIN' *.com | awk '{print $1}' | cut -d '.' -f 1`
set i=1
echo -n "" > allopt.xyz; echo -n "" > allopt.log

while( $i <= $#name )
  @ count=0
  if ( -e $name[$i]_PARAM.rrm == 1 ) then
    @ j = 1
    while( $count <= 4 )
      if (`sed -n ${j}p $name[$i]_PARAM.rrm | grep "\-\-\-\-"` != "")  @ count = $count + 1
      @ j = $j + 1
      if ( $count == 5 )  set atominfo=`sed -n ${j}p $name[$i]_PARAM.rrm`
    end
    if ( `grep "Therm" $name[$i].log` != "" ) then
      if ( `grep "replacements" $name[$i].log` != "" ) then
        echo "$atominfo[3]\n$name[$i] (G =`grep -A15 "replacements" $name[$i].log | grep Free | cut -d "=" -f 2`)" >> allopt.xyz
        echo "# $name[$i]" >> allopt.log
        grep -A15 "replacements" $name[$i].log >> allopt.log
        echo "$name[$i]\t`grep -A15 "replacements" $name[$i].log | grep Free | cut -d "=" -f 2`"
      else 
        echo "$atominfo[3]\n$name[$i] (G =`grep -A15 "Therm" $name[$i].log | grep Free | cut -d "=" -f 2`)" >> allopt.xyz
        echo "# $name[$i]" >> allopt.log
        grep -A15 "Therm" $name[$i].log >> allopt.log
        echo "$name[$i]\t`grep -A15 "Therm" $name[$i].log | grep Free | cut -d "=" -f 2`"
      endif
    else
      @ tmp_energy_line = $atominfo[3] + 1
      set el_ene=`grep -A$tmp_energy_line "Optimized" $name[$i].log | grep ENERGY`
      if ( `grep "TV" $name[$i]_PARAM.rrm` != "" ) then
        @ totalatomnum = $atominfo[3] + $atominfo[5] - 3
      else
        @ totalatomnum = $atominfo[3] + $atominfo[5]
      endif
      echo "$totalatomnum\n$name[$i] ( $el_ene )" >> allopt.xyz
      echo "# $name[$i]\n$el_ene" >> allopt.log
      echo "$name[$i]\t$el_ene"
    endif
    echo "$atominfo[1] $atominfo[2]"
    if ( `grep "Optimized" $name[$i].log` != "" ) then
      grep "Optimized" -A$atominfo[3] $name[$i].log | grep -v "Op" >> allopt.xyz
      if ( `grep "Frozen" $name[$i].com` != "" ) then
        grep "Frozen" -A$atominfo[5] $name[$i].com | grep -v "Frozen" | grep -v "TV" >> allopt.xyz
      endif
    endif 
    echo "" >> allopt.xyz; echo "" >> allopt.log
  endif
  @ i = $i + 1
end




