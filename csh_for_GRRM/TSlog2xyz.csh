#!/bin/csh

set logs=`ls -l | grep -v EG | grep -E "LUP[0-9]_TS[0-9].log" | awk '{print $9}'`
set atomnum=5
set surfaceinfile="surfaceinfile.log"
echo ok
echo "" > allEG.xyz

@ i=1
while( $i <= $#logs )
  echo "" >> allEG.xyz 
  grep "INITIAL STRUCTURE" -A${atomnum} $logs[$i] >> allEG.xyz 
  sed -i "s/INITIAL STRUCTURE/${atomnum}\n[TS] $logs[$i]/g" allEG.xyz
  if( $#surfaceinfile > 0 ) then
    cat $surfaceinfile >> allEG.xyz
  endif
  echo "" >> allEG.xyz
  
  set optlinenum=`grep -e "Optimized" -e "exceed" -n $logs[$i] | cut -d ":" -f 1`
  set irclinenum=`grep -e "FORWARD" -e "BACKWARD" -n $logs[$i] | cut -d ":" -f 1`
  
  if( $#optlinenum > 0 ) then
    if( $optlinenum[1] < $irclinenum[2] ) then
      echo "${atomnum}\n[FW]" >> allEG.xyz
    else
      echo "${atomnum}\n[BW]" >> allEG.xyz
    endif
    @ tmpnum1 = $optlinenum[1] + 1
    @ tmpnum2 = $optlinenum[1] + ${atomnum}
    sed -n ${tmpnum1},${tmpnum2}p $logs[$i] >> allEG.xyz
    if( $#surfaceinfile > 0 ) then
      cat $surfaceinfile >> allEG.xyz
    endif
    echo "" >> allEG.xyz
  endif
  
  if( $#optlinenum > 1 ) then
    echo "${atomnum}\n[BW]" >> allEG.xyz
    @ tmpnum1 = $optlinenum[2] + 1
    @ tmpnum2 = $optlinenum[2] + ${atomnum}
    sed -n ${tmpnum1},${tmpnum2}p $logs[$i] >> allEG.xyz
    if( $#surfaceinfile > 0 ) then
      cat $surfaceinfile >> allEG.xyz
    endif
    echo "" >> allEG.xyz
  endif
  
  @ i++
end


