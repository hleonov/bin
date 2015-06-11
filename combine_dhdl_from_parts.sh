#!/bin/bash 

#combine the dhdl output for crooks phase runs with multiple parts

for lmb in 0 1
do
   cd lambda$lmb/morphes
   for fr in frame*
   #`seq 0 1 399`
   do 
#      echo $n
      cd $fr
      cat *.part*.xvg | grep -v "[#@&%]" | sort -nu > dhdl.xvg      
      cd ../
   done
   cd ../../
done
