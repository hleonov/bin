#!/bin/bash

echo "USAGE: run_fma.sh 1.TRAJ.XTC 2.NR_MIN_COMP 3.DATA.XVG 4.TOPOL.TPR 5.INDEX.NDX 6.FIT_GROUP"

xtc=$1
max_comp=$2
xvg_file=$3
tpr_file=$4
ndx_file=$5
gr=$6

#echo "$tpr_file"
#exit
for i in `seq 1 $max_comp`; do
   ~/Programs/g_fma -f $xtc -s $tpr_file -n $ndx_file \
      -y $xvg_file -fit $gr -analysis $gr -dim $i \
      -o FMAmodel_$i.xvg -v FMAvec_$i.trr -ew FMAvec_ew_$i.trr
      
      paste $xvg_file  FMAmodel_$i.xvg | \
      awk '{if ($3==1) print $2,$5}' | ~/Programs/pls/regression \
      | grep CORR | awk -v ind=$i '{print ind, $5}' >> reg_model.dat  

      paste  $xvg_file FMAmodel_$i.xvg | \
      awk '{if ($3==2) print $2,$5}' | ~/Programs/pls/regression \
      | grep CORR | awk -v ind=$i '{print ind, $5}' >> reg_crossVal.dat
done
