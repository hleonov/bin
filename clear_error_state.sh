#!/bin/bash
module load sge
for i in `qstat | grep 'Eqw' | awk '{print $1}'`
do
   qmod -c $i
done
