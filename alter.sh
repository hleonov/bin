#!/bin/bash

for i in `qstat | grep "A19" | grep "qw " | awk '{printf "%d\n", $1}'`; 
   do  
   qalter -pe  *_fast 48 $i
   done
