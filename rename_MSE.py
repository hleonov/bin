import sys, os
import re
from glob import glob
import pyutil


infile = sys.argv[1]

pyutil.run_command("mv %s %s" % (infile, "mse_"+infile))

fp = open(infile, 'w')

lines = open("mse_"+infile).readlines()

for line in lines:
   match=re.search("^HETATM(.+)MSE", line)
   if (match):
      line = line.replace("HETATM", "ATOM  ")
   match2=re.search("^ATOM(.+)MSE", line)      
   if (match2): 
      line = line.replace("MSE ", "MET ")            
      if ("SE" in match2.group(0)):
         line = line.replace("SE", "SD")
         
   print >>fp, line,
