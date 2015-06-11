import sys, os
import re
from glob import glob
import pyutil


infile = sys.argv[1]

pyutil.run_command("mv %s %s" % (infile, "dbl_"+infile))

fp = open(infile, 'w')

lines = open("dbl_"+infile).readlines()

for line in lines:
   if (len(line)>16 and line[16] == "A"):
      line = line[0:16]+" "+line[17:]
   elif (len(line)>16 and re.search("[B-Z]", line[16])):
      line=""         
   
   if (line != ""):
      print >>fp, line,
