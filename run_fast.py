import sys, os
from glob import glob
import re
import optparse
import pyutil

def find_limits(base, dirs, beg, end):
   for d in dirs:
      match = re.search("%s(\d+)" % base, d)
      if (match):
         num=int(match.group(1))
         if (beg is ""):
            beg=num
         if (end is ""):
          end=num
         if (num<beg):
          beg=num
         elif (num>end):
          end=num        
   return beg, end
   
usage="Usage: %prog [options] <dir_basename>"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-b", dest="b", metavar="int", default="", help="Start from specified run number")
parser.add_option("-e", dest="e", metavar="int", default="", help="End at specified run number")
parser.add_option("--mdp", dest="mdp", metavar="FILE", help="mdp file for fast transition")
parser.add_option("--top", dest="top", metavar="FILE", help="topology file for fast transition")
parser.add_option("--job", dest="job", metavar="FILE", help="job file for fast transition")

#parser.add_option("--pdb", dest="pdb", metavar="FILE", help="initial pdb file for fast transition")
(options, args) = parser.parse_args()
print options
print args

dirs = glob("%s*" %args[0])
print dirs
if (options.b is "" or options.e is ""):
   (b,e) = find_limits(args[0], dirs, "","")
if (options.b is ""): 
   options.b=b
if (options.e is ""):
   options.e=e
for d in dirs: 
   match = re.search("%s(\d+)" % args[0], d)
   if match:
      run=int(match.group(1))
      if (options.b <= run and options.e >= run):
         print "In Run number %s " %run
         os.chdir(d)
         pyutil.run_command("grompp -f ../%s -c %s.pdb  -p ../%s -maxwarn 1" %(options.mdp, d, options.top))
         pyutil.run_command("qsub ../%s %s" %(options.job, d))
         os.chdir("..")
