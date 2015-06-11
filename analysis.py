import sys, os
import optparse
import pyutil
from glob import glob
import re
import string

#trjconv
def trjconv(infile, outfile, begin=0, tpr=None, ndx=None ,special="", printf=None):
   #basic command (-f, -o) + begin + special
   trj_com = "trjconv -f %s -o %s -b %s %s"  %(infile, outfile, begin, special)
   #add param
   if (tpr is not None):
      trj_com += " -s %s" %tpr
   if (ndx is not None):
      trj_com += " -n %s" %ndx
   if (printf is not None):
      trj_com = "printf \"%s\" | " %(printf) + trj_com
   print "Executing trjconv: %s" %trj_com
   pyutil.run_command(trj_com) 


def rmsd(infile, tpr):
   command = "printf \"3 3\n\" | g_rms -f %s -o rmsd_md.xvg -s %s" %(infile, tpr)
   pyutil.run_command(command)
   
# Concatenates part of the trajectory in a given directory 
# Meaning md.part000x.xtc/edr/trr
def concat_parts(base):
   files = glob("%s.part*.xtc" %base)   
   print files
   xtc_files = ""
   for f in files:
      xtc_files = xtc_files + f + " "     
   edr_files = xtc_files.replace("xtc", "edr")
   if (base == "traj"):
      edr_files = edr_files.replace(base, "ener")
   cat_com = "trjcat -f "+xtc_files+" -o %s.xtc" %base
   print cat_com
   if (not os.path.exists(base+".xtc")):
      pyutil.run_command(cat_com) 
   cat_com2 = cat_com.replace("xtc", "trr")
   if (not os.path.exists(base+".trr")):
      pyutil.run_command(cat_com2) 
      
   ene_com = "eneconv -f "+edr_files+" -o %s.edr" %opt.base_file
   print ene_com
   if (not os.path.exists(base+".edr")):
      pyutil.run_command(ene_com)

#deletes part000x files with the prefix $base
def rm_parts(base):
   try:
      pyutil.run_command("rm %s.part*.*" %base)
   except:
      print "Nothing to delete"



# **********************************************************************************

usage = "Usage: %prog [options]"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-f", dest="base_file", 
                  help="Simulations base name for trajectory to concat/analyze")
parser.add_option("--concat", dest="concat", action="store_true", 
                  help="Concatenate parts")
parser.add_option("--rm_parts", dest="rm_parts", action="store_true", 
                  help="Remove 'part' files")
parser.add_option("--pbc",dest="pbc", action="store_true", 
                  help="Run trjconv -pbc mol -ur compact")
parser.add_option("--tpr", dest="tpr", metavar="FILE", default="topol.tpr",
                  help="Reference .tpr file")
parser.add_option("--rmsd", dest="rmsd", action="store_true", 
                  help="Compute rmsd from trajectory")                 
(opt, args) = parser.parse_args()

print opt, args

if (opt.concat):
   concat_parts(opt.base_file)

if (opt.rm_parts):
   rm_parts(opt.base_file)

if (opt.pbc):      
   trjconv(opt.base_file, "md_pbc_mol.xtc", tpr=opt.tpr ,special="-pbc mol -ur compact", printf="0\n")
if (opt.rmsd):
   rmsd(opt.base_file, tpr=opt.tpr)
