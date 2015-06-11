import pymacs, sys
from pymacs.builder import *

mid_res = sys.argv[1]

if (mid_res == None):
   print "Usage: python make_tripeptide.py <1-letter resname>"
   exit()
   
gpg = build_chain("G%sG" % mid_res ,hydrogens = True)
gpg.add_nterm_cap()
gpg.add_cterm_cap()
gpg.write("g%sg.pdb" % mid_res) 
