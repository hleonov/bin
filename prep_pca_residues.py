import sys,os
import re
import pyutil
import optparse

def read_resi(rfile):
   lines=open(rfile).readlines()
   res_list=[]
   for line in lines:
      res_list = res_list + line.split()
   return res_list

#### Main ####
usage="Usage: %prog <res_list> <pdb_ref> <out_ndx> [options]"
parser = optparse.OptionParser(usage=usage)
#parser.add_option("-h","--help", action="help")
parser.add_option("-a", dest="atoms", default="bb", help="atoms to include (bb,ca,bb+cb,sc,all)")

(options, args) = parser.parse_args()
print options
print args
print options.atoms

atoms="bb"
if (options.atoms == "bb"):
   atoms="MainChain"
elif (options.atoms == "ca"):
   atoms="C-alpha"
elif (options.atoms == "bb+cb"):
   atoms="MainChain+Cb"
elif (options.atoms == "sc"):
   atoms = "SideChain"
elif (options.atoms == "all"):
   atoms = "System"

print atoms

r_list = read_resi(args[0])
sel_res = ""
sel_name = ""
for res in r_list:
   sel_res = sel_res + "r %d |" % int(res) 
   sel_name = sel_name + "r_%d_" % int(res)
sel_res = sel_res[0:-1]
sel_name = sel_name[0:29]
print sel_res
print sel_name
ndx_com = "printf \'%s\n\"%s\" & \"%s\" & \"Protein\"\nq\n\' | make_ndx -f %s -o %s" % (sel_res, sel_name, atoms, args[1], args[2])
pyutil.run_command(ndx_com)


