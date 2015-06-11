import sys,os
import re
import pyutil
import string 

if (len(sys.argv)<3):
   print "Usage: python %s infile outfile newName" %sys.argv[0]
   exit

infile = sys.argv[1]
outfile = sys.argv[2]
newName = sys.argv[3]

out = open(outfile, 'w')
itp_lines = open(infile).readlines()

directive = ""
for line in itp_lines:
   #find category
   if (line[:1] == "["):
      match = re.search("\[\s+(\S+)\s+\]", line)
      directive = match.group(1)
      print >>out, line,
   elif (line[:1] == ";" or line.strip() == "" or line[:1] == "#"):
      print >>out, line,
   #change molecule name so it does not collide with old topology
   elif (directive == "moleculetype"):
      match = re.search("(\S+)\s+(\d+)", line)
      if (match):
         molname = match.group(1)
         print >>out, newName+"   "+match.group(2)
         print "Changing molecule name from "+molname+" to "+newName
      else:
         print >>out, line,
   #switch B states of atoms
   elif (directive == "atoms"):
      cols = line.split()
      if (len(cols) > 8):         
         #switch atom names
         cols[1], cols[8] = cols[8], cols[1]
         #switch charges
         cols[6], cols[9] = cols[9], cols[6]
         #switch masses
         cols[7], cols[10] = cols[10], cols[7]
         new_cols = string.join(cols, "\t\t")
         print >>out, new_cols
      else:
         print "Warning: No B states for atom:"
         print line,
         print >>out, line,               

   elif (directive == "bonds"):
      cols = line.split()
      if (len(cols)>5):
         #switch bond param
         cols[3], cols[5] = cols[5], cols[3]
         cols[4], cols[6] = cols[6], cols[4]
         new_cols = string.join(cols, "\t\t")
         print >>out, new_cols
      else:
         print "Warning: No B states for bond:"
         print line,
         print >>out, line,
   elif (directive == "angles"):
      cols = line.split()
      if (len(cols)>6):
         #switch angle param
         cols[4], cols[6] = cols[6], cols[4]
         cols[5], cols[7] = cols[7], cols[5]
         new_cols = string.join(cols, "\t\t")
         print >>out, new_cols
      else:
         print "Warning: No B states for angle:"
         print line,
         print >>out, line,
# Note: this only takes care of dihedrals with c0..c5 parameters.
   elif (directive == "dihedrals"):
      line = line.rstrip()
      line = line.rstrip(';')
      print line
      cols = line.split()
      if (len(cols)>11):
         #switch dihedral param
         cols[5], cols[11] = cols[11], cols[5]
         cols[6], cols[12] = cols[12], cols[6]
         cols[7], cols[13] = cols[13], cols[7]
         cols[8], cols[14] = cols[14], cols[8]
         cols[9], cols[15] = cols[15], cols[9]
         cols[10], cols[16] = cols[16], cols[10]
         new_cols = string.join(cols, "\t\t")
         print >>out, new_cols+"\n"
      else:
         print "Warning: No B states for dihedral:"
         print line,
         print >>out, line
   else:
      print >>out, line,   
         
