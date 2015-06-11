# This script removes double entries from a Gromacs index file

import sys, os, re

infile  = sys.argv[1]
outfile = sys.argv[2]
all_lines = open(infile).readlines()

out_fh = open(outfile, 'w')

seen = []
write = 1
total = 0
reduced = 0

for line in all_lines:
   match = re.search("\[\s(\S+)\s\]", line)   
   if (match is not None):
      total = total+1
      ndx_group = match.group(1)
      if (ndx_group not in seen):
         seen.append(ndx_group)
         write = 1
         reduced = reduced + 1
      else:
         write = 0
      if (write == 1):
         print >>out_fh, line,
   else:
      if (write == 1):
         print >>out_fh, line,

out_fh.close()
print "Total groups found in input: %s" % total
print "Unique groups found in input: %s" % reduced
      
