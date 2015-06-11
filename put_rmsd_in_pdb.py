import os, sys
import re
import string

#reads and rmsd list for residues and puts the rmsd in the pdb file.
#there is one rmsd number for an entire residue. it is put in the B-factor 
#coloumn of the entire residue.

def read_rmsd(infile):
   rms_dic = {}
   rmsd_lines = open(infile).readlines()
   for line in rmsd_lines:
      entr = line.split()
      rms_dic[int(entr[0])] = float(entr[1])
   return rms_dic

def put_rmsd(pdb, rms_dic, outfile):
   out = open(outfile, 'w')
   pdb_lines = open(pdb).readlines()
   for line in pdb_lines:
      if (line[:4] not in ['ATOM']):
         print >>out, line,
      else:
         resi = int(line[22:26])
         if resi in rms_dic:
            print >>out, line[:61],"%2.3f" %rms_dic[resi]
         else:
            print >>out, line,
   
   
pdb_file = sys.argv[1]
rmsd_file = sys.argv[2]
outfile  = sys.argv[3]

rms_dic = read_rmsd(rmsd_file)
put_rmsd(pdb_file, rms_dic, outfile)


