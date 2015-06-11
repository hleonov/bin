from ngmx import *
import operator
from string import upper
from pyutil import run_command
import sys, re

if (len(sys.argv) < 2):
   print "Usage: reorder_atp_atoms.py <ff> <pdb-file>"
   exit()
ff = sys.argv[1] #"charmm" #amber
pdb_file = sys.argv[2] # "/home/hleonov/Research/abc/topol_check/umb/pure_distance_RC/short_charmm_from_C3/3A_ref.pdb"


if (ff == "amber"):
    itp_file = "/home/hleonov/Data/topology/ATP_amber.itp"
elif (ff == "charmm"):
    itp_file = "/home/hleonov/Data/topology/atp_charmm36.itp"

ordered=[]
lines = open(itp_file, 'r').readlines()[8:52]
for line in lines:
    ordered.append(line.split()[4]) 
    
d = {}
for i,v in enumerate(ordered):
    d[v.upper()] = i    
d['Mg'] = d['MG']
sorted_d = sorted(d.items(), key=operator.itemgetter(1))
print sorted_d

#for i in sorted_d:
#    print i, sorted_d
#len(d.keys())



#pdb_file = "/home/hleonov/Research/abc/PDB_atps/update_Dec_2014/3a5m.sm.pdb"


#res = get_ipython().getoutput(u'egrep "^(HETATM|ATOM)" $pdb_file')
res, err = run_command("egrep \"^(HETATM|ATOM)\" "+pdb_file)
#print res
lines = res.split("\n")
new = {}

trans = {"O5'": 'O5*', "C5'": 'C5*',"C4'": 'C4*',"O4'": 'O4*', "C3'": 'C3*',"O3'": 'O3*', "C2'": 'C2*',"O2'": 'O2*', "C1'": 'C1*',"O1'": 'O1*', "H5'1" : 'H50', 'H51' : "H5'2" , "H5'2" : 'H51', "H61" : 'H60', "H62" : 'H61', "H4'" : 'H40', "H1'" : 'H10', "H8" : 'H80', "H2''" : 'H20',  "H3'" : 'H30', "H3T" : "H3'", "H3T" : "H3*","Mg" : 'MG', "1H5'" : 'H50', "H2'" : "H2*" }#, "2H5'" : 'H51' } #"MG" : 'Mg',

if (ff =="charmm"): 
    trans = dict(zip(trans.values(),trans))

print "trans: ",trans
for line in lines[:-1]:
    atom = line.split()[2]
    
    #replace resid with generic number - only needed when using pdb2gmx to have MG be a part of the ATP molecule rather than a new residue.
    line = line[0:22]+"1   "+line[27::]
    print "\nprocessing atom: %s ==>  " %atom,
    if (trans.has_key(str(atom))):
        #print atom
        if (atom.startswith("H")): 
           print "assigning to ", trans[atom], d[trans[atom]]
           new[d[trans[atom]]] = line
        elif (atom.startswith("M")):
           #rename the residue name 
           print "assigning to ", trans[atom], d[trans[atom]]
           new[d[trans[atom]]] = line[0:17]+"ATP"+line[20::]
        else: 
           print "assigning to ", trans[atom], d[trans[atom]]
           new[d[trans[atom]]] = line.replace(atom, trans[atom])
           
    elif(d.has_key(str(atom))):
        print "assigning to ", atom, d[atom]
        new[d[atom]] = line   
        if (atom.startswith("M")):
           print "assigning to ", atom, d[atom]
           new[d[atom]] = line[0:17]+"ATP"+line[20:-1]
    else: #exceptions
       #skip unreadable hydrogens
       if (atom.startswith("H")):
          print "Warning: failed to find a hydrogen atom: ", atom
       #specific case where oxygens are ambigious
       else:
          match = re.search("^(O..)", atom)
          if (match):
             print "assigning to ", match.group(1), d[match.group(1)]
             new[d[match.group(1)]] = line
       print "Error: failed to find atom: ", atom
    #print d[atom]
    #print atom, line, 
#print len(res)
fp = open(pdb_file+".ordered.pdb", 'w')
for k in sorted(new):
    print >>fp, new[k]
fp.close()


# In[ ]:



