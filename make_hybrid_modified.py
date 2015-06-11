import sys, os
from pymacs import *
from pymacs.forcefield import *
from pymacs.ndx import *

is_Set = 11

def do_log(fp, s):
    l = "make_hybrid__log_> "+s
    print >>sys.stderr, l
    print >>fp, l


def make_pairs(m1, m2, plist = None, grps = None):

    pairs = []
    if plist:
        for n1, n2 in plist:
            print n1,n2
            a1 = m1.fetch_atoms(n1)[0]
            a2 = m2.fetch_atoms(n2)[0]
            pairs.append( (a1, a2))
        return pairs
    if grps:
        lst1 = m1.get_by_id(grps[0])
        lst2 = m2.get_by_id(grps[1])
        for atom in lst1:
            mi = .05 # nm
            keep = None
            for at in lst2:
                d = atom-at
                if d < mi:
                    keep = at
                    mi = d
            if keep is not None:
                pairs.append( (atom, keep) )
        return pairs
        
    for atom in m1.atoms:
        mi = .05 # nm
        keep = None
        for at in m2.atoms:
            d = atom-at
            if d < mi:
                keep = at
                mi = d
        if keep is not None:
            pairs.append( (atom, keep) )
    return pairs

def assign_ff(model, itp):
    for i, atom in enumerate(model.atoms):
        at = itp.atoms[i]
        atom.atomtype = at.atomtype
        atom.cgnr = at.cgnr
        atom.q = at.q
        atom.m = at.m
        atom.atomtypeB = at.atomtypeB
        atom.qB = at.qB
        atom.mB = at.mB
        
def get_ff_entry(ids, itp, what = 'bond'):        
    if what == 'bond':
        for b in itp.bonds:
            if (b[0] == ids[0] and b[1] == ids[1]) or \
               (b[1] == ids[0] and b[0] == ids[1]):
                return b[3:]
    elif what == 'angle':
        for b in itp.angles:
            if (b[0] == ids[0] and b[1] == ids[1] and b[2] == ids[2]) or \
               (b[2] == ids[0] and b[1] == ids[1] and b[0] == ids[2]):
                return b[4:]
    elif what == 'dihedral':
        if (ids[4] == 3):
            for b in itp.dihedrals:
                if (b[0] == ids[0] and b[1] == ids[1] and b[2] == ids[2] and b[3] == ids[3]) or \
                   (b[3] == ids[0] and b[2] == ids[1] and b[1] == ids[2] and b[0] == ids[3]):
                    return b[5:]
        elif (ids[4] == 1):
            for b in itp.dihedrals:
                sum = 0
                for b1 in range(0,4):
                    for b2 in range(0,4):
                        if (b[b1] == b[b2]):
                            sum += 1
                            break
                if (sum == 4):
                    return b[5:]
    return None

def read_pairs_file(fn):
    l = open(fn).readlines()
    plst = []
    for line in l:
        plst.append(line.split())
    return plst

################################################################################33

desc=()

# define input/output files

files= [
    (efSTO, "-l1", "ligand1.pdb", ffREAD ),
    (efSTO, "-l2", "ligand2.pdb", ffREAD ),
    (efITP, "-itp1", "lig1.itp", ffREAD ),
    (efITP, "-itp2", "lig2.itp", ffREAD ),
    (efNDX, "-n1", "scaffold1", ffOPTRD ),
    (efNDX, "-n2", "scaffold2", ffOPTRD ),
    (efDAT, "-pairs", "pairs", ffOPTRD ),
    (efSTO, "-o", "merged.pdb", ffWRITE ),
    (efITP, "-oitp", "merged.itp", ffWRITE ),
    (efLOG, "-log", "hybrid.log", ffWRITE ),
    ]

# define options

options=[]

# pass options, files and the command line to pymacs

args=parse_common_args(sys.argv,files,options, desc)

if args['-pairs']['flag'] == 11:
    read_pairs_from_file = True
else:
    read_pairs_from_file = False

if args['-n1']['flag'] == 11:
    read_from_idx = True
else:
    read_from_idx = False

logfile = open(args['-log']['fns'],'w')

if read_from_idx and read_pairs_from_file:
    do_log(logfile, "Error: Can either read a pair list or scaffold index files!")
    do_log(logfile,"Exiting!")
    sys.exit(1)
    

do_log(logfile,'Reading ligand 1 from: "%s"' % args['-l1']['fns'])
do_log(logfile,'Reading ligand 2 from: "%s"' % args['-l2']['fns'])

m1 = Model().read(args['-l1'])
m2 = Model().read(args['-l2'])

do_log(logfile,'Reading itp file 1 from: "%s"' % args['-itp1']['fns'])
do_log(logfile,'Reading itp file 2 from: "%s"' % args['-itp2']['fns'])

itp1 = ITPFile(args['-itp1']['fns'])
itp2 = ITPFile(args['-itp2']['fns'])

do_log(logfile,"Assigning forcefield parameters....")
assign_ff(m1, itp1)
assign_ff(m2,itp2)
do_log(logfile,"Making pairs.....")

if read_pairs_from_file:
    do_log(logfile,'Reading file with atom pairs: "%s"' % args['-pairs']['fns'])  
    plst = read_pairs_file(args['-pairs']['fns'])
else:
    plst = None

if read_from_idx:
    do_log(logfile,'Reading scaffold index file: "%s"' % args['-n1']['fns'])  
    grp1 = IndexFile(args['-n1']['fns']).dic['scaffold']
    do_log(logfile,'Reading scaffold index file: "%s"' % args['-n2']['fns'])  
    grp2 = IndexFile(args['-n2']['fns']).dic['scaffold']
    # now we add all atoms with bonds to scaffold atoms
    for b in itp1.bonds:
        if b[0] in grp1 and b[1] not in grp1:
            grp1.append(b[1])
            do_log(logfile,'Adding atom %s to scaffold 1' % m1.atoms[b[1]-1].name)  
        elif b[1] in grp1 and b[0] not in grp1:
            grp1.append(b[0])
            do_log(logfile,'Adding atom %s to scaffold 1' % m1.atoms[b[0]-1].name)  
    for b in itp2.bonds:
        if b[0] in grp2 and b[1] not in grp2:
            grp2.append(b[1])
            do_log(logfile,'Adding atom %s to scaffold 2' % m2.atoms[b[1]-1].name)  
        elif b[1] in grp2 and b[0] not in grp2:
            grp2.append(b[0])
            do_log(logfile,'Adding atom %s to scaffold 2' % m2.atoms[b[0]-1].name)
    grps = [grp1, grp2]
else:
    grps = None



pairs = make_pairs(m1, m2, plst,grps)
morphsA = map(lambda p: p[1], pairs)
morphsB = map(lambda p: p[0], pairs)
dumsA = []
for atom in m2.atoms:
    if atom not in morphsA:
        dumsA.append(atom)
dumsB = []
for atom in m1.atoms:
    if atom not in morphsB:
        dumsB.append(atom)
do_log(logfile, "Generated %d atom-atom pairs" % len(pairs))
do_log(logfile,"Dummies in state A: %d" % len(dumsA))
do_log(logfile,"Dummies in state B: %d" % len(dumsB))



do_log(logfile,"Making B-states....")
for a1, a2 in pairs:
    a1.atomtypeB = a2.atomtype
    a1.nameB = a2.name
    a1.qB = a2.q
    a1.mB = a2.m
    a1.idB = a2.id
    do_log(logfile, "Atom....: %4d  %12s | %6.2f | %6.2f -> %12s | %6.2f | %6.2f" %\
           (a1.id, a1.atomtype, a1.q, a1.m, a1.atomtypeB, a1.qB, a1.mB))

for atom in dumsA:
    atom.id_old = atom.id
    atom.nameB = atom.name
    atom.name = 'D'+atom.name
    atom.atomtypeB = atom.atomtype
    atom.atomtype = 'DUM_'+atom.atomtype
    atom.qB = atom.q
    atom.q = 0
    atom.mB = atom.m
    atom.m = 1.
    m1.residues[0].append(atom)
    do_log(logfile, "Dummy...: %4d  %12s | %6.2f | %6.2f -> %12s | %6.2f | %6.2f" %\
           (atom.id, atom.atomtype, atom.q, atom.m, atom.atomtypeB, atom.qB, atom.mB))

for atom in dumsB:
    atom.atomtypeB = 'DUM_'+atom.atomtype
    atom.qB = 0
    atom.mB = 1.
    do_log(logfile, "Dummy...: %4d  %12s | %6.2f | %6.2f -> %12s | %6.2f | %6.2f" %\
           (atom.id, atom.atomtype, atom.q, atom.m, atom.atomtypeB, atom.qB, atom.mB))

id_dicAB = {}
id_dicBA = {}
for atom in m1.atoms:
    if hasattr(atom,"idB"):
        id_dicAB[atom.id] = atom.idB
        id_dicBA[atom.idB] = atom.id
    if hasattr(atom,"id_old"):
        id_dicAB[atom.id] = atom.id_old
        id_dicBA[atom.id_old] = atom.id
        
do_log(logfile, "Generating bonded parameters....")    
# go over bonds
newbonds = []

for b in itp1.bonds:
    id1 = b[0]
    id2 = b[1]
    a1 = m1.atoms[id1-1]
    a2 = m1.atoms[id2-1]
    bOk = False
    if hasattr(a1,"idB") and hasattr(a2,"idB"):
        idB1 = a1.idB
        idB2 = a2.idB
        entr = get_ff_entry([idB1, idB2], itp2, what= 'bond')
        if entr is not None:
            newbonds.append (b+entr)
            bOk = True
        else:
            bOk = False
    elif a1.atomtypeB[:3] == 'DUM' or a2.atomtypeB[:3] == 'DUM':
        entr = get_ff_entry([a1.id, a2.id], itp1, what= 'bond')
        if entr is not None:
            newbonds.append (b+entr)
            bOk = True
        else:
            bOk = False
    else:
        newbonds.append(b)
        bOk = True

    if not bOk:
        do_log(logfile, "Error: Something went wrong while assigning bonds!")
        do_log(logfile, "A-> Atom1: %d-%s Atom2: %d-%s" %(a1.id, a1.name, a2.id, a2.name))
        do_log(logfile, "B-> Atom1: %d-%s Atom2: %d-%s" %(a1.idB, a1.nameB, a2.idB, a2.nameB))
        do_log(logfile,"Exiting....")
        sys.exit(1)


# angles
newangles = []
for b in itp1.angles:
    id1 = b[0]
    id2 = b[1]
    id3 = b[2]
    a1 = m1.atoms[id1-1]
    a2 = m1.atoms[id2-1]
    a3 = m1.atoms[id3-1]
    bOk = False
    if hasattr(a1,"idB") and hasattr(a2,"idB") and hasattr(a3,"idB"):
        idB1 = a1.idB
        idB2 = a2.idB
        idB3 = a3.idB
        entr = get_ff_entry([idB1, idB2, idB3], itp2, what= 'angle')
        if entr is not None:
            newangles.append (b+entr)
            bOk = True
        else:
            bOk = False
    elif a1.atomtypeB[:3] == 'DUM' or \
             a2.atomtypeB[:3] == 'DUM' or \
             a3.atomtypeB[:3] == 'DUM':
        entr = get_ff_entry([a1.id, a2.id, a3.id], itp1, what= 'angle')
        if entr is not None:
            newangles.append (b+entr)
            bOk = True
        else:
            bOk = False
    else:
        newangles.append(b)
        bOk = True

    if not bOk:
        do_log(logfile, "Error: Something went wrong while assigning angles!")
        do_log(logfile, "A-> Atom1: %d-%s Atom2: %d-%s Atom3: %d-%s" \
               %(a1.id, a1.name, a2.id, a2.name, a3.id, a3.name))
        do_log(logfile, "B-> Atom1: %d-%s Atom2: %d-%s Atom3: %d-%s" \
               %(a1.idB, a1.nameB, a2.idB, a2.nameB, a3.idB, a3.nameB))
        do_log(logfile,"Exiting....")
        sys.exit(1)

# dihedrals
newdihedrals = []
for b in itp1.dihedrals:
    id1 = b[0]
    id2 = b[1]
    id3 = b[2]
    id4 = b[3]
    dih_type = b[4];
    a1 = m1.atoms[id1-1]
    a2 = m1.atoms[id2-1]
    a3 = m1.atoms[id3-1]
    a4 = m1.atoms[id4-1]
    bOk = False
    if hasattr(a1,"idB") and hasattr(a2,"idB") and \
           hasattr(a3,"idB") and hasattr(a4,"idB"):
        idB1 = a1.idB
        idB2 = a2.idB
        idB3 = a3.idB
        idB4 = a4.idB
        entr = get_ff_entry([idB1, idB2, idB3, idB4, dih_type], itp2, what= 'dihedral')
        if entr is not None:
            newdihedrals.append (b+entr)
            bOk = True
    elif a1.atomtypeB[:3] == 'DUM' or \
             a2.atomtypeB[:3] == 'DUM' or \
             a3.atomtypeB[:3] == 'DUM' or \
             a4.atomtypeB[:3] == 'DUM':

        entr = get_ff_entry([a1.id, a2.id, a3.id, a4.id, dih_type], itp1, what= 'dihedral')
        if entr is not None:
            newdihedrals.append (b+entr)
            bOk = True
        else:
            bOk = False
    else:
        newdihedrals.append(b)
        bOk = True
            
    if not bOk:
        do_log(logfile, "Error: Something went wrong while assigning dihedrals!")
        do_log(logfile, "A-> Atom1: %d-%s Atom2: %d-%s Atom3: %d-%s Atom3: %d-%s" \
               %(a1.id, a1.name, a2.id, a2.name, a3.id, a3.name, a4.id, a4.name))
        do_log(logfile, "B-> Atom1: %d-%s Atom2: %d-%s Atom3: %d-%s Atom3: %d-%s" \
               %(a1.idB, a1.nameB, a2.idB, a2.nameB, a3.idB, a3.nameB, a4.idB, a4.nameB))
        do_log(logfile,"Exiting....")
        sys.exit(1)

# now we all parameter for pairs
# let's go for the dummies

for b in itp2.bonds:
    newid1 = id_dicBA[b[0]]
    newid2 = id_dicBA[b[1]]
    a1 = m1.atoms[newid1-1]
    a2 = m1.atoms[newid2-1]
    if a1.atomtype.startswith('DUM') or \
       a2.atomtype.startswith('DUM'):
        newbonds.append( [newid1, newid2, 1, b[3], b[4], b[3], b[4]] )

        
for b in itp2.angles:
    newid1 = id_dicBA[b[0]]
    newid2 = id_dicBA[b[1]]
    newid3 = id_dicBA[b[2]]
    a1 = m1.atoms[newid1-1]
    a2 = m1.atoms[newid2-1]
    a3 = m1.atoms[newid3-1]
    if a1.atomtype.startswith('DUM') or \
       a2.atomtype.startswith('DUM') or \
       a3.atomtype.startswith('DUM'):
        newangles.append( [newid1, newid2, newid3, 1, b[4], b[5], b[4], b[5]] )

for b in itp2.dihedrals:
    newid1 = id_dicBA[b[0]]
    newid2 = id_dicBA[b[1]]
    newid3 = id_dicBA[b[2]]
    newid4 = id_dicBA[b[3]]
    a1 = m1.atoms[newid1-1]
    a2 = m1.atoms[newid2-1]
    a3 = m1.atoms[newid3-1]
    a4 = m1.atoms[newid4-1]
    if a1.atomtype.startswith('DUM') or \
       a2.atomtype.startswith('DUM') or \
       a3.atomtype.startswith('DUM') or \
       a4.atomtype.startswith('DUM'):
        newdihedrals.append( [newid1, newid2, newid3, newid4, b[4]] + b[5:] + b[5:] )

# make pairs
newpairs = []
pp = []
for p in itp1.pairs:
    newpairs.append( p )
    pp.append( (p[0],p[1]) )

for p in itp2.pairs:
    newid1 = id_dicBA[p[0]]
    newid2 = id_dicBA[p[1]]
    if (newid1, newid2) not in pp and \
       (newid2, newid1) not in pp:
        newpairs.append([ newid1, newid2, 1] )

do_log(logfile, "Generating new itp file")    
        
newitp = ITPFile()
newitp.atoms = m1.atoms
for i, atom in enumerate(newitp.atoms):
    atom.cgnr = i +1
newitp.bonds = newbonds
newitp.pairs = newpairs
newitp.angles = newangles
newitp.dihedrals = newdihedrals
do_log(logfile, 'Writing new itp file: "%s"' % args['-oitp']['fns'])    
newitp.write(args['-oitp']['fns'])
do_log(logfile, 'Writing new structure file: "%s"' % args['-o']['fns'])    
#m1.write(args['-o'])
do_log(logfile, 'Writing dummy forcefield file: "%s"' % ('ff'+args['-oitp']['fns']))    

fp = open('ff'+args['-oitp']['fns'],'w')
dd = []
print >>fp, '[ atomtypes ]'
for atom in m1.atoms:
    if atom.atomtype.startswith('DUM') and atom.atomtype not in dd:
        print >>fp, '%8s %12.6f %12.6f %3s %12.6f %12.6f' % \
              (atom.atomtype, 0, 0, 'A',0,0)
        dd.append(atom.atomtype)
    elif atom.atomtypeB.startswith('DUM') and atom.atomtypeB not in dd:
        print >>fp, '%8s %12.6f %12.6f %3s %12.6f %12.6f' % \
              (atom.atomtypeB, 0, 0, 'A',0,0)
        dd.append(atom.atomtypeB)
m1.write(args['-o'])
thanx()
