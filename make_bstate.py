#!/usr/bin/python
#############################################################
#                                                           #
#           THIS FILE IS PART OF PYMACS                     #
#   WRITTEN BY DANIEL SEELIGER (dseelig@gwdg.de)            #
#                   2006-2008                               #
#    COMPUTATIONAL BIOMOLECULAR DYNAMICS GROUP              #
#  MAX-PLANCK-INSTITUTE FOR BIOPHYSICAL CHEMISTRY           #
#                   GOETTINGEN                              #
#                                                           #
#############################################################
__doc__="""
Program to add B-states in gromacs topology files
for free energy simulations.
"""

import sys,os
from pymacs import *
from pymacs.parser import *
from pymacs import forcefield
from pymacs import mutdb
from pymacs import cpp
import copy, time

NUC_MUTS = ['DAT','DAC','DAG','DCT','DCG','DCA',
            'DTA','DTG','DTC','DGA','DGC','DGT',
            'RAU','RAC','RAG','RCU','RCG','RCA',
            'RUA','RUG','RUC','RGA','RGC','RGU',
            ]

def do_log(fp, s):
    l = 'make_Bstate__log_> %s' % s
    print >>sys.stderr, l
    print >>fp, l
    



res_Q = {'K':1,
         'R':1,
         'Z':1,
         'D':-1,
         'E':-1
         }


def sum_charge_of_states(rlist):
    qA = []
    qB = []
    for r in rlist:
        qa = 0
        qb = 0
        for atom in r.atoms:
            qa+=atom.q
            if atoms_morphe([atom]):
                qb+=atom.qB
            else:
                qb+=atom.q
        qA.append(qa)
        qB.append(qb)
    return qA, qB

class ForceFieldError(Exception):
    def __init__(self, s):
        self.s = s
    def __str__(self):
        return repr(self.s)

def is_perturbed_residue(r):
    if atoms_morphe(r.atoms): return True
    return False

def last_perturbed_atom(logfile,r):

    max_order = 0
    last_atom = None
    for atom in r.atoms:
        if atoms_morphe([atom]) and atom.name not in ['N','CA','C','O','H']:
            if not atom.atomtype.startswith('DUM') and not atom.atomtypeB.startswith('DUM'):
                last_atom = atom
##             if atom.order > max_order and atom.symbol!='H':
##                 last_atom = atom
##                 max_order = atom.order
##     if last_atom == None:
##         for atom in r.atoms:
##             if atoms_morphe([atom]):
##                 if atom.order > max_order:
##                     last_atom = atom
##                     max_order = atom.order
    if last_atom == None:
        do_log(logfile,'Error: Could not find a perturbed atom to put rest charges on !')
        sys.exit()
    return last_atom

def write_atoms(fp,logfile,model,charges='AB', dummy_chargeA = 'on',\
                dummy_chargeB = 'on', vdw='AB', scale_mass=True, tqB = 0.,\
                full_morphe = False):
    global QA
    global QB
    QA = 0
    QB = 0
    
    for r in model.residues:
        if is_perturbed_residue(r):
            target_chargeB = tqB.pop(0)
            do_log(logfile,'Making target charge %g for residue %s' % (round(target_chargeB,5), r.resname))
            for atom in r.atoms:
                if atoms_morphe([atom]):
                    if charges == 'AB':      # we move the charges from state A to state B
                        atom.qqA = atom.q
                        atom.qqB = atom.qB
                        if not full_morphe and (atom.q*atom.qB < 0 or atom.atomtype!=atom.atomtypeB):   # we change a charge from + to - or vice versa
                            atom.qqB = 0
                            atom.to_be_morphed = True
                        else:
                            atom.qqB = atom.qB
                    elif charges == 'AA':        # we keep the charges
                        atom.qqA = atom.q
                        atom.qqB = atom.q

                    elif charges == 'BB':        # take charges of state B
                        if not full_morphe:
                            if hasattr(atom,"contQ"):
                                atom.qqA = atom.contQ
                                atom.qqB = atom.qqA
                            if hasattr(atom,"to_be_morphed"): # this a big q morphe. has been set to zero before
                                if vdw == 'BB':
                                    atom.qqA = 0
                                    atom.qqB = atom.qB
                                elif vdw == 'AB':
                                    atom.qqA = 0
                                    atom.qqB = 0
                            elif not hasattr(atom,"contQ") and not hasattr(atom,"to_be_morphed") :
                                atom.qqA = atom.qB
                                atom.qqB = atom.qB
                        else:
                            atom.qqA = atom.qB
                            atom.qqB = atom.qB
                    if atom.atomtype.startswith('DUM') or atom.atomtypeB.startswith('DUM'):
                        if dummy_chargeA == 'off':
                            atom.qqA = 0.
                        if dummy_chargeB == 'off':
                            atom.qqB = 0.

                else:
                    atom.qqA = atom.q
                    atom.qqB = atom.q
            qA_tot = sum(map(lambda a: a.qqA, r.atoms))
            qB_tot = sum(map(lambda a: a.qqB, r.atoms))
            if qB_tot != target_chargeB:
                do_log(logfile,'State B has total charge of %g' % round(qB_tot,5))
                do_log(logfile,'Applying charge correction to ensure integer charges')
                latom = last_perturbed_atom(logfile,r)
                do_log(logfile,'Selecting atom %d-%s (%s) as perturbed atom with highest order' % (latom.id,latom.name, latom.resname))
                newqB = latom.qqB-(qB_tot-target_chargeB)
                do_log(logfile,'Changing chargeB of atom %s from %g to %g' % (latom.name, latom.qqB,newqB))
                latom.qqB = newqB
                qB_tot = sum(map(lambda a: a.qqB, r.atoms))
                do_log(logfile,'New total charge of B-state is %g' % round(qB_tot,5))
            else:
                do_log(logfile,'No corrections applied to ensure integer charges')
#    sys.exit()
    print >>fp,'\n [ atoms ]'
    print >>fp, ';   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB'
    al = model.atoms
    for atom in al:
        if atoms_morphe([atom]):

            if vdw == 'AB':
                atA = atom.atomtype
                atB = atom.atomtypeB
                mA = atom.m
                mB = atom.mB
            elif vdw == 'AA':
                atA = atom.atomtype
                atB = atom.atomtype
                mA = atom.m
                mB = atom.m
            elif vdw == 'BB':
                atA = atom.atomtypeB
                atB = atom.atomtypeB
                mA = atom.mB
                mB = atom.mB
            if scale_mass:
                if atA.startswith('DUM'):
                    mA = 1.
                if atB.startswith('DUM'):
                    mB = 1.
            if hasattr(atom,"qqB"):
                qqB = atom.qqB
                if hasattr(atom,"contQ") and not full_morphe:
                    qqA = atom.contQ
                else:
                    qqA = atom.qqA
            else:
                qqA = atom.q
                qqB = atom.qB
            print >>fp , '%6d%11s%7d%7s%7s%7d%11.6f%11.4f%11s%11.6f%11.4f' % \
                  (atom.id, atA, atom.resnr, atom.resname, atom.name, \
                   atom.cgnr, qqA, mA, atB, qqB, mB)
            QA+=qqA
            QB+=qqB
        else:
            print >>fp , '%6d%11s%7d%7s%7s%7d%11.6f%11.4f' % \
                  (atom.id, atom.atomtype, atom.resnr, atom.resname, atom.name, \
                   atom.cgnr, atom.q, atom.m)
            QA+=atom.q
            QB+=atom.q
    # write qB of latom to qA
    if not full_morphe:
        try:
            latom.contQ = latom.qqB
        except:
            pass

def write_bonds(fp,bonds, state = 'AB', shrink_bonds = False):
    print >>fp,'\n [ bonds ]'
    print >>fp, ';  ai    aj funct            c0            c1            c2            c3'
    for b in bonds:
        if len(b) == 3:
            print >>fp, '%6d %6d %6d' % (b[0].id, b[1].id, b[2])
        else:
            lA = b[3][1]
            kA = b[3][2]
            lB = b[4][1]
            kB = b[4][2]
            if shrink_bonds:
                if b[0].atomtype.startswith('DUM') or \
                   b[1].atomtype.startswith('DUM'): lA*=.5
                if (b[0].atomtypeB is not None and b[0].atomtypeB.startswith('DUM')) or \
                   (b[1].atomtypeB is not None and b[1].atomtypeB.startswith('DUM')): lB*=.5
            if state == 'AB':
                print >>fp, '%6d %6d %6d %14.6f %14.6f %14.6f %14.6f' % \
                      (b[0].id, b[1].id, b[2],lA,kA, lB, kB)
            elif state == 'AA':
                print >>fp, '%6d %6d %6d %14.6f %14.6f %14.6f %14.6f' % \
                      (b[0].id, b[1].id, b[2],lA, kA, lA, kA)
            elif state == 'BB':
                print >>fp, '%6d %6d %6d %14.6f %14.6f %14.6f %14.6f' % \
                      (b[0].id, b[1].id, b[2],lB, kB, lB, kB)

def write_pairs(fp,pairs):
    print >>fp,'\n [ pairs ]'
    print >>fp, ';  ai    aj funct            c0            c1            c2            c3'
    for p in pairs:
        print >>fp, '%6d %6d %6d' % (p[0].id, p[1].id, p[2])
        


def write_angles(fp,angles, state='AB'):
    print >>fp,'\n [ angles ]'    
    print >>fp, ';  ai    aj    ak funct            c0            c1            c2            c3'
    for ang in angles:
        if len(ang) == 4:
            print >>fp, '%6d %6d %6d %6d' % (ang[0].id, ang[1].id, ang[2].id,ang[3])
        else:
            try:
                if state == 'AB':
                    print >>fp, '%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f ; %s %s %s' % \
                          (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[4][1], \
                           ang[4][2], ang[5][1], ang[5][2], ang[0].name, ang[1].name, ang[2].name)
                elif state == 'AA':
                    print >>fp, '%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f ; %s %s %s' % \
                          (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[4][1], \
                           ang[4][2], ang[4][1], ang[4][2], ang[0].name, ang[1].name, ang[2].name)
                elif state == 'BB':
                    print >>fp, '%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f ; %s %s %s' % \
                          (ang[0].id, ang[1].id, ang[2].id,ang[3], ang[5][1], \
                           ang[5][2], ang[5][1], ang[5][2], ang[0].name, ang[1].name, ang[2].name)
            except:
                print 'Error: Trouble with writing this angle!' 
                print ang[0].resname
                print ang[0].name, ang[1].name, ang[2].name
                print ang[0].atomtype, ang[1].atomtype, ang[2].atomtype
                print ang[0].atomtypeB, ang[1].atomtypeB, ang[2].atomtypeB
                print ang[3:]
                sys.exit(1)
                
def write_dihedrals(fp, dihedrals, state='AB'):
    print >>fp,'\n [ dihedrals ]'    
    print >>fp,';  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5'
    for d in dihedrals:
        if len(d) == 5:
            print >>fp, "%6d %6d %6d %6d %4d" % ( d[0].id, d[1].id, d[2].id, d[3].id, d[4])
        elif len(d) == 6:
            print >>fp, "%6d %6d %6d %6d %4d %s" % ( d[0].id, d[1].id, d[2].id, d[3].id, d[4], d[5])
        elif len(d) == 7:
            A, B = check_case(d[:4])
            ast = d[5]
            bs = d[6]
            if ast == None or bs == None:
                print d[0].name, d[1].name, d[2].name, d[3].name, d[0].atomtype, d[1].atomtype, d[2].atomtype, d[3].atomtype, d[0].atomtypeB, d[1].atomtypeB, d[2].atomtypeB, d[3].atomtypeB
                print d[0].type, d[1].type, d[2].type, d[3].type, d[0].typeB, d[1].typeB, d[2].typeB, d[3].typeB
            if ast == 'NULL':
                if d[4] == 3: # Ryckaert-Bellemans
#                    d[5] = [0,0,0,0,0,0]
                    ast = ' '.join(["%g" % x for x in [0,0,0,0,0,0]])
#                    ast = ' '.join(["%g" % x for x in d[5]]) 
                elif d[4] == 1:
#                    d[5] = [0,0,0]
                    ast = ' '.join(["%g" % x for x in [0,0,0]])
#                    ast = ' '.join(["%g" % x for x in d[5]]) 
            elif ast != 'NULL' and hasattr(ast,"append"):
                ast = ' '.join(["%g" % x for x in d[5][1:]])
            if bs == 'NULL':
                if d[4] == 3:
#                    d[6] = [0,0,0,0,0,0]
                    bs = ' '.join(["%g" % x for x in [0,0,0,0,0,0]]) 
#                    bs = ' '.join(["%g" % x for x in d[6]]) 
                elif d[4] == 1:
#                    d[6] = [0,0,0]
                    bs = ' '.join(["%g" % x for x in [0,0,0]])
#                    bs = ' '.join(["%g" % x for x in d[6]])
            elif bs !='NULL' and hasattr(bs,"append"):
                bs = ' '.join(["%g" % x for x in d[6][1:]])
            if state == 'AB':
                print >>fp, "%6d %6d %6d %6d %4d %s %s ; %s %s %s %s %s %s %s %s (%s->%s)" % \
                      ( d[0].id, d[1].id, d[2].id, d[3].id, d[4], ast, bs, d[0].name,d[1].name,d[2].name,d[3].name, \
                        d[0].type,d[1].type,d[2].type,d[3].type,A,B)
            elif state == 'AA':
                print >>fp, "%6d %6d %6d %6d %4d %s %s ; %s %s %s %s %s %s %s %s (%s->%s)" % \
                      ( d[0].id, d[1].id, d[2].id, d[3].id, d[4], ast, ast, d[0].name,d[1].name,d[2].name,d[3].name, \
                        d[0].type,d[1].type,d[2].type,d[3].type, A,B)
            elif state == 'BB':
                print >>fp, "%6d %6d %6d %6d %4d %s %s ; %s %s %s %s %s %s %s %s (%s->%s)" % \
                      ( d[0].id, d[1].id, d[2].id, d[3].id, d[4], bs, bs, d[0].name,d[1].name,d[2].name,d[3].name, \
                        d[0].type,d[1].type,d[2].type,d[3].type, A,B)

def find_predefined_dihedrals(dihedrals, rlist, rdic, dih_lib):
    for r in rlist:
        idx = r.id - 1
        dih = rdic[r.resname][3]
        imp = rdic[r.resname][2]
        for d in dih:
            al = []
            for name in d[:4]:
                if name.startswith('+'):
                    next = r.chain.residues[idx+1]
                    if name == '+H' and next.resname=='PRO':
                        name = '+CD'
##                     try:
##                         atom = next.fetch(name[1:])[0]
##                     except:
##                         print next, name, d
##                         sys.exit()
                    al.append(atom)
                elif name.startswith('-'):
                    prev = r.chain.residues[idx-1]
                    atom = prev.fetch(name[1:])[0]
                    al.append(atom)
                else:
                    atom = r.fetch(name)[0]
                    al.append(atom)
            for dx in dihedrals:
                func = dx[4]
                if (dx[0].id == al[0].id and \
                   dx[1].id == al[1].id and \
                   dx[2].id == al[2].id and \
                   dx[3].id == al[3].id) or \
                   (dx[0].id == al[3].id and \
                    dx[1].id == al[2].id and \
                    dx[2].id == al[1].id and \
                    dx[3].id == al[0].id):
                    A,B =  check_case(al[:4])
                    if d[4] == 'default':
                        if A=='AAAA':
                            astate = get_dihedral_param(al[0].type,al[1].type,al[2].type,al[3].type,dih_lib, func)
                            if astate is None:
                                print 'Error: Dihedral parameters (state A) not found for:'
                                print al[0].resname,al[1].resname,al[2].resname,al[3].resname
                                print al[0].atomtype,al[1].atomtype,al[2].atomtype,al[3].atomtype, func
                                print al[0].type,al[1].type,al[2].type,al[3].type, func
                                print al[0].atomtypeB,al[1].atomtypeB,al[2].atomtypeB,al[3].atomtypeB, func
                                print al[0].typeB,al[1].typeB,al[2].typeB,al[3].typeB, func
                                print dx
                                sys.exit(1)
                        else:   # contains dummies
                            astate = d[5]
                    else:
                        astate = d[4]
                    if d[5] == 'default':
                        if B=='AAAA':
#                            print al[0].typeB,al[1].typeB,al[2].typeB,al[3].typeB
                            bstate = get_dihedral_param(al[0].typeB,al[1].typeB,al[2].typeB,al[3].typeB,dih_lib, func)
                            if bstate is None:
                                print 'Error: Dihedral parameters (state B) not found for:'
                                print al[0].resname,al[1].resname,al[2].resname,al[3].resname
                                print al[0].atomtype,al[1].atomtype,al[2].atomtype,al[3].atomtype, func
                                print al[0].type,al[1].type,al[2].type,al[3].type, func
                                print al[0].atomtypeB,al[1].atomtypeB,al[2].atomtypeB,al[3].atomtypeB, func
                                print al[0].typeB,al[1].typeB,al[2].typeB,al[3].typeB, func
                                print dx
                                sys.exit(1)
                        else:
                            bstate = d[4]
                    else:
                        bstate = d[5]
                    if len(dx) > 5:
                        dx[5] = astate
                    else: dx.append(astate)
                    if len(dx) > 6:
                        dx[6] = bstate
                    else: dx.append(bstate)

        for d in imp:
            al = []
            for name in d[:4]:
                if name.startswith('+'):
                    next = r.chain.residues[idx+1]
                    atom = next.fetch(name[1:])[0]
                    al.append(atom)
                elif name.startswith('-'):
                    prev = r.chain.residues[idx-1]
                    atom = prev.fetch(name[1:])[0]
                    al.append(atom)
                else:
                    atom = r.fetch(name)[0]
                    al.append(atom)
            for dx in dihedrals:
                func = dx[4]
                if (dx[0].id == al[0].id and \
                   dx[1].id == al[1].id and \
                   dx[2].id == al[2].id and \
                   dx[3].id == al[3].id) or \
                   (dx[0].id == al[3].id and \
                    dx[1].id == al[2].id and \
                    dx[2].id == al[1].id and \
                    dx[3].id == al[0].id):
                    A,B =  check_case(al[:4])

                    if d[4] == 'default':
                        if A=='AAAA':       # it's a real improper in state A
                            astate = get_dihedral_param(al[0].type,al[1].type,al[2].type,al[3].type,dih_lib, func)
                            if astate is None:
                                print 'Error: Improper parameters (state A) not found for:'
                                print al[0].resname,al[1].resname,al[2].resname,al[3].resname
                                print al[0].atomtype,al[1].atomtype,al[2].atomtype,al[3].atomtype, func
                                print al[0].type,al[1].type,al[2].type,al[3].type, func
                                print al[0].atomtypeB,al[1].atomtypeB,al[2].atomtypeB,al[3].atomtypeB, func
                                print al[0].typeB,al[1].typeB,al[2].typeB,al[3].typeB, func
                                print dx
                                sys.exit(1)
                        else:   # contains dummies -> we take the bstate
                            astate = 'bstate'
                    elif d[4] == 'undefined': # no improper potential in state A
                        if A == 'AAAA':     
                            astate = [1, 0,0,2]    # zero potential
                        else:  # contains dummies -> we simply take the bstate
##                             print 'taking bstate', al[0].name, al[1].name, al[2].name, al[3].name, d[4],d[5], A, B, func
##                             print al[0].atomtype, al[1].atomtype, al[2].atomtype, al[3].atomtype
##                             print al[0].atomtypeB, al[1].atomtypeB, al[2].atomtypeB, al[3].atomtypeB
                            astate = 'bstate'
                    else:                   # defined improper
                        astate = d[4]
                        
                    if d[5] == 'default':
                        if B=='AAAA':      # it's a real improper in state B
                            bstate = get_dihedral_param(al[0].typeB,al[1].typeB,al[2].typeB,al[3].typeB,dih_lib, func)
                            if bstate is None:
                                print 'Error: Dihedral parameters (state B) not found for:'
                                print al[0].resname,al[1].resname,al[2].resname,al[3].resname
                                print al[0].atomtype,al[1].atomtype,al[2].atomtype,al[3].atomtype, func
                                print al[0].type,al[1].type,al[2].type,al[3].type, func
                                print al[0].atomtypeB,al[1].atomtypeB,al[2].atomtypeB,al[3].atomtypeB, func
                                print al[0].typeB,al[1].typeB,al[2].typeB,al[3].typeB, func
                                print dx
                                sys.exit(1)
                            if astate == 'bstate':
                                astate = bstate
                        else:  # contains dummies -> we take the astate
                            bstate = astate
                    elif d[5] == 'undefined':
                        if B=='AAAA':    # no improper potential in state B
                            bstate = [1, 0,0,2]    # zero potential
                        else:
                            bstate = astate   # take A state
                    else:                   # take predefined value
                        bstate = d[5]
                        if astate == 'bstate': astate = bstate

                    # write values to dx
                    if len(dx) > 5:
                        dx[5] = astate
                    else: dx.append(astate)
                    if len(dx) > 6:
                        dx[6] = bstate
                    else: dx.append(bstate)
    return dihedrals


def get_bond_param(type1,type2,bond_lib):
    for entr in bond_lib:
        if (type1==entr[0] and type2==entr[1]) or \
           (type2==entr[0] and type1==entr[1]):
            return entr[2:]
    return None


def get_angle_param(type1,type2,type3,ang_lib):
    for entr in ang_lib:
        if (type1 == entr[0] and \
            type2 == entr[1] and \
            type3 == entr[2]) or \
           (type1 == entr[2] and \
            type2 == entr[1] and \
            type3 == entr[0]):
            return entr[3:]
    return None 

def get_dihedral_param(type1,type2,type3,type4,dih_lib, func):
    for entr in dih_lib:
        if (type1 == entr[0] and \
            type2 == entr[1] and \
            type3 == entr[2] and \
            type4 == entr[3] and func==entr[4]) or \
           (type1 == entr[3] and \
            type2 == entr[2] and \
            type3 == entr[1] and \
            type4 == entr[0] and func==entr[4]):
            return entr[4:]
    for entr in dih_lib:
        if ('X' == entr[0] and \
            type2 == entr[1] and \
            type3 == entr[2] and \
            type4 == entr[3] and func==entr[4]) or \
           (type1 == entr[0] and \
            type2 == entr[1] and \
            type3 == entr[2] and \
            'X' == entr[3] and func==entr[4]):
            return entr[4:]
        if ('X' == entr[3] and \
            type2 == entr[2] and \
            type3 == entr[1] and \
            type4 == entr[0] and func==entr[4]) or \
           (type1 == entr[3] and \
            type2 == entr[2] and \
            type3 == entr[1] and \
            'X' == entr[0] and func==entr[4]):
            return entr[4:]
    for entr in dih_lib:
        if ('X' == entr[0] and \
            type2 == entr[1] and \
            type3 == entr[2] and \
            'X' == entr[3] and func==entr[4]) or \
           ('X' == entr[3] and \
            type2 == entr[2] and \
            type3 == entr[1] and \
            'X' == entr[0] and func==entr[4]):
            return entr[4:]
    for entr in dih_lib:
        if ('X' == entr[0] and \
            'X' == entr[1] and \
            type3 == entr[2] and \
            type4 == entr[3] and func==entr[4]) or \
           (type1 == entr[3] and \
            type2 == entr[2] and \
            'X' == entr[1] and \
            'X' == entr[0] and func==entr[4]):
            return entr[4:]
    # check whether we have a predefined dihedral
    
    return None 
        

def topline2atom(line):
    entr = line.split()
    idx = int(entr[0])
    atomtype = entr[1]
    resnr = int(entr[2])
    resname = entr[3]
    name = entr[4]
    cgnr = int(entr[5])
    q = float(entr[6])
    m = float(entr[7])
    try:
        atomtypeB = entr[8]
        qB = float(entr[9])
        mB = float(entr[10])
    except:
        atomtypeB = None
        qB = None
        mB = None
    a = Atom(id=idx,atomtype=atomtype,\
        resnr = resnr, resname = resname,\
        name = name, cgnr = cgnr, q = q, \
        m = m, atomtypeB = atomtypeB, \
        qB = qB, mB = mB)
    return a

def read_top_header(lines):
    ret = []
    for line in lines:
        if not line.strip().startswith('[ atoms ]'):
            ret.append(line.rstrip())
        else:
            break
    return ret

def write_header(fp,lines):
    for line in lines:
        print >>fp, line

def write_footer(fp,lines):
    for line in lines:
        print >>fp, line
        
def read_top_footer(lines):
    for line in lines:
        if line.strip().startswith('#ifdef POSRES'):
            idx = lines.index(line)
            return [l.rstrip() for l in lines[idx:]]

        
def read_atoms(lines):
    lst = readSection(lines,'[ atoms ]','[')
    al = []
    for line in lst:
        a = topline2atom(line)
        al.append(a)
    return al

def read_bonds(lines, al):
    lst = readSection(lines,'[ bonds ]','[')
    bonds = []
    for line in lst:
        idx = [int(x) for x in line.split()]
        bonds.append([al[idx[0]-1], al[idx[1]-1], idx[2]])
    return bonds

def read_pairs(lines, al):
    lst = readSection(lines,'[ pairs ]','[')
    pairs = []
    for line in lst:
        idx = [int(x) for x in line.split()]
        pairs.append([al[idx[0]-1], al[idx[1]-1], idx[2]])
    return pairs

def read_angles(lines, al):
    lst = readSection(lines,'[ angles ]','[')
    angles = []
    for line in lst:
        idx = [int(x) for x in line.split()]
        angles.append([al[idx[0]-1], al[idx[1]-1], al[idx[2]-1], idx[3]])
    return angles

def read_dihedrals(lines, al):
    starts = []
    dih = []
    for i, line in enumerate(lines):
        if line.strip().startswith('[ dihedrals ]'):
            starts.append(i)
    for s in starts:
        lst = readSection(lines[s:],'[ dihedrals ]','[')
        for line in lst:
            entr = line.split()
            idx = [int(x) for x in entr[:4]]
            
            func = int(entr[4])
            try:
                rest = entr[5]
            except:
                rest = ''
            dih.append([al[idx[0]-1],al[idx[1]-1],al[idx[2]-1],\
                        al[idx[3]-1],func,rest])
    return dih

def assign_fftypes(al, types):
    for atom in al:
        atom.type = types[atom.atomtype][0]
        if atom.atomtypeB is not None:
            atom.typeB = types[atom.atomtypeB][0]
        else:
            atom.typeB = types[atom.atomtype][0]

def find_bonded_entries(bonds, bond_lib):
    
    for b in bonds:
        a1,a2,func = b
        A, B = check_case([a1,a2])
        if a1.atomtypeB is not None or a2.atomtypeB is not None:
            astate = None
            bstate = None
            if A == 'AA' and B == 'AA': # we need A and B state
                astate = get_bond_param(a1.type,a2.type,bond_lib)
                bstate = get_bond_param(a1.typeB,a2.typeB,bond_lib)
                b.extend([astate,bstate])
            elif 'D' in A and B=='AA':
                bstate = get_bond_param(a1.typeB,a2.typeB,bond_lib)
                astate = bstate
                b.extend([astate,bstate])
            elif 'D' in B and A=='AA':
                astate = get_bond_param(a1.type,a2.type,bond_lib)
                bstate = astate
                b.extend([astate,bstate])
            # catch errors
            elif 'D' in B and 'D' in A:
                print 'Error: fake bond %s-%s (%s-%s -> %s-%s)' % (a1.name,a2.name,a1.atomtype,a2.atomtype, a1.atomtypeB,a2.atomtypeB)
                sys.exit(1)
            if astate == None:
                print 'Error: No bond entry found for astate %s-%s (%s-%s -> %s-%s)' % (a1.name,a2.name,a1.atomtype,a2.atomtype, a1.type,a2.type)
                print A, B
                sys.exit(1)
            if bstate == None:
                print 'Error: No bond entry found for bstate %s-%s (%s-%s -> %s-%s)' % (a1.name,a2.name,a1.atomtypeB,a2.atomtypeB, a1.typeB,a2.typeB)
                print A, B
                sys.exit(1)
    return bonds

def find_angle_entries(ang, ang_lib):

    ret = []
    for a in ang:
        a1,a2,a3,func = a
        astate = None
        bstate = None
        if a1.atomtypeB is not None or \
           a2.atomtypeB is not None or \
           a3.atomtypeB is not None:
            A, B = check_case([a1,a2,a3])
            if 'D' in A and 'D' in B: # fake angle
                astate = [1,0,0]
                bstate = [1,0,0]
                a.extend([astate,bstate])
            elif A == 'AAA' and 'D' in B: # take astate
                astate = get_angle_param(a1.type,a2.type,a3.type,ang_lib)
                bstate = astate
                a.extend([astate,bstate])
            elif 'D' in A and B == 'AAA': # take bstate
                bstate = get_angle_param(a1.typeB,a2.typeB,a3.typeB,ang_lib)
                astate = bstate
                a.extend([astate,bstate])
            elif A == 'AAA' and B == 'AAA':
                if a1.atomtypeB != a1.atomtype or \
                   a2.atomtypeB != a2.atomtype or \
                   a3.atomtypeB != a3.atomtype:
                    astate = get_angle_param(a1.type,a2.type,a3.type,ang_lib)
                    bstate = get_angle_param(a1.typeB,a2.typeB,a3.typeB,ang_lib)
                    a.extend([astate, bstate])
                else:
                    astate = get_angle_param(a1.type,a2.type,a3.type,ang_lib)
                    bstate = astate
                    a.extend([astate, bstate])
            if astate is None:
                print 'Error: No angle entry (state A) found for:'
                print A, B
                print a1.resname, a2.resname, a3.resname
                print a1.name, a2.name, a3.name
                print a1.atomtype, a2.atomtype, a3.atomtype
                print a1.type, a2.type, a3.type
                print a1.atomtypeB, a2.atomtypeB, a3.atomtypeB
                print a1.typeB, a2.typeB, a3.typeB
                sys.exit(1)
            if bstate is None:
                print 'Error: No angle entry (state B) found for:'
                print a1.resname, a2.resname, a3.resname
                print a1.name, a2.name, a3.name
                print a1.atomtype, a2.atomtype, a3.atomtype
                print a1.type, a2.type, a3.type
                print a1.atomtypeB, a2.atomtypeB, a3.atomtypeB
                print a1.typeB, a2.typeB, a3.typeB
                sys.exit(1)
                                
    return ang


def check_case(atoms):
    
    A = ''
    B = ''
    for a in atoms:
        if a.atomtype.startswith('DUM'): A += 'D'
        else: A += 'A'
        if a.atomtypeB is not None:
            if a.atomtypeB.startswith('DUM'): B += 'D'
            else: B += 'A'
        else: B += 'A'
    
    return A, B

def atoms_morphe(atoms):
    for atom in atoms:
        if atom.atomtypeB is not None and (atom.q!=atom.qB or atom.m != atom.mB): return True
    return False

def types_morphe(atoms):
    for atom in atoms:
        if atom.atomtypeB is not None and atom.atomtype != atom.atomtypeB: return True
    return False

def find_dihedral_entries(dih,dih_lib):

    for d in dih:
        if len(d) == 6:
            a1,a2,a3,a4, func, val = d
            if atoms_morphe([a1,a2,a3,a4]):
                A,B= check_case(d[:4])
                if A!='AAAA' and B!='AAAA':
                    # these are fake dihedrals
                    d[5] = 'NULL'
                    d.append('NULL')
                else:
                    if A == 'AAAA' and B!='AAAA':
                        astate = get_dihedral_param(a1.type,a2.type,a3.type,a4.type,dih_lib, func)
                        bstate = astate
                    elif B == 'AAAA' and A!='AAAA':
                        bstate = get_dihedral_param(a1.typeB,a2.typeB,a3.typeB,a4.typeB,dih_lib, func)
                        astate = bstate
                    elif A=='AAAA' and B=='AAAA':
                        if val=='':
                            astate = get_dihedral_param(a1.type,a2.type,a3.type,a4.type,dih_lib, func)
                        else:
                            astate = val
                        if types_morphe([a1,a2,a3,a4]):
                            bstate = get_dihedral_param(a1.typeB,a2.typeB,a3.typeB,a4.typeB,dih_lib, func)
                        else:
                            bstate = astate
 
                    if astate == None :
                        print 'Error: No dihedral angle found (state A: predefined state B) for:' 
                        print a1.resname, a2.resname, a3.resname, a4.resname
                        print a1.name, a2.name, a3.name, a4.name, func
                        print a1.atomtype, a2.atomtype, a3.atomtype, a4.atomtype
                        print a1.type, a2.type, a3.type, a4.type
                        print a1.atomtypeB, a2.atomtypeB, a3.atomtypeB, a4.atomtypeB
                        print a1.typeB, a2.typeB, a3.typeB, a4.typeB
                        print d
                        sys.exit(1)
                        
                    if bstate == None :
                        print 'Error: No dihedral angle found (state B: predefined state A) for:' 
                        print a1.resname, a2.resname, a3.resname, a4.resname
                        print a1.name, a2.name, a3.name, a4.name, func
                        print a1.atomtype, a2.atomtype, a3.atomtype, a4.atomtype
                        print a1.type, a2.type, a3.type, a4.type
                        print a1.atomtypeB, a2.atomtypeB, a3.atomtypeB, a4.atomtypeB
                        print a1.typeB, a2.typeB, a3.typeB, a4.typeB
                        print d
                        sys.exit(1)
                        
                    d[5] = astate
                    d.append(bstate)

    return dih
                   
def find_mut_residues(m, ff = 'amber03'):
    fname = 'ff'+ff+'.mtp'
    rdic = {}
    rlist = []
    for r in m.residues:
        if r.resname[1]=='2' or r.resname in NUC_MUTS:
            rlist.append(r)
            entr = mutdb.read_mtp_entry(r.resname,fname)
            rdic[r.resname] = entr
            mol = entr[0]
            atomnames = map(lambda a: a.name, mol.atoms)
            atoms = r.fetchm(atomnames)
            for i, atom in enumerate(mol.atoms):
                atoms[i].atomtypeB = atom.atomtypeB
                atoms[i].qB = atom.qB
                atoms[i].mB = atom.mB
    return rlist,rdic


def get_extra_dihedrals(logfile, rlist, ff = 'amber99sb'):
    extra_dih = []
    for r in rlist:
        if r.resname in ['DAT','DAC','DGC','DGT']:
            do_log(logfile, "Adding extra dihedrals for residue %d-%s" % (r.id, r.resname))
            alist = r.fetchm(['C1\'','N9','C8','DC2'])
            extra_dih.append( alist )
            alist = r.fetchm(['C1\'','N9','C4','DC6'])
            extra_dih.append( alist )
        elif r.resname in ['DTA','DTG','DCG','DCA']:
            do_log(logfile, "Adding extra dihedrals for residue %d-%s" % (r.id, r.resname))
            alist = r.fetchm(['C1\'','N1','C6','DC4'])
            extra_dih.append( alist )
            alist = r.fetchm(['C1\'','N1','C2','DC8'])
            extra_dih.append( alist )
    return  extra_dih

def write_extra_dihedrals(fp, extra_dih, dihtype, stateA, stateB):
    for dih in extra_dih:
        if hasattr(stateA,"append"):
            astate = ' '.join([str(x) for x in stateA])
        else:
            astate = stateA
        if hasattr(stateB,"append"):
            bstate = ' '.join([str(x) for x in stateB])
        else:
            bstate = stateB
        
        print >>fp, "%6d %6d %6d %6d %4d %s %s ; extra dihedral" % \
              (dih[0].id, dih[1].id, dih[2].id, dih[3].id, dihtype, astate, bstate)
    
def get_extra_angles(logfile, rlist, ff = 'amber99sb'):
    extra_ang = []
    for r in rlist:
        if r.resname in ['DAT','DAC','DGC','DGT']:
            do_log(logfile, "Adding extra angles for residue %d-%s" % (r.id, r.resname))
            alist = r.fetchm(['C8','N9','DC6'])
            extra_ang.append( alist )
            alist = r.fetchm(['C8','N9','DC6'])
            extra_ang.append( alist )
        elif r.resname in ['DTA','DTG','DCG','DCA']:
            do_log(logfile, "Adding extra angles for residue %d-%s" % (r.id, r.resname))
            alist = r.fetchm(['C6','N1','DC8'])
            extra_ang.append( alist )
            alist = r.fetchm(['C6','N1','DC8'])
            extra_ang.append( alist )
    return  extra_ang
        
def write_extra_angles(fp, extra_ang, angtype, stateA, stateB):
    for ang in extra_ang:
        if hasattr(stateA,"append"):
            astate = ' '.join([str(x) for x in stateA])
        else:
            astate = stateA
        if hasattr(stateB,"append"):
            bstate = ' '.join([str(x) for x in stateB])
        else:
            bstate = stateB
        
        print >>fp, "%6d %6d %6d %4d %s %s ; extra angle" % \
              (ang[0].id, ang[1].id, ang[2].id, angtype, astate, bstate)
    
    


def main(argv):

    desc = ('')

    files= [
        (efTOP, "-p", "topol.top", ffREAD ),
        (efITP, "-itp", "topol_A.itp", ffOPTRD ),
        (efTOP, "-o", "newtop.top", ffOPTWR ),
        (efITP, "-oitp", "newtop.itp", ffOPTWR ),
        (efLOG, "-log", "bstate.log", ffWRITE ),
        (efDAT, "-qtrace", "qtrace.dat", ffWRITE ),
        ]
    
    options = [
        ('-ff',FALSE,etSTR,'force field (oplsaa, amber03 or amber99sb)','amber99sb'),
        ('-split',FALSE,etBOOL,'Write splitted topologies for vdw and q morphes',False),
        ('-eval',FALSE,etBOOL,'Evaluate #define directives',False),
        ('-shrink',FALSE,etBOOL,'Shrink bonds of dummy atoms',False),
        ('-scale_mass',FALSE,etBOOL,'Scale down dummy masses',True),
        ]
    
    args=parse_common_args(argv,files,options, desc)
    scale_m = args['-scale_mass']['value']
    bShrink = args['-shrink']['value']
    #ff = args['-ff']['value']
    ff= "amber99sb"
    topfile = "topol.top"
    #topfile = args['-p']['fns']
    logfile = open(args['-log']['fns'],'w')
    qfile = open(args['-qtrace']['fns'],'w')
    
    bEval = args['-eval']['value']
    do_log(logfile, "Log file opened at %s" % time.asctime())
    do_log(logfile, "Command: %s" % ' '.join(sys.argv))


    
    if bEval:
        proc = cpp.PreProcessor()
        proc(topfile)
        cpp_dic = proc.cpp_namespace

    if ff not in ['oplsaa','amber03','amber99sb']:
        raise ForceFieldError('Sorry only oplsaa, amber03 and amber99sb are supported!')
    

    if '-itp' in argv:
        outfile = args['-oitp']['fns']
        infile = args['-itp']['fns']

    else:
        infile = topfile
        outfile = args['-o']['fns']

    do_log(logfile,'Reading Topology "%s"' % infile)
    lines = open(infile).readlines()

    header = read_top_header(lines)
    footer = read_top_footer(lines)

    l = kickOutComments(lines,';#')
    al = read_atoms(l)
    m = Model(atoms = al)
    rlist, rdic = find_mut_residues(m, ff = ff)
    do_log(logfile,"Adding B-state parameters for these residues:")

    for r in rlist:
        do_log(logfile," --> %d-%s" % (r.id, r.resname))

    bonds = read_bonds(l, al)
    pairs = read_pairs(l, al)
    ang = read_angles(l, al)
    dih = read_dihedrals(l,al)
    do_log(logfile,'Reading force field parameters from "%s"' % topfile)
    types, bond_lib, ang_lib, dih_lib = forcefield.read_ff(topfile,ff=ff)
    
    assign_fftypes(al,types)
    bonds = find_bonded_entries(bonds,bond_lib)
    angles = find_angle_entries(ang, ang_lib)
    dih = find_predefined_dihedrals(dih,rlist,rdic,dih_lib)
    dih = find_dihedral_entries(dih,dih_lib)
    if bEval:
        for d in dih:
            if len(d) == 7:
                if not hasattr(d[5],"append"):
                    if cpp_dic.has_key(d[5]):
                        d[5] = cpp_dic[d[5]]
                if not hasattr(d[6],"append"):
                    if cpp_dic.has_key(d[6]):
                        d[6] = cpp_dic[d[6]]

    qtot_A, qtot_B = sum_charge_of_states(rlist)
    qAsave = copy.deepcopy(qtot_A)
    qBsave = copy.deepcopy(qtot_B)
    do_log(logfile,"Making topology for full switch")
    do_log(logfile,'Writing "%s"' % outfile)
    fp = open(outfile,'w')
    write_header(fp,header)
    write_atoms(fp,logfile,m,scale_mass=scale_m,full_morphe=True,tqB = qBsave)
    do_log(logfile,"Total charge of topology:")
    do_log(logfile,"State A -> %4.2f"%QA)
    do_log(logfile,"State B -> %4.2f"%QB)
    print >>qfile," complete = %d" %( int(round((QB-QA),0)))
    write_bonds(fp,bonds, shrink_bonds = bShrink)
    write_pairs(fp,pairs)
    write_angles(fp,angles)

#    extra_ang = get_extra_angles(logfile,rlist,ff)
#    if extra_ang:
#        write_extra_angles(fp, extra_ang, 1,[0,800],[0,800])
    write_dihedrals(fp,dih)
    extra_dih = get_extra_dihedrals(logfile,rlist,ff)
    if extra_dih:
        write_extra_dihedrals(fp, extra_dih, 1, [180,40,2],[180,40,2])
    write_footer(fp,footer)
    fp.close()
#    qAsave = copy.deepcopy(qtot_A)
    qAsave = []
    for i in range(len(qtot_A)):
        print qtot_A[i], qtot_B[i]
        if round(qtot_B[i],2) == round(qtot_A[i],2):
            qAsave.append(qtot_A[i])
        else:
            qAsave.append(0)
    qBsave = copy.deepcopy(qtot_B)
#    qtot_A, qtot_B = sum_charge_of_states(rlist)
    do_log(logfile,"Full morphe done")
    do_log(logfile,"----------------------------------------")

    
    if args['-split']['value']:
        do_log(logfile,"Splitting topology into three components")
        do_log(logfile,"Topology 1 (qqoff): Charges of dummy atoms are removed, LJ parameters are unchanged")
        qqoff = outfile.split('.')[0]+'_qqoff.'+outfile.split('.')[1]
        fp = open(qqoff,'w')
        do_log(logfile,'Writing "%s"' % qqoff)
        write_header(fp,header)
        contQ = copy.deepcopy(qAsave)
        write_atoms(fp,logfile,m,charges = 'AB', dummy_chargeB='off',vdw='AA',scale_mass=scale_m,tqB = qAsave)
        do_log(logfile,"Total charge of topology:")
        do_log(logfile,"State A -> %4.2f"%QA)
        do_log(logfile,"State B -> %4.2f"%QB)
        print >>qfile," qqoff = %d" %( int(round((QB-QA),0)))
        do_log(logfile,"----------------------------------------")
        write_bonds(fp,bonds, state = 'AA', shrink_bonds = bShrink)
        write_pairs(fp,pairs)
        write_angles(fp,angles, state='AA')
        write_dihedrals(fp,dih, state='AA')
        if extra_dih:
            write_extra_dihedrals(fp, extra_dih, 1, [180,40,2],[180,40,2])
            
#        if extra_ang:
#            write_extra_angles(fp, extra_ang, 1,[0,800],[0,800])
        write_footer(fp,footer)
        fp.close()
        vdw = outfile.split('.')[0]+'_vdw.'+outfile.split('.')[1]
        fp = open(vdw,'w')
        do_log(logfile,"Topology 2 (vdw): Charges are unchanged, LJ parameters are switched to B-states")
        do_log(logfile,'Writing "%s"' % vdw)
        do_log(logfile,"----------------------------------------")
#        qAsave = copy.deepcopy(qBsave)
#        print qAsave
#        qAsave = []
#        for i in range(len(qtot_A)):
#            qAsave.append(0)
        qBsave = copy.deepcopy(qtot_B)
        write_header(fp,header)
        write_atoms(fp,logfile,m,charges = 'BB' ,dummy_chargeA = 'off', dummy_chargeB = 'off',vdw='AB',scale_mass=scale_m,tqB = contQ)
        write_bonds(fp,bonds, state = 'AB', shrink_bonds = bShrink)
        write_pairs(fp,pairs)
        write_angles(fp,angles, state = 'AB')

#        if extra_ang:
#            write_extra_angles(fp, extra_ang, 1,[0,800],[0,800])

        write_dihedrals(fp,dih, state = 'AB')
        if extra_dih:
            write_extra_dihedrals(fp, extra_dih, 1, [180,40,2],[180,40,2])
        write_footer(fp,footer)
        do_log(logfile,"Total charge of topology:")
        do_log(logfile,"State A -> %4.2f"%QA)
        do_log(logfile,"State B -> %4.2f"%QB)
        print >>qfile," vdw = %d" % ( int(round((QB-QA),0)))
        do_log(logfile,"----------------------------------------")
        fp.close()
        qqon = outfile.split('.')[0]+'_qqon.'+outfile.split('.')[1]
        fp = open(qqon,'w')
        do_log(logfile,"Topology 3 (qqon): Charges are morphed to their B-state values, LJ parameters are unchanged")
        do_log(logfile,'Writing "%s"' % qqon)
        write_header(fp,header)
        write_atoms(fp,logfile,m,charges = 'BB',dummy_chargeA = 'off', \
                    vdw = 'BB' ,scale_mass=scale_m, tqB = qtot_B)
        do_log(logfile,"Total charge of topology:")
        do_log(logfile,"State A -> %4.2f" %QA)
        do_log(logfile,"State B -> %4.2f"%QB)
        print >>qfile," qqon = %d" % ( int(round((QB-QA),0)))
        do_log(logfile,"----------------------------------------")
        write_bonds(fp,bonds, state = 'BB', shrink_bonds = bShrink)
        write_pairs(fp,pairs)
        write_angles(fp,angles, state = 'BB')
        
#        if extra_ang:
#            write_extra_angles(fp, extra_ang, 1,[0,800],[0,800])
        write_dihedrals(fp,dih, state = 'BB')
        if extra_dih:
            write_extra_dihedrals(fp, extra_dih, 1, [180,40,2],[180,40,2])

        write_footer(fp,footer)
        fp.close()
        do_log(logfile,"Splitted topologies done")

    do_log(logfile,"----------------------------------------")
    do_log(logfile,"Successful termination")
    do_log(logfile,"----------------------------------------")
    thanx()
        
if __name__=='__main__':
    main(sys.argv)
