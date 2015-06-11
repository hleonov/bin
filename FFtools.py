import sys,os,copy
#from scipy import *
#from random import gauss as gaul
#from gaussfit import *

def ReadAtp(atp):
    '''
    Reads a GROMACS forcefield ATP file'''

    ifile=open(atp,'r').readlines()
    fftypes={}
    for line in ifile:
        if line[0]!=';':
            col=line.split()
            fftypes.update({col[0]:col[1]})
    return fftypes

def ReadNonbonded(itp,heavy=0):
    '''
    Reads a GROMACS forcefield nonbonded ITP file. The heavy flag set to 1 will yield the heavy water part of the forcefield instead of the standard.'''

    ifile=open(itp,'r').readlines()
    lstart=0
    lend=0
    for i in range(len(ifile)):
        if '[ atomtypes ]' in ifile[i]:
            lstart=i+1
        if ifile[i][0]=='[' and '[ atomtypes ]' not in ifile[i]:
            lend=i-1
    if not lstart:
        print '\n ERROR: Using wrong file?'
        sys.exit(1)
    if not lend:
        lend=len(ifile)

    fftypes={}
    hflag=False
    for i in range(lstart,lend):
        if '#ifdef HEAVY_H' in ifile[i] and not heavy:
            hflag=True
        if '#else' in ifile[i] and not heavy:
            hflag=False
        if '#else' in ifile[i] and heavy:
            hflag=True
        if '#endif' in ifile[i]:
            hflag=False
        if ifile[i][0]!=';' and ifile[i][0]!='#' and not hflag:
            col=ifile[i].split()
            ### If statement due to missing QM integers in the amber ff port ###
            ### Can be removed is that changes somewhen! ###
            if '.' in col[2]:
                fftypes.update({col[0]:{'type':col[1],
                                        'qmtype':'',
                                        'mass':col[2],
                                        'charge':col[3],
                                        'ptype':col[4],
                                        'sigma':col[5],
                                        'epsilon':col[6]}})
            else:
                fftypes.update({col[0]:{'type':col[1],
                                        'qmtype':col[2],
                                        'mass':col[3],
                                        'charge':col[4],
                                        'ptype':col[5],
                                        'sigma':col[6],
                                        'epsilon':col[7]}})
    return fftypes

def ReadBonded(itp):
    '''
    Reads a GROMACS forcefield bonded ITP file.'''

    ifile=open(itp,'r').readlines()
    sections=[]
    seclines=[]
    bonds={}
    angles={}
    propers={}
    impropers={}
    constraints={}
    ### Finding sections in file ###
    for i in range(len(ifile)):
        if ifile[i][0]=='[':
            sections.append(ifile[i].strip())
            seclines.append(i)
    seclines.append(len(ifile)+1)
    bcnt=0
    acnt=0
    pcnt=0
    icnt=0
    dihflag=0
    ### Parsing the sections for entries ###
    for i in range(len(sections)):
        if sections[i]=='[ bondtypes ]':
            for j in range(seclines[i]+1,seclines[i+1]-1):
                if ifile[j].strip() and ifile[j][0]!=';':
                    c=ifile[j].split()
                    bonds.update({bcnt:{'ai':c[0],
                                        'aj':c[1],
                                        'func':c[2],
                                        'bond':c[3],
                                        'fc':c[4]}})
                    bcnt+=1
            #print ifile[seclines[i]+1:seclines[i+1]-1]
        if sections[i]=='[ constrainttypes ]':
            ### NOT IMPLEMENTED YET ###
            pass
        if sections[i]=='[ angletypes ]':
            for j in range(seclines[i]+1,seclines[i+1]-1):
                if ifile[j].strip() and ifile[j][0]!=';':
                    c=ifile[j].split()
                    angles.update({acnt:{'ai':c[0],
                                         'aj':c[1],
                                         'ak':c[2],
                                         'func':c[3],
                                         'angle':c[4],
                                         'fc':c[5]}})
                    acnt+=1
        if sections[i]=='[ dihedraltypes ]':
            ### DEFINE ENTRIES AND FUNC TYPE 2 ARE NOT HANDLED YET ###
            dihflag+=1
            for j in range(seclines[i]+1,seclines[i+1]-1):
                if ifile[j].strip() and ifile[j][0]!=';' and ifile[j][0]!='#':
                    c=ifile[j].split()
                    if len(c)>3 and c[4]=='3':
                        if dihflag%2:
                            propers.update({pcnt:{'ai':c[0],
                                                  'aj':c[1],
                                                  'ak':c[2],
                                                  'al':c[3],
                                                  'func':c[4],
                                                  'c0':c[5],
                                                  'c1':c[6],
                                                  'c2':c[7],
                                                  'c3':c[8],
                                                  'c4':c[9],
                                                  'c5':c[10]}})
                            pcnt+=1
                        else:
                            impropers.update({icnt:{'ai':c[0],
                                                    'aj':c[1],
                                                    'ak':c[2],
                                                    'al':c[3],
                                                    'func':c[4],
                                                    'c0':c[5],
                                                    'c1':c[6],
                                                    'c2':c[7],
                                                    'c3':c[8],
                                                    'c4':c[9],
                                                    'c5':c[10]}})
                            icnt+=1
                    if len(c)>3 and c[4]=='1':
                        if dihflag%2:
                            propers.update({pcnt:{'ai':c[0],
                                                  'aj':c[1],
                                                  'ak':c[2],
                                                  'al':c[3],
                                                  'func':c[4],
                                                  'c0':c[5],
                                                  'c1':c[6],
                                                  'c2':c[7],
                                                  'c3':'',
                                                  'c4':'',
                                                  'c5':''}})
                            pcnt+=1
                        else:
                            impropers.update({icnt:{'ai':c[0],
                                                    'aj':c[1],
                                                    'ak':c[2],
                                                    'al':c[3],
                                                    'func':c[4],
                                                    'c0':c[5],
                                                    'c1':c[6],
                                                    'c2':c[7],
                                                    'c3':'',
                                                    'c4':'',
                                                    'c5':''}})
                            icnt+=1
                            
    return bonds,angles,propers,impropers
