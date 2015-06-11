#!/usr/bin/env python
#AUTHOR: Martin Hoefling and Andrea Vaiana
#EMAIL: mhoefli@gwdg.de avaiana@gwdg.de
#DATE: 29. Jun 2010
#ARCH: any
#REQUIRES: python 
#STATUSCOMMENTS: Seems to work.
#DESCRIPTION: Allows calculation of type D vsites (Gromacs 4.0 manual pp. 128), e.g. CH3 or NH3 groups.
#USAGE: Just fire it up without parameters :-)
#EXAMPLE: see usage...

###

from __future__ import division
from math import pi,cos,sqrt
import sys

debug=False

if len(sys.argv)!=7:
    print "Usage: %s ffanythingbon.itp CB CT HC MCH3A 12"%sys.argv[0]
    print "    First (1) atom is the connector, the second (2) connects to three hydrogen (3) atoms."
    print "    (4) is the type of virtual site (usually depends on the angle)"
    print "    (5) is the mass of the hydrogen connecting atom: C:12 N:14"
    sys.exit(1)


filename=sys.argv[1]

connector=sys.argv[2]
metcar=sys.argv[3]
hydrogen=sys.argv[4]
consttype=sys.argv[5]
try:
    mymass=float(sys.argv[6])
except ValueError:
    print "Please insert a correct mass"
    sys.exit(1)

print connector, metcar, hydrogen
print "Reading in parameters..."

fh=open(filename)

inbonds=False
inangles=False
inconstraints=False
inheavy=False

alpha=False
cthc=False
ctcx=False
vsvs=False

for line in fh:
    if line.startswith("["):
        secname=line.split()[1]
        if secname=="bondtypes":
            inbonds=True
            inangles=False
            inconstraints=False
    
        elif secname=="angletypes":
            inangles=True
            inbonds=False
            inconstraints=False
    
        elif secname=="constrainttypes":
            inconstraints=True
            inbonds=False
            inangles=False
    
        else:
            inbonds=False
            inangles=False
            inconstraints=False


    elif line.strip()=="":
        continue
    
    elif line.startswith("#ifdef HEAVY_H"):
        inheavy=True

    elif line.startswith("#else") or line.startswith("#endif"):
        inheavy=False
    
    elif inheavy:
        continue
    
    elif line.strip().startswith(";"):
        continue

    else:
        if inbonds:
            spl=line.split()
            name1=spl[0]
            name2=spl[1]

            if (name1==connector and name2==metcar) or (name1==metcar and name2==connector):
                ctcx=float(spl[3])
                print "Bondlength between %s and %s is %6.4f"%(name1,name2,ctcx)
    
            elif (name1==hydrogen and name2==metcar) or (name1==metcar and name2==hydrogen):
                cthc=float(spl[3])
                print "Bondlength between %s and %s is %6.4f"%(name1,name2,cthc)

        if inangles:
            spl=line.split()
            name1=spl[0]
            name2=spl[1]
            name3=spl[2]

            if name2==metcar:
                if (name1==connector and name3==hydrogen) or (name1==hydrogen and name3==connector):
                    alpha=float(spl[4])
                    print "Angle between %s,%s and %s is %6.2f"%(name1,name2,name3,alpha)

        if inconstraints:
            spl=line.split()
     #       print line
            name1=spl[0]
            name2=spl[1]
            if name1==consttype and name2==consttype:
                vsvs=float(spl[3])
                print "Constraint between constrainttype %s is %6.4f"%(name1,vsvs)



if not alpha:
    print "couldn't read angle"
    sys.exit(1)
    
if not cthc:
    print "couldn't read in distance from \"carbon to hydrogen\""
    sys.exit(1)

if not ctcx:
    print "couldn't read in distance from \"connector atom to carbon\""
    sys.exit(1)

if not vsvs:
    print "coudn't read in v-site distance"
    sys.exit(1)

alphaprime=180.-alpha
cta=cos(alphaprime*pi/180.)*cthc
ctb=(3.024/(mymass+3.024))*cta
ctpb=ctb+ctcx
ctpvs=sqrt(ctpb**2+(vsvs/2)**2)

if debug:
    print "alpha'",alphaprime
    print "CT-A",cta
    print "CT-B",ctb
    print "CX-B",ctpb
print
print "Put the following in your [ constrainttypes ] section:"
print "########################"
print
print "%s %s 2 %8.6f"%(consttype,connector,ctpvs)
print
print "########################"


