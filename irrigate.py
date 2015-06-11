#!/usr/bin/python

from sys import argv
from random import random
execfile("/home/daniel/proj/hook/src/pdb.py")

infile = argv[1]
outfile = argv[2]
resol = float(argv[3])

protein = Helix()
protein.load(infile)

def frange(fr, to, step):
	l = int((to-fr)/step)+1
	return [fr+x*step for x in range(l)]

prot_only = [x for x in protein if x.aa not in ['SOL ', 'DMP ']]
minx = min([x.atom.x for x in prot_only])
miny = min([x.atom.y for x in prot_only])
minz = min([x.atom.z for x in prot_only])
maxx = max([x.atom.x for x in prot_only])
maxy = max([x.atom.y for x in prot_only])
maxz = max([x.atom.z for x in prot_only])

x_vals = [x for x in frange(minx + resol, maxx - resol, resol)]
y_vals = [y for y in frange(miny + resol, maxy - resol, resol)]
z_vals = [z for z in frange(minz + resol, maxz - resol, resol)]

res = []
print "Now in main loop (", len(x_vals), " iterations)"
for i in x_vals:
	print str(x_vals.index(i))
	for j in y_vals:
		for k in z_vals:
			add = True
			mindist = 99999
			for a in protein:
				c = a.atom.c
				if i<=c[0]<=i+resol and j<=c[1]<=j+resol and k<=c[2]<=k+resol:
					add=False
					break
			if add:
				for a in protein:
					if a.aa == 'DMP ':
						d = dist([i+.5*resol,j+.5*resol,k+.5*resol], a.atom.c)
						if d < 2* resol:
							add = False
						#	print "Declined ", i, j, k, " : Near Lipid."	 
							break
					if a.aa == 'SOL ':
						d = dist([i+.5*resol,j+.5*resol,k+.5*resol], a.atom.c)
						if d < resol:
							add = False
						#	print "Declined ", i, j, k, " : Near Solvent."
							break
			if add:
				res.append([i,j,k])
				
#Choose only those results which have a neighbor in the results
goodres = []
print "Calculating neighbors"
for r in res:
	d = [dist(r, x) for x in res]
	close = [x for x in d if x <= sqrt(3)*resol]
	if len(close) > 2:
		goodres.append(r)
	#else:
	#	print "Declined ", r[0], r[1], r[2], " : Not enough neighbors."

h2o = Helix()
h2o.copy(protein[-3:])
print "adding " + str(len(goodres)) + " water molecules"
for i in goodres:
    orix = i[0] + resol/3 + resol/3 * random()
    oriy = i[1] + resol/3 + resol/3 * random()
    oriz = i[2] + resol/3 + resol/3 * random()			  
    rx   = 2 * pi * random()
    ry   = 2 * pi * random()
    rz   = 2 * pi * random()
    old_O = map(None, h2o[0].atom.c)
    h2o.locate_origin(old_O)
    h2o.rotx(rx) ; h2o.roty(ry) ; h2o.rotz(rz)
    h2o.locate_origin([-orix, -oriy, -oriz])
    for j in h2o:
    	j.res += 1
    protein.append(h2o)
    
protein.renum()
protein.save(outfile)

print "Program terminated successfully."
