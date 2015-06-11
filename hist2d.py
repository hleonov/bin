from math import *
import pylab as p
import matplotlib.pyplot as plt
import numpy as np
import sys, os
import pyutil
import matplotlib.cm as cm

infile = sys.argv[1]
out    = sys.argv[2]
labelx = sys.argv[3]
labely = sys.argv[4]
title =  sys.argv[5]
data=pyutil.load_xvg_Data(infile)


H, yedges, xedges = np.histogram2d(data[:,1], data[:,0], range=[[0, 1],[0,1]], bins=(300, 300))

extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
print extent

#cutoff=50

#T = H > cutoff
#H = H - (H * T) +(cutoff * T)

plt.imshow(H, origin='lower',extent=extent, interpolation='nearest', cmap=cm.Paired , figure=p.figure(facecolor='w',figsize=(14,10), dpi=300) )
plt.colorbar()
plt.title(title, size=20)
plt.xlabel(labelx)
plt.ylabel(labely)
plt.savefig(out)
