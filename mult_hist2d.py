from math import *
import pylab as p
import matplotlib.pyplot as plt
import numpy as np
import sys, os
import pyutil
import matplotlib.cm as cm

def read_param(pfile):
   plots = []
   legends = []
   titles = []
   lines = open(pfile).readlines()
   for line in lines:
      if (line.strip()):   
         parts = line.split("=")
         parts[0] = parts[0].strip()
         parts[1] = parts[1].strip()
         if ('files_p' in parts[0]):
            plots.append(parts[1].split())
         if ('legend_p' in parts[0]):  
            legends.append(parts[1].split())
         if ('titles' in parts[0]):
            titles=parts[1].split()
   return plots, legends, titles

usage = "Usage: %prog -f base_file -o output [options]\nLabels are taken from the relative path's first word"
parser = optparse.OptionParser(usage=usage)
#parser.add_option("-d", dest="d", default=".", help="directory to look for the files. default: current")
#parser.add_option("-f", dest="base", help="mandatory base name for xy scatter data")
parser.add_option("-p", dest="param", help="mandatory parameter file for multiple xy scatter data")
parser.add_option("-o", dest="out", help="mandatory Output file for the scatterplot")
parser.add_option("-x", dest="xlab", default="" , help="Label for x-axis. default: \"\"")
parser.add_option("-y", dest="ylab", default="" ,help="Label for y-axis. default: \"\"")
parser.add_option("--nr", dest="nrows", metavar="int", help="Use for a strict number of rows")
parser.add_option("--nc", dest="ncols", metavar="int", default="2" ,help="Default number of cols for given data will be 2")
parser.add_option("--axs", dest="axis_size", default="16" ,metavar="int",help="Axis label size. default: 16")
parser.add_option("--ts", dest="tick_size", default="14" ,metavar="int",help="Tick label size. default: 14")
parser.add_option("--fig_w", dest="fig_width", default=24.8 ,metavar="float",help="Figure width. default: 24.8")
parser.add_option("--fig_h", dest="fig_height", default=35.2 ,metavar="float",help="Figure height. default: 35.2")
parser.add_option("--xmin", dest="xmin", help="x axis min. limit")
parser.add_option("--xmax", dest="xmax", help="x axis max. limit")
parser.add_option("--ymin", dest="ymin", help="y axis min. limit")
parser.add_option("--ymax", dest="ymax", help="y axis max. limit")
parser.add_option("--lmin", dest="lmin", help="color bar min. limit")
parser.add_option("--lmax", dest="lmax", help="color bar max. limit")
parser.add_option("--cm", dest="cmap", default="Paired", help="colormap")
mandatories = ["param", "out"]
(opt,args) = parser.parse_args()

for m in mandatories:
    if not opt.__dict__[m]:
        print "mandatory option is missing - %s\n"%m
        parser.print_help()
        exit(-1)

plots, legends, titles = read_param(opt.param)

fig=pl.figure(facecolor='w', figsize=(float(opt.fig_width), float(opt.fig_height)), dpi=300)   #size of A4 page in pixels for 300dpi
print "Using the following files for the plots:"
plotN = 1


data=pyutil.load_xvg_Data(infile)
#range=[[-100, 100],[-100,100]],
H, xedges, yedges = np.histogram2d(data[:,1], data[:,0], bins=(300, 300))

extent = [yedges[0], 1, xedges[0], xedges[-1]]
print extent

#extent = [-500, 500, -100, 100]
#plt.imshow(H, extent=extent, interpolation='nearest', cmap=cm.Paired , figure=p.figure(figsize=(14,10), dpi=300) )
plt.imshow(H, extent=extent, interpolation='nearest', cmap=cm.Paired , figure=fig )
plt.colorbar()
#plt.xlabel('PA_O3A_PB_Mg dihedral')
#plt.ylabel('PG_O3B_PB_Mg dihedral')
plt.savefig(out)
