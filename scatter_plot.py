#make a scatter plot out of a list of files

import sys, os
import re
from glob import glob
import subprocess
import numpy
import pyutil
import pylab as pl
import optparse
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from DavidLib import load_xvg_Data, smooth

colors=['k', 'g','r','b','m','c','y'] #put 'k' at front
#lab = ["WT", "D19", "R19", "A19"]

usage = "Usage: %prog -f base_file -o output [options]\nLabels are taken from the relative path's first word"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-d", dest="d", default=".", help="directory to look for the files. default: current")
parser.add_option("-f", dest="base", help="mandatory base name for xy scatter data. If exists, 3rd col will determine size of datapoint (e.g for frequency)")
parser.add_option("-o", dest="out", help="mandatory Output file for the scatterplot")
parser.add_option("-x", dest="xlab", default="ev1" , help="Label for x-axis. default: \"ev1\"")
parser.add_option("-y", dest="ylab", default="ev2" ,help="Label for y-axis. default: \"ev2\"")
parser.add_option("-m", dest="marker", default="." ,help="Marker symbol. default: \".\"")
parser.add_option("-s", dest="msize", default="2" ,help="Marker size. default: 2")
parser.add_option("-g", dest="lmsize", default=6 ,help="Legend Marker size. default: 6")
parser.add_option("-l", dest="linewidth", default="0" ,help="Linewidth. default: 0")
parser.add_option("--ts", dest="tick_size", default="10" ,metavar="int",help="Tick label size. default: 10")
parser.add_option("--tick", dest="tick_space", help="New definiton for major tick spacing");
parser.add_option("--alpha", dest="alpha", default="1.0", help="make color transparent. range: 0.0-1.0, default: 1.0 (opaque)")
parser.add_option("--mult", dest="mult", default=100, help="Multiply 3rd column by this value for visible data")

mandatories = ["base", "out"]
(opt,args) = parser.parse_args()

for m in mandatories:
    if not opt.__dict__[m]:
        print "mandatory option is missing - %s\n"%m
        parser.print_help()
        exit(-1)

print opt.d, opt.base

if (opt.d is not "."):
   lst = glob("%s/%s" %(opt.d, opt.base))
else:
   lst = glob("*%s*" %opt.base)
lst.sort()
print lst
if (len(lst)>7):
   print "Too many data series, current color map will fail";
   exit()

print "Using these files for the plot:"
fig=pl.figure(figsize=(14, 10))
ax1=fig.add_subplot(111)
for i in (range(len(lst))):
   data=[]
   lab = (lst[i].split("_"))[0].split("/")[0]
   print lst[i]+" label: "+lab+" color "+colors[i]
   data = load_xvg_Data(lst[i]) #already replaces @ and & to comments
   #pl.plot(data[:,0],data[:,1], c=colors[i], marker=opt.marker, markersize=int(opt.msize), linewidth=int(opt.linewidth), label=lab)
   try: 
   #opt.mult*data[:,2] 
      ax1.scatter(data[:,0],data[:,1], s=opt.mult*data[:,2], c=colors[i],facecolors='none', edgecolors=colors[i], marker="o",linewidth=None, alpha=float(opt.alpha), label=lab)
      print "Third coloumn will be used for marker size"
   except IndexError:
      print "Normal XY scatter, no third coloumn"
      #ax1.scatter(data[:,0],data[:,1], s=opt.msize, c=colors[i],marker="o", linewidth=None, alpha=float(opt.alpha), label=lab)
      pl.plot(data[:,0],data[:,1], c=colors[i], marker=opt.marker, markersize=int(opt.msize), linewidth=int(opt.linewidth), label=lab)

ax1.grid(True, which='both')
if (opt.tick_space is not None):
   print opt.tick_space
   ax1.xaxis.set_major_locator(MultipleLocator(int(opt.tick_space)))
   ax1.yaxis.set_major_locator(MultipleLocator(int(opt.tick_space)))

ax1.tick_params(which='both', labelsize=opt.tick_size)
pl.legend(loc='best', scatterpoints=1, markerscale=int(opt.lmsize), numpoints=1)
#pl.plot([0,0.4],[0,0.4], c="blue")
#pl.xlim(xmin=0)
#pl.ylim(ymin=0)
ax1.tick_params(labelsize=20)
pl.xlabel(opt.xlab, size=20)   
pl.ylabel(opt.ylab, size=20)
pl.savefig(opt.out)
 
#pyutil.scatter_plot_xy("2d_1-2_%s.xvg" %name , "2d_1-2_%s.pdf" %name, "ev1","ev2")
#pyutil.scatter_plot_xy("2d_1-3_%s.xvg" %name,  "2d_1-3_%s.pdf" %name, "ev1", "ev3")
