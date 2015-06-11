#make a scatter plot out of a list of files

import sys, os
import re
import math
from glob import glob
import subprocess
import numpy as np
import pyutil
import pylab as pl
import optparse
from DavidLib import load_xvg_Data, smooth
import matplotlib.cm as cm
import matplotlib.colors as clrs

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
      
#colors=['k','r','g','b','m','c','y', 'lightblue']
colors=['g','r','k','b','m','c','y', 'lightblue']

usage = "Usage: %prog -f base_file -o output [options]\nLabels are taken from the relative path's first word"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-f", dest="file", default=None, help="one data file to process instead of listing them in a .param file")
parser.add_option("-p", dest="param", default=None, help="mandatory parameter file for multiple xy scatter data")
parser.add_option("-o", dest="out", help="mandatory Output file for the scatterplot")
parser.add_option("-x", dest="xlab", default="ev1" , help="Label for x-axis. default: \"ev1\"")
parser.add_option("-y", dest="ylab", default="ev2" ,help="Label for y-axis. default: \"ev2\"")
parser.add_option("-m", dest="marker", default="." ,help="Marker symbol. default: \".\"")
parser.add_option("-s", dest="msize", default="2" ,help="Marker size. default: 2")
parser.add_option("-l", dest="linewidth", default="0" ,help="Linewidth. default: 0")
parser.add_option("--nr", dest="nrows", metavar="int", help="Use for a strict number of rows")
parser.add_option("--nc", dest="ncols", metavar="int", default="2" ,help="Default number of cols for given data will be 2")
parser.add_option("--axs", dest="axis_size", default="20" ,metavar="int",help="Axis label size. default: 20")
parser.add_option("--ts", dest="tick_size", default="20" ,metavar="int",help="Tick label size. default: 20")
parser.add_option("--fig_w", dest="fig_width", default=24.8 ,metavar="float",help="Figure width. default: 24.8")
parser.add_option("--fig_h", dest="fig_height", default=35.2 ,metavar="float",help="Figure height. default: 35.2")
parser.add_option("--color", dest="color", action="store_true", help="color scatter by 3rd col value")
parser.add_option("--nolegend", dest="nolegend", action="store_true", help="Do not Display legend")
parser.add_option("--onelegend", dest="onelegend", action="store_true", help="Display legend once")
parser.add_option("--xmin", dest="xmin", help="x axis min. limit")
parser.add_option("--xmax", dest="xmax", help="x axis max. limit")
parser.add_option("--ymin", dest="ymin", help="y axis min. limit")
parser.add_option("--ymax", dest="ymax", help="y axis max. limit")
parser.add_option("--lmin", dest="lmin", help="color bar min. limit")
parser.add_option("--lmax", dest="lmax", help="color bar max. limit")

mandatories = ["out"]
(opt,args) = parser.parse_args()

for m in mandatories:
    if not opt.__dict__[m]:
        print "mandatory option is missing - %s\n"%m
        parser.print_help()
        exit(-1)
plots = []
legends = []
titles = []
if (opt.param == None): 
   if (opt.file == None):
      print "Error: Need an input file, either one file with -f, or multiple described in -p file.param"
      exit()
   if (not os.path.exists(opt.file)):
      print "Error: File in (-f) option (%s) does not exist. "%opt.file
      exit()
   plots.append([opt.file])
   legends.append([""])
   titles=opt.file[:-4]
else: 
   if (opt.param[-3:] == "xvg"):
      plots.append(opt.param.split())
      legends.append("1".split()) 
      titles.append(opt.param[:-4].split())
   else:
      plots, legends, titles = read_param(opt.param)


#if (len(plots) >8 ):
#   print "Too many figures per one page: "+len(plots)
#   exit()

fig=pl.figure(figsize=(float(opt.fig_width), float(opt.fig_height)), dpi=300)   #size of A4 page in pixels for 300dpi
print "Using the following %d file(s) for the plots:" %(len(plots))
plotN = 1
max_rows = math.ceil(float(len(plots))/float(opt.ncols))
#iterate over number of plots

for p in (range(len(plots))):
#for lst in plots:
   print "Plot "+str(plotN)
   print plots[p]
   print len(plots[p])
   if (len(plots[p])>8):
      print "Too many data series, current color map will fail";
      exit()
   if (len(plots[p]) > 2 and opt.color):
      print "Scatter will be colored according to 3rd coloumn, use only one series per plot"      
      exit(-1)
      
   #ax1=fig.add_subplot(max_rows,opt.ncols,plotN, title=titles[p])
   #ax1.set_title(titles[p], size=16)
   pl.subplots_adjust(hspace=0.5)
   ax1=pl.subplot(max_rows,opt.ncols,plotN, title=titles[p])
   ax1.set_title(titles[p], size=opt.axis_size)

   #iterate over number of sets shown on one plot
   for i in (range(len(plots[p]))):
      data=[]
      #lab = (lst[i].split("_"))[0].split("/")[0]
      if (len(legends) > 0):
         print str(plots[p][i]) +" label: "+str(legends[p][i])
      else:
         print str(plots[p][i])
      data = load_xvg_Data(plots[p][i]) #already replaces @ and & to comments
      if (opt.color):
         cmap=cm.rainbow 
         vmax = max(data[:,2])
         vmin = min(data[:,2])        
         if (opt.lmax):
            vmax = opt.lmax
         if (opt.lmin):
            vmin = opt.lmin
         if ((not opt.onelegend) and (not opt.nolegend)):
            result = pl.scatter(data[:,0],data[:,1], c=data[:,2], s=float(opt.msize), linewidth=int(opt.linewidth), cmap=cmap, vmin=vmin, vmax=vmax,label=legends[p][i])      
         else:
            result = pl.scatter(data[:,0],data[:,1], c=data[:,2], s=float(opt.msize), linewidth=int(opt.linewidth), cmap=cmap, vmin=vmin, vmax=vmax)      
         cbar = pl.colorbar(result, shrink=0.9)
                #cbar = pl.colorbar(result, shrink=0.7, ticks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
                # cbar.ax.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
                #axc, kw = mp.colorbar.make_axes(ax1, shrink=0.8)
                #cbar = pl.colorbar.Colorbar(axc, result)
         pl.legend(loc='upper right',  scatterpoints=1, fontsize=20)
      else:    
         #ax1.scatter(data[:,0],data[:,1], s=1, c=colors[i],marker="o", linewidths=None)
         if (len(legends) > 0):
            if ( i == 3):
               opt.msize = 8
            pl.plot(data[:,0],data[:,1], marker=opt.marker, c=colors[i], markersize=float(opt.msize), linewidth=int(opt.linewidth),label=legends[p][i])
         else:
            pl.plot(data[:,0],data[:,1], marker=opt.marker, c=colors[i], markersize=float(opt.msize), linewidth=int(opt.linewidth))               
      pl.grid(True)         
      ax1.tick_params(labelsize=opt.tick_size)      
      if (opt.xmin):
         pl.xlim(xmin=float(opt.xmin))
      if (opt.xmax):         
         pl.xlim(xmax=float(opt.xmax))
      if (opt.ymin):
         pl.ylim(ymin=float(opt.ymin))
      if (opt.ymax):
         pl.ylim(ymax=float(opt.ymax))

   plotN = plotN+1
   #pl.xlabel(opt.xlab, fontsize=opt.axis_size)   
   #pl.ylabel(opt.ylab, fontsize=opt.axis_size)   
   if ((not opt.color) and (len(plots[p]) > 1)  and not opt.onelegend and not opt.nolegend):
      pl.legend(loc='best', numpoints=1, markerscale=2)  
      

if (opt.onelegend):
   pl.legend(loc='best', numpoints=1, markerscale=2)        
# Set common labels
fig.text(0.5, 0.05, opt.xlab, ha='center', va='center', fontsize=opt.axis_size)
fig.text(0.06, 0.5, opt.ylab, ha='center', va='center', rotation='vertical', fontsize=opt.axis_size)

       
#pl.grid(True)
pl.subplots_adjust(hspace=0.35)
pl.savefig(opt.out)
 
#pyutil.scatter_plot_xy("2d_1-2_%s.xvg" %name , "2d_1-2_%s.pdf" %name, "ev1","ev2")
#pyutil.scatter_plot_xy("2d_1-3_%s.xvg" %name,  "2d_1-3_%s.pdf" %name, "ev1", "ev3")
