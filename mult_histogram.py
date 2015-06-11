#make a scatter plot out of a list of files

import sys, os
import re, string
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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
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
      
colors=['k','r','g','b','m','c','y', 'lightblue']

usage = "Usage: %prog -f base_file -o output [options]\nLabels are taken from the relative path's first word"
parser = optparse.OptionParser(usage=usage)
#parser.add_option("-d", dest="d", default=".", help="directory to look for the files. default: current")
parser.add_option("-p", dest="param", help="mandatory parameter file for multiple xy scatter data")
parser.add_option("-o", dest="out", help="mandatory Output file for the scatterplot")
parser.add_option("-x", dest="xlab",  help="Label for x-axis. ")
parser.add_option("-y", dest="ylab", help="Label for y-axis. ")
parser.add_option("-m", dest="marker", default="." ,help="Marker symbol. default: \".\"")
parser.add_option("-s", dest="msize", default="2" ,help="Marker size. default: 2")
parser.add_option("-l", dest="linewidth", default="0" ,help="Linewidth. default: 0")
parser.add_option("--nr", dest="nrows", metavar="int", help="Use for a strict number of rows")
parser.add_option("--nc", dest="ncols", metavar="int", default="2" ,help="Default number of cols for given data will be 2")
parser.add_option("-a", dest="alpha", metavar="float", default=0.75 ,help="Alpha factor - for filling transparency. default: 0.75")
parser.add_option("--hist", dest="histtype", default="step" ,help="Histogram type. Default: step")
parser.add_option("-b", dest="bins", metavar="int", default=200 ,help="Number of bins. Default: 200")
parser.add_option("-n", dest="norm", metavar="int", default=2, help="Normalize histograms. 0=no/1=len(x)*dBin/2=len(x)")
parser.add_option("-L", dest="lower", metavar="float", help="Lower limit for the bins & x axis")
parser.add_option("-U", dest="upper", metavar="float", help="Upper limit for the bins & x axis")
parser.add_option("--axs", dest="axis_size", default="16" ,metavar="int",help="Axis label size. default: 16")
parser.add_option("--ts", dest="tick_size", default="14" ,metavar="int",help="Tick label size. default: 14")
parser.add_option("--fig_w", dest="fig_width", default=24.8 ,metavar="float",help="Figure width. default: 24.8")
parser.add_option("--fig_h", dest="fig_height", default=35.2 ,metavar="float",help="Figure height. default: 35.2")



mandatories = ["param", "out"]
(opt,args) = parser.parse_args()

for m in mandatories:
    if not opt.__dict__[m]:
        print "mandatory option is missing - %s\n"%m
        parser.print_help()
        exit(-1)
plots = []
legends = []
titles = []
Xmin = 1000000
Xmax = -1000000

if (opt.param[-3:] == "xvg"):
   plots.append(opt.param.split())
   legends.append("1".split()) 
   titles.append(opt.param[:-4].split())
else:
   plots, legends, titles = read_param(opt.param)

majorLocator   = MultipleLocator(0.1)
majorFormatter = FormatStrFormatter('%1.1f')
minorLocator   = MultipleLocator(3)

#if (len(plots) >8 ):
#   print "Too many figures per one page: "+len(plots)
#   exit()


print "Using the following files for the plots:"
plotN = 1
Ymax=0
max_rows = math.ceil(float(len(plots))/float(opt.ncols))
axes_list=[]
general_ax = ""
axes_array=[]

#iterate over number of plots
print len(plots)
fig=pl.figure(figsize=(float(opt.fig_width), float(opt.fig_height)), dpi=300)   #size of A4 page 
for p in (range(len(plots))):
   print "Plot "+str(plotN)
   print plots[p]
   print len(plots[p])
   if (len(plots[p])>8):
       print "Too many data series, current color map will fail";
       exit()
       
   ax1=fig.add_subplot(max_rows,opt.ncols,plotN, title=titles[p])
   ax1.set_title(titles[p], size=16)
   pl.subplots_adjust(hspace=0.5)
   axes_list.append(ax1)
   ax1.tick_params(labelsize=opt.tick_size)         
      
    #### Handle axes, first is defined, the others share its X axis.
   #if (general_ax is "" and len(plots)>1):
   #   general_ax=fig.add_subplot(max_rows,1,p+1)
   #   ax1=general_ax      
   #else:
   #   ax1=fig.add_subplot(max_rows,1,p+1), sharex=general_ax)
   
   
   #iterate over number of sets shown on one plot
   for i in (range(len(plots[p]))):      
      data=[]
      #lab = (lst[i].split("_"))[0].split("/")[0]
      if (len(legends) > 0):
         print str(plots[p][i]) +" label: "+str(legends[p][i])
      else:
         print str(plots[p][i])         
     
      #### Hide tick labels from second histogram and on
      #if (i!=(len(plots[p])-1)):
      #   pl.setp( ax1.get_xticklabels(), visible=False)

      #### load data and define x limits
      data = load_xvg_Data(plots[p][i]) #already replaces @ and & to comments
     
      print len(data[:,1])
      if (opt.lower == None and Xmin>data[:,1].min() ):
         Xmin = data[:,1].min()
         #opt.lower = data[:,1].min()
      if (opt.lower != None):
         Xmin = float(opt.lower)
      if (opt.upper == None and Xmax<data[:,1].max() ):
         Xmax = data[:,1].max()
         #opt.upper = data[:,1].max()   
      if (opt.upper != None):
         Xmax = float(opt.upper)
      print Xmin, Xmax
      #### histogram and normalize according to: 2=use 1/data length as weights 
      if (opt.norm == 2):
         w = pl.ones(len(data[:,1]))*100/len(data[:,1]) 
         try:
            pl.hist(data[:,1], int(opt.bins), range=(Xmin, Xmax), weights=w ,histtype=opt.histtype, facecolor=colors[i], alpha=opt.alpha, label=legends[p][i])       
         except ValueError: 
            print len(data[:,1]), int(opt.bins), Xmin, Xmax, colors[i], legends[p][i]
      else:
         pl.hist(data[:,1], int(opt.bins), range=(Xmin, Xmax), normed=opt.norm, histtype=opt.histtype, facecolor=colors[i], alpha=opt.alpha, label=legends[p][i])     
      
      #### Save the max. Y value over all histograms, so an equal scale can be used for all Y axes
      if (Ymax < pl.gca().get_ylim()[1]):
         Ymax=pl.gca().get_ylim()[1]
      pl.grid(True)   
      pl.legend(loc=1, numpoints=1, markerscale=4, fontsize=18) #location=1=upper right 
      pl.xlabel(opt.xlab, fontsize=opt.axis_size)   
      pl.ylabel(opt.ylab, fontsize=opt.axis_size)   
      
   
   

  #ax1.xaxis.set_major_locator(majorLocator)
  # ax1.xaxis.set_major_formatter(majorFormatter)  
  # ax1.xaxis.set_minor_locator(minorLocator)
   
   pl.grid(True)         
   ax1.tick_params(labelsize=opt.tick_size)              
   plotN = plotN+1
#### Change the scale for all Y axes to the largest      
#for r in (range(len(plots))):
#   pl.sca(axes_list[r])
#   pl.ylim(0,Ymax)


#pl.legend(loc='best', numpoints=1, markerscale=18, fontsize=24)  
               
#pl.grid(True)
pl.subplots_adjust(hspace=0.25)
pl.savefig(opt.out)

