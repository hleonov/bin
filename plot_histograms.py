#make a plot of histograms out of a list of files

import sys, os
import re, string
import math
from glob import glob
import subprocess
import numpy
import pyutil
import pylab as pl
import optparse
from DavidLib import load_xvg_Data, smooth
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
      
colors=['k','r','g','b','m','c','y']

usage = "Usage: %prog -f base_file -o output [options]\nLabels are taken from the relative path's first word"
parser = optparse.OptionParser(usage=usage)
#parser.add_option("-d", dest="d", default=".", help="directory to look for the files. default: current")
#parser.add_option("-f", dest="base", help="mandatory base name for xy scatter data")
parser.add_option("-p", dest="param", help="mandatory parameter file for xy data")
parser.add_option("-o", dest="out", help="mandatory Output file for the histogram")
parser.add_option("-x", dest="xlab", default="Distance (nm)" , help="Label for x-axis. default: \"Distance(nm)\"")
parser.add_option("-y", dest="ylab", default="Count" ,help="Label for y-axis. default: \"Count\"")
parser.add_option("-a", dest="alpha", metavar="float", default=0.75 ,help="Alpha factor - for filling transparency. default: 0.75")
parser.add_option("--hist", dest="histtype", default="stepfilled" ,help="Histogram type. Default: stepfilled")
parser.add_option("-b", dest="bins", metavar="int", default=200 ,help="Number of bins. Default: 200")
parser.add_option("-n", dest="norm", metavar="int", default=2, help="Normalize histograms. 0=no/1=len(x)*dBin/2=len(x)")
parser.add_option("-L", dest="lower", metavar="float", help="Lower limit for the bins & x axis")
parser.add_option("-U", dest="upper", metavar="float", help="Upper limit for the bins & x axis")
parser.add_option("--nc", dest="ncols", metavar="int", default="1" ,help="Default number of cols for given data will be 1")
parser.add_option("--axs", dest="axis_size", default="14" ,metavar="int",help="Axis label size. default: 14")
parser.add_option("--ts", dest="tick_size", default="10" ,metavar="int",help="Tick label size. default: 10")
parser.add_option("--fig_w", dest="fig_width", default=24.8 ,metavar="float",help="Figure width. default: 24.8")
parser.add_option("--fig_h", dest="fig_height", default=35.2 ,metavar="float",help="Figure height. default: 35.2")
#parser.add_option("-l", dest="linewidth", default="0" ,help="Linewidth. default: 0")
mandatories = ["param", "out"]
(opt,args) = parser.parse_args()
print opt
for m in mandatories:
    if not opt.__dict__[m]:
        print "mandatory option is missing - %s\n"%m
        parser.print_help()
        exit(-1)

plots, legends, titles = read_param(opt.param)

#if (len(plots) >8 ):
#   print "Too many figures per one page: "+len(plots)
#   exit()
majorLocator   = MultipleLocator(0.1)
#majorLocator   = MultipleLocator(1)
majorFormatter = FormatStrFormatter('%1.1f')
minorLocator   = MultipleLocator(3)

print "Using the following files for the plots:"
plotN = 1
Ymax=0
max_rows = math.ceil(float(len(plots))/float(opt.ncols))
axes_list=[]
#iterate over number of figures (each would contain multiple histograms)
for p in (range(len(plots))):
   fig=pl.figure(figsize=(float(opt.fig_width), float(opt.fig_height)), dpi=300)   #size of A4 page in pixels for 300dpi
   fig.suptitle(string.join( titles[p].split("_"), " "))
   print "Plot %d %s" %(plotN, plots[p])
   max_rows = len(plots[p])
   print "num of hists:" +str(max_rows)
   if (max_rows>7):
      print "Too many data series, current color map will fail";
      exit()
   
   #iterate over number of sets shown on one plot
   general_ax = ""
   axes_array=[]
   for i in (range(len(plots[p]))):
      data=[]
      #lab = (lst[i].split("_"))[0].split("/")[0]
      print str(plots[p][i])+" "+str(i) +" label: "+str(legends[p][i])

      #### Handle axes, first is defined, the others share its X axis.
      if (general_ax is ""):
         general_ax=fig.add_subplot(max_rows,1,i+1)
         ax1=general_ax      
      else:
         ax1=fig.add_subplot(max_rows,1,i+1, sharex=general_ax)
      axes_list.append(ax1)
      ax1.tick_params(labelsize=opt.tick_size)         

      #### Hide tick labels from second histogram and on
      if (i!=(len(plots[p])-1)):
         pl.setp( ax1.get_xticklabels(), visible=False)

      #### load data and define x limits
      data = load_xvg_Data(plots[p][i]) #already replaces @ and & to comments   
      print len(data[:,1])
      if (opt.lower == None):
         opt.lower = data[:,1].min()
      if (opt.upper == None):
         opt.upper = data[:,1].max()   
      opt.lower = float(opt.lower)
      opt.upper = float(opt.upper)
      
      #### histogram and normalize according to: 2=use 1/data length as weights 
      if (opt.norm == 2):
         w = pl.ones(len(data[:,1]))*100/len(data[:,1]) 
#         try:
         pl.hist(data[:,1], int(opt.bins), range=(opt.lower, opt.upper), weights=w ,histtype=opt.histtype, facecolor=colors[i], alpha=opt.alpha, label=legends[p][i])       
 #        except ValueError: 
 #           print len(data[:,1]), int(opt.bins), opt.lower, opt.upper, colors[i], legends[p][i]
      else:
         pl.hist(data[:,1], int(opt.bins), range=(opt.lower, opt.upper), normed=opt.norm, histtype=opt.histtype, facecolor=colors[i], alpha=opt.alpha, label=legends[p][i]) 
      
      #### Save the max. Y value over all histograms, so an equal scale can be used for all Y axes
      if (Ymax < pl.gca().get_ylim()[1]):
         Ymax=pl.gca().get_ylim()[1]
      pl.grid(True)   
      pl.legend(loc=1, numpoints=1, markerscale=4) #location=1=upper right
      
   #### Change the scale for all Y axes to the largest   
   for i in (range(len(plots[p]))):
      pl.sca(axes_list[i])
      pl.ylim(0,Ymax)
 
   ax1.xaxis.set_major_locator(majorLocator)
   ax1.xaxis.set_major_formatter(majorFormatter)  
   ax1.xaxis.set_minor_locator(minorLocator)
    
   pl.xlim(opt.lower, opt.upper)
   pl.xlabel(opt.xlab, fontsize=opt.axis_size)   
   #pl.ylabel(opt.ylab, fontsize=opt.axis_size)
 
#pl.grid(True)
   pl.subplots_adjust(hspace=0.2)
#   pl.savefig("fig"+str(plotN)+"_"+opt.out)
   pl.savefig(opt.out)
   plotN = plotN+1
