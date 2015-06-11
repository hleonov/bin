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
import matplotlib.gridspec as gridspec
import matplotlib.colorbar
#from scipy.optimize import curve_fit
from scipy import stats

#def func(x, a, b, c):
#    return a * np.exp(-b * x) + c
def func(x, a, b, c, d):
    return a*np.exp(-c*(x-b))+d  
    
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
parser.add_option("-x", dest="xlab", default="ev1" , help="Label for x-axis. default: \"ev1\"")
parser.add_option("-y", dest="ylab", default="ev2" ,help="Label for y-axis. default: \"ev2\"")
parser.add_option("-m", dest="marker", default="." ,help="Marker symbol. default: \".\"")
parser.add_option("-s", dest="msize", default="2" ,help="Marker size. default: 2")
parser.add_option("-l", dest="linewidth", default="0" ,help="Linewidth. default: 0")
parser.add_option("--nr", dest="nrows", metavar=int, help="Use for a strict number of rows")
parser.add_option("--nc", dest="ncols", metavar=int, default=2 ,help="Default number of cols for given data will be 2")
parser.add_option("--axs", dest="axis_size", default="20" ,metavar="int",help="Axis label size. default: 20")
parser.add_option("--ts", dest="tick_size", default="20" ,metavar="int",help="Tick label size. default: 20")
parser.add_option("--fig_w", dest="fig_width", default=24.8 ,metavar="float",help="Figure width. default: 24.8")
parser.add_option("--fig_h", dest="fig_height", default=35.2 ,metavar="float",help="Figure height. default: 35.2")
parser.add_option("--color", dest="color", action="store_true", help="color scatter by 3rd col value")
parser.add_option("--nolegend", dest="nolegend", action="store_true", help="Do not Display legend")
parser.add_option("--onelegend", dest="onelegend", action="store_true", help="Display legend once")
parser.add_option("--reg", dest="linearReg", action="store_true", help="Perform linear regression on scatter data and plot it")
parser.add_option("--xmin", dest="xmin", help="x axis min. limit")
parser.add_option("--xmax", dest="xmax", help="x axis max. limit")
parser.add_option("--ymin", dest="ymin", help="y axis min. limit")
parser.add_option("--ymax", dest="ymax", help="y axis max. limit")
parser.add_option("--lmin", dest="lmin", help="color bar min. limit")
parser.add_option("--lmax", dest="lmax", help="color bar max. limit")
parser.add_option("--hist", dest="hist", default=None, help="Display histogram: None, histX, histY, histXY")

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
if (opt.param[-3:] == "xvg"):
   plots.append(opt.param.split())
   legends.append("1".split()) 
   titles.append(opt.param[:-4].split())
else:
   plots, legends, titles = read_param(opt.param)

reg_output_fp = open("linreg_output.dat",'w')

#if (len(plots) >8 ):
#   print "Too many figures per one page: "+len(plots)
#   exit()

fig=pl.figure(figsize=(float(opt.fig_width), float(opt.fig_height)), dpi=300)   #size of A4 page in pixels for 300dpi
print "Using the following files for the plots:"
plotN = 1
max_rows = int(math.ceil(float(len(plots))/float(opt.ncols)))
print "max_rows: ", max_rows, " max_cols: ", opt.ncols
#iterate over number of plots
print len(plots)
rowN = colN = 0
for p in (range(len(plots))):
   print "Plot "+str(plotN)
   print plots[p]
   print len(plots[p])
   if (len(plots[p])>8):
      print "Too many data series, current color map will fail";
      exit()
   if (len(plots[p]) > 2 and opt.color):
      print "Scatter will be colored according to 3rd coloumn, use only one series per plot"      
      exit(-1)
   
   rowN = p / int(opt.ncols)
   colN = p % int(opt.ncols)
#   print "rowN: ", rowN, " colN: ", colN      
   
   #old method
   #ax1=fig.add_subplot(max_rows,opt.ncols,plotN, title=titles[p])
   #ax1.set_title(titles[p], size=16)
   #pl.subplots_adjust(hspace=0.5)
      
   outer_grid = gridspec.GridSpec(max_rows, int(opt.ncols))   
   
   if (opt.hist == "histX"):
       print "in histX"
       inner_grid = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_grid[p], height_ratios=[1,3], hspace=0)
       ax1 = pl.Subplot(fig, inner_grid[1])  #scatter in second cell
       ax2 = pl.Subplot(fig, inner_grid[0], sharex=ax1)  #hist in first cell
       ax2.set_title(titles[p], size=opt.axis_size)         
   elif (opt.hist == "histY"):
       print "in histY"
       inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer_grid[p], width_ratios=[4,1], wspace=0)
       ax1 = pl.Subplot(fig, inner_grid[0])   #scatter in first cell (col)
       ax2 = pl.Subplot(fig, inner_grid[1], sharey=ax1)   #hist in second cell
       ax1.set_title(titles[p], size=opt.axis_size)
   elif (opt.hist == "histXY"):    
       print "in histXY" 
       inner_grid = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer_grid[p], height_ratios=[1,3], width_ratios=[4,1], wspace=0,hspace=0)
       ax1 = pl.Subplot(fig, inner_grid[2])  #scatter in 3rd cell
       ax2 = pl.Subplot(fig, inner_grid[0], sharex=ax1)  #hist in first cell
       ax2.set_title(titles[p], size=opt.axis_size)
       ax3 = pl.Subplot(fig, inner_grid[3], sharey=ax1)   #hist in 4t cell
       #newax = pl.Subplot(fig, inner_grid[1])
   else:
      ax1 = pl.subplot2grid((max_rows, opt.ncols), (rowN,colN))
      ax1.set_title(titles[p], size=opt.axis_size)
      
   #ax1=pl.subplot(max_rows,opt.ncols,plotN, title=titles[p])

   #iterate over number of sets shown on one plot
   for i in (range(len(plots[p]))):
      data=[]
      #lab = (lst[i].split("_"))[0].split("/")[0]
      if (len(legends) > 0):
         print str(plots[p][i]) +" label: "+str(legends[p][i])
      else:
         print str(plots[p][i])
      data = load_xvg_Data(plots[p][i]) #already replaces @ and & to comments
            
      #take the log of the y axis
      #data[:,1] = np.log(data[:,1])
      
      if (opt.color):
         cmap=cm.rainbow 
         vmax = max(data[:,2])
         vmin = min(data[:,2])        
         if (opt.lmax):
            vmax = opt.lmax
         if (opt.lmin):
            vmin = opt.lmin
#         if ((not opt.onelegend) and (not opt.nolegend)):
         result = ax1.scatter(data[:,0],data[:,1], c=data[:,2], s=float(opt.msize), linewidth=int(opt.linewidth), cmap=cmap, vmin=vmin, vmax=vmax,label=legends[p][i])      
         
#         else:
#            result = ax1.scatter(data[:,0],data[:,1], c=data[:,2], s=float(opt.msize), linewidth=int(opt.linewidth), cmap=cmap, vmin=vmin, vmax=vmax)      
                       
#           newax = fig.add_axes(ax1.get_position())
         if (opt.hist == "histY"):
            newax, kw = matplotlib.colorbar.make_axes_gridspec(ax2,  fraction=0.1, shrink=0.9, orientation="vertical")                        
         elif (opt.hist == "histXY"):
            newax, kw = matplotlib.colorbar.make_axes_gridspec(ax3,  fraction=0.1, shrink=0.9, orientation="vertical")                        
         else: 
            newax, kw = matplotlib.colorbar.make_axes_gridspec(ax1,  fraction=0.1, shrink=0.9, orientation="horizontal")                        

         cbar = pl.colorbar(result,  cax=newax, **kw)
         #if  ((not opt.onelegend) and (not opt.nolegend)):                
         #   ax1.legend(loc='upper left',  scatterpoints=1, fontsize=20)
      else:    
         if (len(legends) > 0):
            #if ( i == 3):
               #opt.msize = 12
            ax1.plot(data[:,0],data[:,1], marker=opt.marker, c=colors[i], markersize=float(opt.msize), linewidth=int(opt.linewidth),label=legends[p][i])
         else:
            ax1.plot(data[:,0],data[:,1], marker=opt.marker, c=colors[i], markersize=float(opt.msize), linewidth=int(opt.linewidth))               
      ax1.grid(True)         
      ax1.tick_params(labelsize=opt.tick_size)      
     
      n_xmin = float(opt.xmin) if opt.xmin else ax1.get_xlim()[0]
      n_xmax = float(opt.xmax) if opt.xmax else ax1.get_xlim()[1]
      n_ymin = float(opt.ymin) if opt.ymin else ax1.get_ylim()[0]
      n_ymax = float(opt.ymax) if opt.ymax else ax1.get_ylim()[1]
      ax1.set_xlim([n_xmin,n_xmax])
      ax1.set_ylim([n_ymin, n_ymax])
  
   xx = np.linspace(n_xmin, n_xmax, 100)
   dataX = data[:,0]
   dataY = data[:,1]

   #linear regression
   if (opt.linearReg): 
      gradient, intercept, lin_r_value, p_value, std_err = stats.linregress(dataX,dataY)
      print >>reg_output_fp, "gradient: %3.3f, intersection: %3.3f, R^2: %3.3f, P-value: %3.3f\n" %(gradient, intercept, lin_r_value**2,p_value)
      ax1.plot(dataX,gradient*dataX+intercept, 'k', linewidth=2, label=r'linear $R^2$=%1.3f' % (lin_r_value**2))
   
   if  ((not opt.onelegend) and (not opt.nolegend)):                
      ax1.legend(loc='upper left',  scatterpoints=1, fontsize=16)
   
   fig.add_subplot(ax1)   
   if (opt.hist is not None):
      if (opt.hist == "histX"):
         orient='vertical'     
         ax2.hist(data[:,0], bins=120, orientation=orient, histtype='stepfilled', facecolor='blue', alpha=0.75, normed=True)       
      elif (opt.hist == "histY"):
          orient='horizontal'
          ax2.hist(data[:,1], bins=40, orientation=orient, histtype='stepfilled', facecolor='blue', alpha=0.75, normed=True)      
      elif (opt.hist == "histXY"):
         ax2.hist(data[:,0], bins=120, orientation="vertical", histtype='stepfilled', facecolor='blue', alpha=0.75, normed=True) 
         ax3.hist(data[:,1], bins=40, orientation="horizontal", histtype='stepfilled', facecolor='blue', alpha=0.75, normed=True) 
         fig.add_subplot(ax3)
         pl.setp(ax3.get_xticklabels(), visible=False)
         pl.setp(ax3.get_yticklabels(), visible=False)
      
      fig.add_subplot(ax2)
      pl.setp(ax2.get_xticklabels(), visible=False)
      pl.setp(ax2.get_yticklabels(), visible=False)
   #pl.xlabel(opt.xlab, fontsize=opt.axis_size)   
   #pl.ylabel(opt.ylab, fontsize=opt.axis_size)   
   if ((not opt.color) and (len(plots[p]) > 1)  and not opt.onelegend and not opt.nolegend):
      ax1.legend(loc='best', numpoints=1, markerscale=2)  
      

#if (not opt.nolegend):
#   ax1.legend(loc='best', numpoints=1, markerscale=2)        
# Set common labels
fig.text(0.5, 0.05, opt.xlab, ha='center', va='center', fontsize=opt.axis_size)
fig.text(0.06, 0.5, opt.ylab, ha='center', va='center', rotation='vertical', fontsize=opt.axis_size)

       
#pl.grid(True)
pl.subplots_adjust(hspace=0.35)
pl.savefig(opt.out)




######## log of some attempted commands which were not used:
#linear regression with linear matrix equation
#A = np.vstack([data[:,0], np.ones(len(data[:,0]))]).T
#   m, c = np.linalg.lstsq(A, data[:,1])[0]
# ax1.plot(data[:,0],m*data[:,0]+c, 'k', linewidth=2)

  
  #xx = np.linspace(n_xmin, n_xmax, 100)
   #quadratic fitting of the data
   #z3 = np.polyfit(data[:,0],data[:,1], 3)
   #p3 = np.poly1d(z3)
   
   #quadratic fit + r_square calculation based on R^2 = 1 - (SSreg/SStotal)
   #SS is the sum of squares (SStot = Sum[(yi-y_avg)**2], SSres = Sum[(yi-func_i)**2] 
#   (z2,res, _, _, _) = np.polyfit(dataX,dataY, 2, full=True)
#   p2 = np.poly1d(z2)	#coefficents, p2(x_val) = y_val
#   y_avg = np.sum(dataY)/len(dataY)
#   SSres = np.sum( (dataY - p2(dataX))**2)
#   SStot = np.sum((dataY-y_avg)**2)
#   quad_r_value = 1 - (SSres/SStot)
#   tr= np.sum((np.polyval(np.polyfit(dataX,dataY, 2), dataX) - dataY)**2) # another way to calc residuals
   #ax1.plot(xx,p2(xx), 'gray', linewidth=3, label=r'quadratic $R^2$=%1.3f' % (quad_r_value**2))


   #exponential fitting
   #yy = func(xx, 2.5, 1.3, 0.5, 0.1)   
   #yn = yy + 0.2*np.random.normal(size=len(xx))
   #popt, pcov = curve_fit(func, data[:,0],yn)
   #popt, pcov = curve_fit(func, data[:,0],[6,-6,0.9, 0.1])
   #print popt

   #ax1.plot(xx,p3(xx), 'k', linewidth=2)
   #ax1.plot(xx,func(xx,*popt), 'k', linewidth=2)
    
#pyutil.scatter_plot_xy("2d_1-2_%s.xvg" %name , "2d_1-2_%s.pdf" %name, "ev1","ev2")
#pyutil.scatter_plot_xy("2d_1-3_%s.xvg" %name,  "2d_1-3_%s.pdf" %name, "ev1", "ev3")
