from pylab import *
from DavidLib import load_xvg_Data
import matplotlib.cm as cm
import matplotlib.colors as clrs
import optparse
import os
import pyutil
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
            titles.extend(parts[1].split())
   for i in range(len(titles)):
      titles[i] = titles[i].replace("_", " ")             
   print titles
   return plots, legends, titles

usage = "Usage: %prog -f base_file -o output [options]\nLabels are taken from the relative path's first word"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-p", dest="param", help="mandatory parameter file for multiple xy scatter data")
parser.add_option("-o", dest="out", help="mandatory Output file for the scatterplot")
parser.add_option("-x", dest="xlab", default="" , help="Label for x-axis. default: \"\"")
parser.add_option("-y", dest="ylab", default="" ,help="Label for y-axis. default: \"\"")
parser.add_option("-l", dest="linewidth", default="3" ,help="Linewidth. default: 3")
parser.add_option("--nr", dest="nrows", metavar="int", help="Use for a strict number of rows")
parser.add_option("--nc", dest="ncols", metavar="int", default="2" ,help="Default number of cols for given data will be 2")
parser.add_option("--axs", dest="axis_size", default="22" ,metavar="int",help="Axis label size. default: 22")
parser.add_option("--tl", dest="tick_size", default="22" ,metavar="int",help="Tick label size. default: 22")
parser.add_option("--ts", dest="title_size", default="28" ,metavar="int",help="Title size. default: 22")
parser.add_option("--fig_w", dest="fig_width", default=24.8 ,metavar="float",help="Figure width. default: 24.8")
parser.add_option("--fig_h", dest="fig_height", default=35.2 ,metavar="float",help="Figure height. default: 35.2")
parser.add_option("--xmin", dest="xmin", help="x axis min. limit")
parser.add_option("--xmax", dest="xmax", help="x axis max. limit")
parser.add_option("--ymin", dest="ymin", help="y axis min. limit")
parser.add_option("--ymax", dest="ymax", help="y axis max. limit")
parser.add_option("--revX", dest="revX", action="store_true", help = "flip X axis (max=1.0 is assumed)")
parser.add_option("--smooth", dest="smooth", action="store_true", help = "plot raw & smooth")
parser.add_option("--nolegend", dest="nolegend", action="store_true", help = "Do not plot the legend")
mandatories = ["param", "out"]
(opt,args) = parser.parse_args()

for m in mandatories:
    if not opt.__dict__[m]:
        print "mandatory option is missing - %s\n"%m
        parser.print_help()
        exit(-1)

fontsz = 22
m_size = 4 
fidpi  = 300
sl = 50
gauss_deg=20      #frames to smooth over - window size
gamma = 4
nalp = 1.0
if (opt.smooth):
   nalp = 0.3
m_lw   = float(opt.linewidth)/2.0

###matplotlib parameters
rcParams['font.sans-serif']='Arial'
rcParams['legend.fontsize']= fontsz-2
rcParams['font.size']= fontsz
rcParams['lines.markersize']= m_size

plots, legends, titles = read_param(opt.param)

cm.rainbow

figure(facecolor='w', figsize=(float(opt.fig_width), float(opt.fig_height)), dpi=fidpi)
print "Using the following files for the plots:"
plotN = 1
max_rows = math.ceil(float(len(plots))/float(opt.ncols))
print max_rows
#iterate over number of plots
for p in (range(len(plots))):
#for lst in plots:
   print "Plot "+str(plotN)
   print plots[p]
   colors=[]      
   norm = clrs.Normalize(0, len(plots[p])-1)
  
   subplots_adjust(hspace=0.5)
   ax1=subplot(max_rows,opt.ncols,plotN, title=titles[p])
   ax1.set_title(titles[p], size=opt.title_size)
   ax1.tick_params(which='both', labelsize=opt.tick_size)
   print legends
   #iterate over number of sets shown on one plot
   for j in (range(len(plots[p]))):
      print plots[p][j]
      #adjust colors
      i=0
      clr = j+(len(plots[p])*i)
      if (len(plots[p])<8):
         colors=['k','g', 'r', 'b','m','c','y']
      else:   
         try:
            colors.append(cm.rainbow(norm(j+(len(plots[p])*i))))   
         except:
            colors[clr] = cm.rainbow(norm(j+(len(plots[p])*i)))
      if (os.path.exists(plots[p][j])):
         data = load_xvg_Data(plots[p][j])
         data[:,0] = data[:,0]
         if (not opt.revX):
            print data[:,0]
            print clr, colors[clr], p, j
            print legends[p][clr]
            plot(data[:,0], data[:,1], color=colors[clr], lw=opt.linewidth,label=legends[p][clr], alpha=nalp)
         else: 
            plot(1-data[:,0], data[:,1], color=colors[clr], lw=opt.linewidth,label=legends[p][clr], alpha=nalp)
         if (opt.smooth):
            plot(data[:,0], pyutil.smoothListGaussian(data[:,1], degree=gauss_deg, gamma=gamma), alpha=1.0,color=colors[clr], lw=m_lw)#, label=opt['legend'][clr])
      else: 
          print "Error: no such file: ", plots[p][j]

             
   xlabel(opt.xlab, fontsize=opt.axis_size)
   ylabel(opt.ylab, fontsize=opt.axis_size)
   if (opt.xmin):
      xlim(xmin=float(opt.xmin))
      if (opt.xmax):         
         xlim(xmax=float(opt.xmax))
      if (opt.ymin):
         ylim(ymin=float(opt.ymin))
      if (opt.ymax):
         ylim(ymax=float(opt.ymax))
      
   plotN = plotN+1
if (not opt.nolegend): 
   legend(loc='best', fancybox='True')
savefig(opt.out)
         







