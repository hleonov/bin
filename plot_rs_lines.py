from pylab import *
from DavidLib import load_xvg_Data, smooth
import matplotlib.cm as cm
import matplotlib.colors as clrs
import pyutil
import os
##################################################################################
#  File: plot_rs_lines.py
#  Written by: Hadas Leonov
#  A script to make a figure of smoothed lines, with the noisy raw data in the 
#  background. Fits well for xvg files of trajectories with data vs. Time. 
#  Takes a parameter file for info re directories, files, title, legend and axes.
#  See example: V20_I60.param
##################################################################################

def read_param(pfile, opt):
   lines = open(pfile).readlines()
   for line in lines:
      if (line.strip()):   
         parts = line.split("=")
         parts[0] = parts[0].strip()
         parts[1] = parts[1].strip()
#         if (('fname' in parts[0]) or ('title' in parts[0]) or ('label' in parts[0])):
         if (('title' in parts[0]) or ('label' in parts[0])):
            opt[parts[0]] = parts[1]
         else:
            opt[parts[0]] = parts[1].split()
   for prm in opt.items():
      print prm

if (len(sys.argv) < 3):
   print "Usage: python plot_rs_lines.py <param_file> <output_file> [no_raw]"
   exit()

param  = sys.argv[1] 
out    = sys.argv[2]
no_raw = int(sys.argv[3])
print no_raw
opt = {}
fontsz = 20
m_lw   = 1.5
m_size = 4 
fisize = (10,7)
fidpi  = 100
sl = 50
gauss_deg=20      #frames to smooth over - window size
gamma = 4
nalp = 0.3

read_param(param, opt)


###matplotlib parameters
rcParams['font.sans-serif']='Arial'
rcParams['legend.fontsize']= fontsz-2
rcParams['font.size']= fontsz
rcParams['lines.markersize']= m_size


#colors=['k','r','g','b','m','c','y']

cm.rainbow

figure(facecolor='w', figsize=fisize, dpi=fidpi)

norm = clrs.normalize(0, len(opt['dirs'])*len(opt['fname'])-1)
colors=[]      

subplot(1,1,1)
suptitle(opt['title'], fontsize=fontsz+4, y=0.94)
for i in (range(len(opt['dirs']))):
   for j in (range(len(opt['fname']))):
      clr = j+(len(opt['fname'])*i)
      if (len(opt['dirs']) < 6 ):
         colors=['k','r','g','b','m','c','y']
      else: 
         try:
            colors.append(cm.rainbow(norm(j+(len(opt['fname'])*i))))   
         except:
            colors[clr] = cm.rainbow(norm(j+(len(opt['fname'])*i)))
   
      if (os.path.exists(opt['dirs'][i]+"/"+opt['fname'][j])):
         data = load_xvg_Data(opt['dirs'][i]+"/"+opt['fname'][j])
         avg = data.mean(1)
         if (no_raw == 0):
            plot(0.001*data[:,0], 10*data[:,1], color=colors[clr], lw=m_lw, alpha=nalp)
      else:
         print "Error: no such file: ",opt['dirs'][i]+"/"+opt['fname'][j]

for i in (range(len(opt['dirs']))):
   for j in (range(len(opt['fname']))):
      clr = j+(len(opt['fname'])*i)
      if (os.path.exists(opt['dirs'][i]+"/"+opt['fname'][j])):
         data = load_xvg_Data(opt['dirs'][i]+"/"+opt['fname'][j]) #already replaces @ and & to comments
         #smooth by cookbook (normal, hanning, blackman.. )
         #plot(0.001*data[:,0], pyutil.smooth(data[:,1]*10, sl), color=colors[clr], lw=m_lw, label=opt['legend'][clr])
         #Gaussian-smoothing
         plot(0.001*data[:,0], pyutil.smoothListGaussian(data[:,1]*10, degree=gauss_deg, gamma=gamma), color=colors[clr], lw=m_lw, label=opt['legend'][clr])

xlabel(opt['xlabel'])
ylabel(opt['ylabel'])
#ylim(0,8)


##!!! use (1.0, 1) to put it outside the box, upper right
lgd = legend(bbox_to_anchor=(1.0,1), loc='upper left', borderaxespad=0.)
#lgd = legend(bbox_to_anchor=(0.9,0.1), loc='best', borderaxespad=0.)
#show()
savefig(out, bbox_extra_artists=(lgd,), bbox_inches='tight', pad_inches=0.2)

## to put legend inside the box
#legend(loc='best', fancybox='True')
#savefig(out)
