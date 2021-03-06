from pylab import *
from DavidLib import load_xvg_Data, smooth

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
         if (('fname' in parts[0]) or ('title' in parts[0]) or ('label' in parts[0])):
            opt[parts[0]] = parts[1]
         else:
            opt[parts[0]] = parts[1].split()
   for prm in opt.items():
      print prm

param  = sys.argv[1] 
out   = sys.argv[2]

opt = {}
fontsz = 16
m_lw   = 1.5
m_size = 4 
fisize = (10,7)
fidpi  = 100
sl = 50
nalp = 0.3

read_param(param, opt)


###matplotlib parameters
rcParams['font.sans-serif']='Arial'
rcParams['legend.fontsize']= fontsz-2
rcParams['font.size']= fontsz
rcParams['lines.markersize']= m_size


colors=['k','r','g','b','m','c','y']

cm.rainbow

figure(facecolor='w', figsize=fisize, dpi=fidpi)

subplot(1,1,1)
suptitle(opt['title'], fontsize=fontsz+4, y=0.94)
for i in (range(len(opt['dirs']))):
   data = load_xvg_Data(opt['dirs'][i]+"/"+opt['fname'])
   avg = data.mean(1)
   plot(0.001*data[:,0], 10*data[:,1], color=colors[i], lw=m_lw, alpha=nalp)

for i in (range(len(opt['dirs']))):
   data = load_xvg_Data(opt['dirs'][i]+"/"+opt['fname']) #already replaces @ and & to comments
   plot(0.001*data[:,0], smooth(data[:,1]*10, sl), color=colors[i], lw=m_lw, label=opt['legend'][i])
xlabel(opt['xlabel'])
ylabel(opt['ylabel'])
legend(loc='best', fancybox='True')
#show()
savefig(out)
