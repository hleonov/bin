from __future__ import division 
import sys, traceback
from glob import glob
import subprocess
import numpy
import math
import pylab as pl
#from DavidLib import load_xvg_Data, smooth
from cStringIO import StringIO
import re, os

### File IO
def load_xvg_Data(File):
  s = open(File).read().replace('@','#').replace('&','#')
  return pl.loadtxt(StringIO(s))
  
  
## System commands    
def run_command(command,exe='/bin/bash') :
    print command
    p = subprocess.Popen(command, shell=True, executable=exe ,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    #print out,err
    if (err != ""):
        print err
    if (p.wait() != 0):
        print "Failed command: "+command        
        retval = -1
        #exit()
    else:
        retval = 0 

    return out, retval 

#concat ATOM records, not HETATMS    
def concat_pdbs(pdbs, out) :
    print "concat "+str(pdbs)
    print "into "+out
    com = []
    #first pdb - get header, atoms, ter (all but "end")
    com.append("grep -v END %s > %s" % (pdbs[0], out))
    #get rest
    for i in range(1,len(pdbs)):
        com.append("grep ATOM %s >> %s" % (pdbs[i], out))
        com.append("echo TER >> %s" % out)
    com.append("echo END >> %s" % out)
    for k in com:
        run_command(k)


## Math stuff
# normalize values of a list to make its max = normalizeTo
def normList(L, normalizeTo=1):    
    vMax = max(L)
    return [ x/(vMax*1.0)*normalizeTo for x in L]
    
def dist(a,b):
  return sqrt(((a-b)**2).sum())


# Compute probabilities and free energies when drawing an arbitrary cross on a scatter plot
# The cross is defined by T1 (threshold on the X axis), and T2 (y-axis). 
# ref is the reference tile, for the free energy calculation
# The numbering of the tiles:
#  1 \ 2
# -------
#  3 \ 4

def scatter_probabilities(infile, T1, T2, ref=4, Tmp=300):
   k_B = 1.3806504000000001e-23
   ref=int(ref)
   print "Using x threshold: %s\n      y threshold %s\n      temperature %s\n      ref tile: %s" %(T1, T2, Tmp, ref)
   print "Tiles numbering:\n 1 | 2\n ------\n 3 | 4"
   l = open(infile).readlines()
   total=len(l)
   tile_count=[0,0,0,0,0]
   frac = [0,0,0,0,0]
   #free_energy=[0,0,0,0,0]
   for line in l:     
      tiles=[1,2,3,4]
      parts=line.split()
      if (float(parts[0]) < T1):
         tiles=list(set(tiles) & set([1,3]))         
      else:
         tiles=list(set(tiles) & set([2,4]))
      if (float(parts[1]) < T2):
         tiles=list(set(tiles) & set([3,4]))
      else:
         tiles=list(set(tiles) & set([1,2]))
      tile_count[tiles[0]] = tile_count[tiles[0]] +1 
   print "Tile  probability       fraction (m/1-m)"
   for i in range(1,5):
      tile_count[i] = tile_count[i]/len(l) 
     # frac = tile_count[i] / tile_count[ref]
      if (numpy.abs(1-tile_count[i])>0.000001):
         frac[i] = tile_count[i] / (1-tile_count[i])
      else: 
         frac = 1
      #free_energy[i] = -1*k_B*Tmp*numpy.log(frac)      
      print i, " \t\t", tile_count[i], "\t",frac[i]  
   return frac[ref]   
      #print tiles, line,
        

### Signal smoothing
## --> from SciPy.org cookbook
def smooth(x, window_len=11, window='hanning'):
     if x.ndim != 1:
         raise ValueError, "smooth only accepts 1 dimension arrays." 
     if x.size < window_len:
         raise ValueError, "Input vector needs to be bigger than window size." 
     if window_len<3:
         return x 
     if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
         raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'" 
     s = pl.concatenate([ 2*x[0]-x[window_len:1:-1] , x, 2*x[-1]-x[-1:-window_len:-1] ])
     if window == 'flat': #moving average
         w = ones(window_len,'d')
     else:
         w = eval('numpy.'+window+'(window_len)') 
     y=numpy.convolve(w/w.sum(),s,mode='same')
     return y[window_len-1:-window_len+1]
     

### Gaussian data smoothing function copied from the web 
### Gaussiam Kernel is 1/exp((x-mu)**2 / 2*sig**2) or can also
### be implemented as  1/exp(gamma*(x-mu)**2). 
### This means that gamma=(1/2*sig**2) or sig=sqrt[1/(2*gamma)]
### gamma=16 was chosen to have edges at a fraction of 0.02
def smoothListGaussian(data,degree=5, gamma=16):  
    list= data.tolist()
    #pad the list
    list= [list[0]]*(degree-1) + list + [list[-1]]*degree
    #compute weights
    window=degree*2-1  
    weight=numpy.array([1.0]*window)  
    weightGauss=[]  
    for i in range(window):  
        i=i-degree+1  
        frac=i/float(window)  
        gauss=1/(numpy.exp( gamma*(frac)**2))  
        weightGauss.append(gauss)  
    weight=numpy.array(weightGauss)*weight  
    print weight
   
    #smooth
    smoothed=[0.0]*(len(list)-window)  
    for i in range(len(smoothed)):  
        smoothed[i]=sum(numpy.array(list[i:i+window])*weight)/sum(weight)  
    return pl.array(smoothed)
## Plots

def scatter_plot_xy(infile,outfile,cfile, xlabel=" ", ylabel=" ", title= " ", writeOut=1, writeLabels=1, markerSize=12):
   colors=['b','r','g','c','m','k']
   print "called with: "+infile+" "+outfile+" "+cfile
   labels = open(cfile).readlines()
   dic={}
   cur_color=0
   for l in labels:
      #l2=l.replace("\n","")
      dic[l] = None
   for l in dic:
      dic[l] = cur_color
      cur_color += 1   
            
   print dic, len(dic)
   print colors, len(colors)
   try:
      stored = []
      incr = []
      r=0.8
      ang=30
      print "label radius: %2.3f, incremental angle: %d" %(r, ang)
      data=numpy.loadtxt(infile,comments='#') 
      print "#labels: %d, #data points: %d" %(len(labels), len(data))
     # pl.plot(data[:,0],data[:,1],'bo')
      for i in range(len(data)):
         #pl.plot(data[i,0],data[i,1],"bo")
         #print colors, dic[labels[i]]
         pl.plot(data[i,0],data[i,1],"%so" % colors[dic[labels[i]]], markersize=markerSize )
         xy_tup = (data[i][0], data[i][1])                  
         xy_tup2 = (data[i][0]+r*numpy.cos(0), data[i][1]+r*numpy.sin(0))

         for j in range(len(stored)):
            if ( numpy.sqrt((stored[j][0]-xy_tup2[0])**2+(stored[j][1]-xy_tup2[1])**2) < r):
               incr[j] = incr[j]+ang
               xy_tup2 = (data[i][0]+r*numpy.cos(incr[j]), data[i][1]+r*numpy.sin(incr[j]))
              # break
         if xy_tup not in stored:
            stored.append(xy_tup2)
            incr.append(0)
         pl.annotate(str(i+1), xy=xy_tup, xycoords='data', xytext=xy_tup2, textcoords='data', fontsize=7) 
        
      x1,x2=pl.xlim()
      y1,y2=pl.ylim()
      (xmin, xmax) = (data[:,0].min()-1, data[:,0].max()+1)
      (ymin, ymax) = (data[:,1].min()-1, data[:,1].max()+1)
      print xmin, xmax
      print ymin, ymax
      pl.xlim(xmin=min(xmin,ymin), xmax=max(xmax, ymax))
      pl.ylim(ymin=min(xmin,ymin), ymax=max(xmax, ymax))

#      pl.xlim(xmin=min(x1,data[:,0].min())-1, xmax=max(x2,data[:,0].max())+1)
#      pl.ylim(ymin=min(y1,data[:,1].min())-1, ymax=max(y2,data[:,1].max())+1)
      print "pyutil data xlim: %s\tylim: %s" %(pl.xlim() , pl.ylim() )
      
      #temporary, add y=x line
      #my_x=numpy.arange(pl.xlim()[0]-1, pl.xlim()[1]+1)
      #pl.plot(my_x,my_x, "c--")
      
      #pl.show()
      pl.xlabel(xlabel)
      pl.ylabel(ylabel)
      if (writeOut):
         pl.suptitle(title)
         pl.savefig(outfile, format='pdf', orientation='landscape')         
         pl.clf()
     

   except:
      print "Cannot make figure"+infile
      traceback.print_exc(file=sys.stdout)
      exit()

def make_scatter_panel(in_list, outfile, xlabel, ylabel,max_c=2):
   for i in range(len(in_list)):
      data=numpy.loadtxt(in_list[i], comments="#")
      pl.subplot(len(in_list)/max_c,max_c,i)
      l=pl.plot(data[:,0],data[:,1],'bo')
      pl.grid(True)
      pl.xlabel(xlabel)
      pl.ylabel(ylabel)
   pl.savefig(outfile, format='eog')
