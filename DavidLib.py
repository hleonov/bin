from pylab import *
from cStringIO import StringIO
from subprocess import Popen
import re, os


### Classes 
class Sim():
  from subprocess import Popen
  def __init__(self, path, tpr, ndx, xtc=None, edr=None, pdb=None):
    self.path = path
    
    if xtc == None:
      self.xtc = os.path.splitext(tpr)[0]+ '.part0001.xtc'
    else:
      self.xtc = xtc
      
    if edr == None:
      self.edr = os.path.splitext(tpr)[0]+ '.part0001.edr'
    else:
      self.edr = edr
      
    if pdb != None:
      self.pdb = pdb
      self.abs_pdb = os.path.join(self.path, self.pdb)
      
    self.tpr  = tpr
    self.ndx  = ndx
    
    #self.abs = []
    self.abs_tpr = os.path.join(self.path, self.tpr)
    self.abs_ndx = os.path.join(self.path, self.ndx)
    self.abs_xtc = os.path.join(self.path, self.xtc)
    self.abs_edr = os.path.join(self.path, self.edr)
    
    
    #self.ndx_groups = self.ndxes(ndx)
    
  def __repr__(self):
    return self.path + os.path.splitext(self.tpr)
    
  def ndxes(self, ndxfiles):
    ndx_groups = []
    for n in ndxfiles:
      for i, tag in enumerate(re.findall( '\[ (.+) \]', open(n).read() )):
	ndx_groups.append(n, i, tag)
    return ndx_groups


def PlotEV_Structure(X,Y):
   Z, Zx, Zy = histogram2d(X,Y, bins=len(X)/2000.0)
   Zx -= (Zx[1]-Zx[0])*0.5; Zy -= (Zy[1]-Zy[0])*0.5
   contour(Zx[1:],Zy[1:],Z.T,  [1,10,100], linewidths=[2,2,2], linestyles=['dotted','dashdot','solid'], colors='k') 
  
### Systemcalls
def run(command, e=False, o=False):
    if e or o:
	ro,re = Popen(command, shell=True, stdout=-1, stderr=-1).communicate()
	if not e:
	    return ro
	elif not o:
	    return re
	else:
	    return re, ro
    else:
	Popen(command, shell=True, stdout=0, stderr=0).communicate()

### Statistik
def Stat_error(data,err_step=4):
    delta = int(len(data)/err_step)
    means = []
    for i in range(err_step):
	means.append(average(data[delta*i: delta*(i+1)]))
    return std(array(means))#/sqrt(float(err_step)-1.0)
    

### Plotting
def draw_Box(xa,xe,ya,ye, ec='k', fc='none', rundness= 0.0):
  from matplotlib.patches import FancyBboxPatch
  p_bbox = FancyBboxPatch((xa, ya),
                           xe-xa, ye-ya,
                           boxstyle="square,pad=%f"%rundness,
                           ec="k", fc="none", zorder=10.
                           )
  gca().add_patch(p_bbox)
   



### File IO
def load_xvg_Data(File):
  s = open(File).read().replace('@','#').replace('&','#')
  return loadtxt(StringIO(s))

def xvg_Vectors(File): #not tested !!!!
  Dat = load_xvg_Data(File).T
  Time = Dat[0]; Dat = Dat[1:].reshape(len(Dat,-1,3)) #<-- not tested!
  

### Vector operationen
def unitVect(a):
  if len(shape(a)) == 1:
    return a / norm(a)
  if len(shape(a)) == 2:
    return a / outer( sqrt((a**2).sum(1)), ones(shape(a)[1]) )

def dist(a,b):
  return sqrt(((a-b)**2).sum())
  
def angle(a,b):
    return arccos(  sum( unitVect(a)*unitVect(b) , len(shape(a))-1 )  )*180.0/pi #scalar product

def RotMatrix(Vref,Vrot):
  "returns a rotational Matrix nedded to rotate Vrot onto Vref ==> Vref = RotMatrix * Vrot"
  vref =  unitVect(Vref)      # Normalizing Vectors
  vrot =  unitVect(Vrot)
  cp 	= sum(vref*vrot)      #cos phi calculating the angle (directly as the cosine of the angle)
  phi	= arccos(cp)	      # phi between vref and vrot
  sp	= sin(phi)	          # sin phi
  ra 	= unitVect(cross(vrot, vref)) #rotation axis normalized
  D1 = outer(ra, ra)*(1-cp)
  rsp = ra* sp  
  D2 = 	 array([[ cp,		-rsp[2], 	 rsp[1] ],
	        [ rsp[2],	 cp,		-rsp[0] ],
	        [-rsp[1],	 rsp[0],	 cp    ]])
	        
  return D1+D2


### Signal operations
## --> stolen from the web!
def smooth(x, window_len=11, window='hanning'):
     if x.ndim != 1:
         raise ValueError, "smooth only accepts 1 dimension arrays." 
     if x.size < window_len:
         raise ValueError, "Input vector needs to be bigger than window size." 
     if window_len<3:
         return x 
     if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
         raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'" 
     s = concatenate([ 2*x[0]-x[window_len:1:-1] , x, 2*x[-1]-x[-1:-window_len:-1] ])
     if window == 'flat': #moving average
         w = ones(window_len,'d')
     else:
         w = eval(window+'(window_len)') 
     y=convolve(w/w.sum(),s,mode='same')
     return y[window_len-1:-window_len+1]

## Ableitungen
def diff_1(x):
  return x[1:]-x[:-1]

## Autocorrelation
def Autocorrelation(X):
  x = X - mean(X)
  c = correlate(x,x,'full')
  return c[len(c)/2:] /max(c)


### Fitting script
def Gaussian_KBT_fit(Para_Vector, x):
  "Gaussfit mit p[0]=A=Skalirungsfaktor, p[1]=k=Kraftkonstante, p[2]=x0"
  T = 300; kB = 1.3806504e-23; p = Para_Vector
  return p[0]*exp(-0.5/T/kB*p[1]*(x-p[2])**2)
 
 
def fit(function, init_para_Vector, messdaten_x, messdaten_y):
    "Function must look like f([p0,p1,p2...], messdaten_x)"
    from scipy import optimize  # load bib 
    def f_resi(para_Vector):    # define function to calculate the residue
        return messdaten_y - function(para_Vector, messdaten_x)
    return optimize.leastsq(f_resi, init_para_Vector) #return the fit!



### PDB operationen
def read_PDB_Coords(Filename):
  s = open(Filename).read().split("ENDMDL")[0]
  return array(re.findall("ATOM.{28}([-\.0-9]+)  ([-\.0-9]+)  ([-\.0-9]+)",s), dtype = 'float32')
  

def pdb_Vector(start, direction, spacer=0.3, scale = [0,1], Atomtype='V', Residuetype='VEC', Chain=' ', resid=1):
  pdb_format="ATOM  %5d %-4s %3s%2s%4d %11.3f %7.3f %7.3f %5.2f %5.2f\n"
  s = ''
  length = float(  sqrt( (direction**2).sum() )  )
  for i,v in  enumerate(arange(scale[0],scale[1]+.001,spacer/length)[:,newaxis]*direction + start):
    s += pdb_format%(i, Atomtype, Residuetype, Chain, resid, v[0], v[1], v[2], 1.0, 0.0)
  return s

def pdb_Atoms(Atoms, Atomtype='A', Residuetype='HEL', Chain=' ', resid=1):
  pdb_format="ATOM  %5d %-4s %3s%2s%4d %11.3f %7.3f %7.3f %5.2f %5.2f\n"
  s = ''
  for i,a in  enumerate(Atoms):
    s += pdb_format%(i, Atomtype, Residuetype, Chain, resid, a[0], a[1], a[2], 1.0, 0.0)
  return s
