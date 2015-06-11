import pylab as pl
import numpy as np
import pyutil 
from glob import glob
import re

rows = cols = 20
frows = 4 
fcols = 3
v_min = 0

#frame_comp_dist_pls_filter_ew_129_463.dat
files = glob("frame_comp_dist_pls_filter_[0-9]*.dat")
figure = pl.figure(figsize = (20,16))

for k,f in enumerate(files):
   data = pyutil.load_xvg_Data(f)
   contact = re.search("(\d+\_\d+)", f).group(1)
   arr = np.zeros((rows, cols))
   for i in range(len(data)):
      arr[int(data[i,1])-1, int(data[i,0])-1] = data[i,2] 
   #pl.subplot(frows,fcols,k)
   print k, f
   pl.subplot(frows,fcols,k+1)
   pl.title(contact, fontsize=20)
   pl.imshow(arr, interpolation="none")
   pl.colorbar()
pl.savefig("mode_dist.png", dpi=300) 


