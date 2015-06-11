# Initialize
import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
import chimera
import getopt
import re
import math

### CLOSE ALL Before start
rc('close all')

# import objects
#pdb
aqp0 = chimera.openModels.open('/home/rbrione/cazuelas/AQP0/4_Camilo/CCC/ccc_chimera_autmat/2b60.pdb',type="PDB")
rc('background solid white')
rc('set projection orthographic')
rc('turn x 90')
rc('~disp :#0')
rc('represent stick :MC3')
rc('color dark gray :MC3')
rc('color byhet :MC3')
rc('disp :MC3')


#
# open densities
# experimental dmpc
#dmpc_exp_vop = chimera.openModels.open('/home/rbrione/cazuelas/AQP0/4_Camilo/CCC/ccc_chimera_autmat/trim02/volume_diff_masked_vop_trim.mrc')
dmpc_exp_vop3A = chimera.openModels.open('/home/rbrione/cazuelas/AQP0/4_Camilo/CCC/ccc_chimera_autmat/trim03/volume_diff_masked_vop_trim_zone3A_retrim.mrc')
rc('volume #1 step 1')
rc('volume #1 showOutlineBox true')
rc('volume #1 outlineBoxRgb black')
### Set countour level 0
##    rc('volume #1 level %f' %  (l/2.))
##    print 'CUTOFF EXP. MAP %0.1f' % (l/2.)
##levt = 0.5 # exp density cutoff
##rc('volume #1 level %f' % levt)

# opls dmpc
dmpc_opls = chimera.openModels.open('/home/rbrione/cazuelas/AQP0/4_Camilo/CCC/ccc_chimera_autmat/trim03/dmpc_monomer_opls_trim.mrc')
rc('volume #2 step 1')
rc('volume #2 showOutlineBox true')
rc('volume #2 outlineBoxRgb red')
##    print 'CUTOFF EXP. MAP %0.1f' % (l/2.)
##levo = 0.004 # 1 sigma = 0.001
##rc('volume #2 level %f' % levo)

# c36 dmpc
dmpc_c36 = chimera.openModels.open('/home/rbrione/cazuelas/AQP0/4_Camilo/CCC/ccc_chimera_autmat/trim03/dmpc_monomer_c36_trim.mrc')
rc('volume #3 step 1')
rc('volume #3 showOutlineBox true')
rc('volume #3 outlineBoxRgb yellow')
###levc = 0.004 # 1 sigma = 0.00095
##rc('volume #3 level 0.004')

# c36 dmpc tip3p 01
dmpc_c36_tip3p1 = chimera.openModels.open('/home/rbrione/cazuelas/AQP0/4_Camilo/CCC/ccc_chimera_autmat/trim03/dmpc_monomer_c36_tip3p_01_trim.mrc')
rc('volume #4 step 1')
rc('volume #4 showOutlineBox true')
rc('volume #4 outlineBoxRgb green')
###levc = 0.004 # 1 sigma = 0.0011
##rc('volume #4 level 0.004')

# c36 dmpg
dmpg_c36 = chimera.openModels.open('/home/rbrione/cazuelas/AQP0/4_Camilo/CCC/ccc_chimera_autmat/trim03/dmpg_monomer_c36_trim.mrc')
rc('volume #5 step 1')
rc('volume #5 showOutlineBox true')
rc('volume #5 outlineBoxRgb gray')
###levc = 0.006 # 1 sigma = 0.0015
##rc('volume #5 level 0.006')

# c36 popc
popc_c36 = chimera.openModels.open('/home/rbrione/cazuelas/AQP0/4_Camilo/CCC/ccc_chimera_autmat/trim03/popc_monomer_c36_trim.mrc')
rc('volume #6 step 1')
rc('volume #6 showOutlineBox true')
rc('volume #6 outlineBoxRgb cyan')
###levc = 0.005 # 1 sigma = 0.00124
##rc('volume #6 level 0.005')

# c36 sopc
sopc_c36 = chimera.openModels.open('/home/rbrione/cazuelas/AQP0/4_Camilo/CCC/ccc_chimera_autmat/trim03/sopc_monomer_c36_trim.mrc')
rc('volume #7 step 1')
rc('volume #7 showOutlineBox true')
rc('volume #7 outlineBoxRgb magenta')
###levc = 0.004 # 1 sigma = 0.00089
##rc('volume #7 level 0.004')

# c36 sope
sope_c36 = chimera.openModels.open('/home/rbrione/cazuelas/AQP0/4_Camilo/CCC/ccc_chimera_autmat/trim03/sope_monomer_c36_trim.mrc')
rc('volume #8 step 1')
rc('volume #8 showOutlineBox true')
rc('volume #8 outlineBoxRgb blue')
#levc = 0.004
### 1 sigma = 0.000925
##rc('volume #8 level 0.004')



# Cross Correlation Coefficient calculation
# measure correlation  map-model1  map-model2  [ aboveThreshold true|false ] [ rotationAxis axis ] [ angleRange start,end,step ] [ plot true|false ]
#
# Before calculating the CCC the option 'volume level threshold-level' defines above which values to use in the first map.
# In the trimmed vop map values of rho from -4 to 4. Could try a range from 0.0 to 3.0. Where there is significant density in the experimental map.
# rc('volume #1 level 1.0')

# Boxes Zx dimension
nsl = 10  #number of slices
zt = 119. #z bins exp  masked trimmed
zo = 131. #z bins opls 
zc = 105. #z bins c36  
###

print 'experimental slice level ALL'
print 'opls slice level ALL'
print 'c36 slice level ALL sigma'
print '%i number of slices'

for zet in range(0, nsl):
      print 'experimental slice from z=%i to %i ' % ((round(zet*(zt/nsl))), (round((zet+1)*(zt/nsl))))
      print 'opls slice from z=%i , %i' % (round((zo/nsl)*zet), round((zo/nsl)*(zet+1)))
      print 'c36 slice from z=%i , %i' % (round((zc/nsl)*zet), round((zc/nsl)*(zet+1)))
      rc('volume #1 region 0,0,%d,102,102,%d' % ((zt/nsl)*zet, (zt/nsl)*(zet+1)))
      rc('volume #2 region 0,0,%d,116,116,%d' % ((zo/nsl)*zet, (zo/nsl)*(zet+1)))
      rc('volume #3 region 0,0,%d,94,93,%d' % ((zc/nsl)*zet, (zc/nsl)*(zet+1)))
      rc('volume #4 region 0,0,%d,94,93,%d' % ((zc/nsl)*zet, (zc/nsl)*(zet+1)))
      rc('volume #5 region 0,0,%d,94,93,%d' % ((zc/nsl)*zet, (zc/nsl)*(zet+1)))
      rc('volume #6 region 0,0,%d,94,93,%d' % ((zc/nsl)*zet, (zc/nsl)*(zet+1)))
      rc('volume #7 region 0,0,%d,94,93,%d' % ((zc/nsl)*zet, (zc/nsl)*(zet+1)))
      rc('volume #8 region 0,0,%d,94,93,%d' % ((zc/nsl)*zet, (zc/nsl)*(zet+1)))
#      rc('copy file lero%i.png' % zet)
      
      for i in range(1,9):
##    #   rc('measure mapStats #%d' % i) # Map statistics
          for j in range(1,9):
##    #       rc('measure mapStats #%d' % j)
              print '%d,%d' % (i, j)
              rc('measure correlation #%d #%d aboveThreshold false plot true' % (i, j))
##        
##            
##     rc('measure correlation #2 #1 aboveThreshold true plot true')  
         
