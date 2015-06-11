import pyutil, sys

if (len(sys.argv)<4):
   print "Usage: python smooth_xvg.py <infile> <outfile> <win-size>"
   exit()

infile = sys.argv[1]
outfile = sys.argv[2]
winsize = sys.argv[3]
gamma    = sys.argv[4]

data=pyutil.load_xvg_Data(infile)
new_data=data
new_data[:,1] = pyutil.smoothListGaussian(data[:,1], degree=int(winsize), gamma=float(gamma))

fp = open(outfile, 'w')

for i in range(len(new_data)):
   print >>fp, new_data[i,0], new_data[i,1]
fp.close()

