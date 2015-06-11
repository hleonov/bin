import sys, os
import numpy
import string
import optparse
from glob import glob
import pyutil

usage = "Usage: %prog [options]"
parser = optparse.OptionParser(usage=usage)
#parser.add_option("-h", "--help", action="help")
parser.add_option("-d", "--base_dir", dest="base_dir", help="Base directory to look for files")
parser.add_option("-f", "--base_file",dest="base_file", help="Base filename to look for")
parser.add_option("-c", "--col", dest="col",default=1, help="coloum number to start averaging")
parser.add_option("-o", dest="output", help="output file name")
(opt, args) = parser.parse_args()
print args, opt

files = glob("%s*/%s*" %(opt.base_dir, opt.base_file))

print "files to analyze: ", files
data = pyutil.load_xvg_Data(files[0])

print data.shape
prev_3d_array = data[:,:,numpy.newaxis]
#prev_3d_array = numpy.concatenate(tmp, axis=2)
print prev_3d_array.shape
for i in range(1,len(files)):
   print "reading ", files[i]
   data = pyutil.load_xvg_Data(files[i])
   #print data.shape
   tmp = data[ : , : , numpy.newaxis]
   #print tmp.shape
   arr = numpy.concatenate((prev_3d_array,tmp), axis=2)
   prev_3d_array = arr
   print arr.shape

average = arr[:,:,:].mean(axis=2)

numpy.savetxt(opt.output, average, "%8.3f")
#outf = open(opt.output, 'w')
#print >>outf, str(average).replace('[','').replace(']','')
