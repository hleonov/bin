import sys, os
import numpy
import string
import optparse
from glob import glob
import pyutil

usage = "Usage: %prog [options]"
parser = optparse.OptionParser(usage=usage)
parser.addoption("-h", "--help", action="help")
parser.addoption("-d", "--base_dir", dest="base_dir", help="Base directory to look for files")
parser.addoption("-f", "--base_file",dest="base_file", help="Base filename to look for")
parser.addoption("-c", "--col", dest="col",default=1, help="coloum number to start averaging")
(opt, args) = parser.parse_args()
print args, opt

files = glob("%s*/*%s*" %(opt.base_dir, opt.base_file))

print "files to analyze: ", files

data1 = pyutil.load_xvg_Data(files[0])
data2 = pyutil.load_xvg_Data(files[0])
3d_array = nu
for i in (range(1,len(files))):
   data = pyutil.load_xvg_Data(files[i])
   3d_array = numpy.dstack(
