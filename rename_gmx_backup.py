from glob import glob
import pyutil as pu
import optparse
import os, sys

usage = "Usage: %prog [options]"
parser.add_option("-p", dest="prefix", help="prefix")
parser.add_option("-s", dest="suffix", help="suffix")
(opt, args) = parser.parse_args()
print opt, args
print "Searching for files in the form of %s*%s" %(opt.prefix, opt.suffix)
files=glob("%s*%s" %(opt.prefix, opt.suffix))

for f in files:
   pu.run_command("mv \%s\# %s" % (f[0:-1], f[1:-3]))
