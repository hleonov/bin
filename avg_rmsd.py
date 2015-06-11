import sys, os
import numpy
import string
import optparse

usage="Usage: %prog [options] <infile>"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-b", dest="b", metavar="int", default=0, help="Start from time (ps) specified")
parser.add_option("-e", dest="e", metavar="int", default=-1, help="End at time (ps) specified")

(opt, args) = parser.parse_args()
#print opt.b, opt.e
#print args

lines = open(args[0]).readlines()

data = []
time = []

for line in lines:
   if line[0] not in ['@','#']:
      entr = line.split()
      if (  float(entr[0])>=float(opt.b) and (opt.e==-1 or float(entr[0]) <= float(opt.e))  ):
         data.append(float(entr[1]))

print 10*numpy.average(data)
