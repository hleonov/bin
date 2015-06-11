#!/usr/bin/python
from optparse import OptionParser
parser=OptionParser()
parser.add_option("-f", "--file", dest="infile", help="inputfile")
parser.add_option("-c", "--column", dest="col", default=1,  help="column to sort for maxdist")
parser.add_option("-m", "--max", dest="lim", default=5, help="maximum number of lines listed")
(options, args) = parser.parse_args()

col=int(options.col)
inp=open(options.infile, "r")
array=[]
for line in inp:
	if list(line)[0] != "#" and list(line)[0] != "@" and list(line)[0] != "&":
		array.append(str.split(line))

array.sort(key=lambda x: float(x[col]))

for i in range(int(options.lim)):
	print( " ".join(array[len(array)-i-1]) )
