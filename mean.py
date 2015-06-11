#!/usr/bin/python
import numpy as np
import math as m
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
		                  help="reading from file", metavar="FILE")
parser.add_option("-c", "--column", dest="column", help="calculate the average and std from which column")
(options, args) = parser.parse_args()
column=int(options.column)
inf=open(options.filename, "r")
array=[]
for line in inf:
	if(list(line)[0]!="#" and list(line)[0]!="@" and list(line)[0]!="\\"):
		array.append(float(str.split(line)[column]))
print("mean: " + str(np.mean(array)))
print("std: " + str(np.std(array)))
