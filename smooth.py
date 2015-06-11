#!/usr/bin/python

from string import split
from math import *

from sys import argv

def average(l):
    return sum(l)*1./len(l)

def smooth(fin, fout, param=150):
	f = open(fin).readlines()
	while f[0][0] in ['@','#']:
		f = f[1:]
	f = [map(float, split(x)) for x in f]
	smooth = []
	for i in range(len(f)-param):
		smooth.append([f[i][0], average([x[1] for x in f[i:i+param]])])
	fout = open(fout,"w")
	for i in smooth:
		fout.write(str(i[0]) + "   " + str(i[1]) + "\n")
	fout.close()
	print "OK"
	
smooth(argv[1], argv[2], int(argv[3]))
