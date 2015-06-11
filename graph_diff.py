import sys,os
import string

# receives two files with x,y data. X should have the same values.
# substracts the y col for similar x data and writes out
file1 = sys.argv[1]
file2 = sys.argv[2]

def read_file(name):
   l = open(name).readlines()
   dic = {}
   for line in l:
      if line[0] not in ['@','#']:
         entr = line.split()
         dic[float(entr[0])] = float(entr[1])
  # print dic
   return dic
   
data1 = read_file(file1)
data2 = read_file(file2)

for key in (sorted(data1.keys())):
   if key in data2:
      print "%s\t%f" %(key, data1[key]-data2[key])
