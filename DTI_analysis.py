import sys, os
from glob import glob
from numpy import *
import re 
import optparse
import pyutil

def read_data(fn,b=0,e=-1):
    print "b="+str(b)+" e="+str(e)
    l = open(fn).readlines()
    data = []
    for line in l:
        if line[0] not in ['@','#']:
            entr = line.split()            
            if (  float(entr[0])>=float(b) and (e==-1 or float(entr[0]) <= float(e)) and (len(entr)==2)  ):
               # print line
                data.append( float(entr[1] ) )
    print >>sys.stderr, 'Read file:', fn, ' with %d data points' % len(data)
    return data

#given in ps
#if I write every 500 steps (0.5ps) then there are 50000 points in 1000ps
def datapoint_from_time(time):
    return time*50

def block_aver( data, block_size = 1000, offset = 0):
    # make blocks of 1 ns length
    total_time = len(data) / 50.
    print "end time "+str(total_time)
    next_time = block_size
    results = []
    while next_time < total_time:
	     #the offset is only for first window
        beg = datapoint_from_time(offset)
        end = datapoint_from_time(next_time)
        res = average( data[beg:end] )
        #results.append( (offset+block_size*.5, res ) )
        results.append( (next_time, res ) )
        offset = next_time
        next_time += block_size
    return results

def convergence( data, block_size = 1000, offset = 0):
    # make blocks of 1 ns length
    total_time = len(data) / 50.
    next_time = block_size
    results = []
    while next_time < total_time:
        beg = datapoint_from_time(offset)
        end = datapoint_from_time(next_time)
        res = average( data[beg:end] )
        results.append( (next_time, res ) )
        next_time += block_size
    return results

def concat(base) :
   files = glob("%s.part*.xvg" %base)
 #  if (not os.path.exists(base+".xvg")):
   pyutil.run_command("cat %s.part0001.xvg | awk '{if (NF>=2) print $0}' > %s.xvg" %(base, base))
   for f in files:
      if (f != "%s.part0001.xvg" %base):
         print "not first: "+f
         pyutil.run_command("grep -v \"^[@#]\" %s |  awk '{if (NF==2) print $0}' >> %s.xvg" %(f,base))
   if (len(files) > 1 and (not os.path.exists(base+".xtc"))):
      pyutil.run_command("trjcat -f %s.part*.xtc -o %s.xtc" %(base, base))
      pyutil.run_command("trjcat -f %s.part*.trr -o %s.trr" %(base, base))
      pyutil.run_command("eneconv -f %s.part*.edr -o %s.edr" %(base, base))
#========================= MAIN =============================

usage="Usage: %prog [options] <dir_basename>"
parser = optparse.OptionParser(usage=usage)
#parser.add_option("-h","--help", action="help")
parser.add_option("-b", dest="b", metavar="int", default=0, help="Start from time (ps) specified")
parser.add_option("-e", dest="e", metavar="int", default=-1, help="End at time (ps) specified")
parser.add_option("-i", dest="ti_name", default="ti_", help="base name for TI file (default: \"ti_\")")
(opt, args) = parser.parse_args()
print opt.b, opt.e, opt.ti_name
print args

dirs = os.listdir(".")
print opt.ti_name

#an easier way to traverse directories:
# dirs = glob("%s*.*" % args[0])
for d in dirs:
    match = re.search("%s(\d+\.\d+)" % args[0], d)
    if match:
        lmb = match.group(1)
        print "Processing -> directory %s" % (args[0] + lmb) 
        os.chdir(args[0]+lmb)
        concat(opt.ti_name+lmb)    
        data = read_data(opt.ti_name+lmb+'.xvg',b=opt.b,e=opt.e)        
        #data = read_data(ti_name+lmb+'.part0001.xvg',b=options.b,e=options.e)        
        #data = read_data("dhdl.part0001.xvg") 
        fp = open('fe_avg.txt','w')
        print >>fp, average(data), std(data),  len(data)
        fp.close()
        res = block_aver( data, 1000 )
        fp = open('fe_block_1000.txt','w')
        for t, r in res:
            print >>fp,  t, r
        fp.close()
        res = block_aver( data, 2000 )
        fp = open('fe_block_2000.txt','w')
        for t, r in res:
            print >>fp, t, r
        fp.close()

        res = convergence( data, 1000 )
        fp = open('fe_convergence.txt','w')
        for t, r in res:
           print >>fp, t, r
        fp.close()
        print '.............done'
        os.chdir('..')
        
