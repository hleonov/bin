#!/usr/bin/python

#get residue ids of a contact from one line
def get_resids(line):
	parts = str.split(line)
	return (int(parts[2]), int(parts[6]))

#input: 2 tuples of resid1-2, and compare numerically.
#compare resid1 first, then resid2
#returns: 1 if c1 is smaller (comes earlier in contact output sequence) than c2
#	  0 if they are equal
#	  -1 if c1 is larger
def compare_contacts(c1, c2):
	if (c1[0] < c2[0]):
		return 1
	elif (c1[0] > c2[0]):
		return -1
	else:
		if (c1[1] < c2[1]):
			return 1
		elif (c1[1] > c2[1]):
			return -1
		else:
			return 0
			
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--f1", dest="infile1", help="contacts.dat input file 1")
parser.add_option("--f2", dest="infile2", help="contacts.dat input file 2")
parser.add_option("-o", "--output", dest="output", help="output", default="combined_contacts.dat")

(options, args) = parser.parse_args()

if not (options.infile1 and  options.infile2):
	print "please specify the two contact.dat input files"
	print "-c"
	exit(1)
try:
	infile1=open(options.infile1, "r")

except:
	print "Error opening --f1 file " + options.infile1
	exit(1)
try:
	infile2=open(options.infile2, "r")	
except:
	print "Error opening --f2 file " + options.infile2
	exit(1)

outfile=open(options.output, "w")	
contacts=[]
line1=go1=1	#go=proceed to next line
line2=go2=1	#line==true if there are still more lines
c1_res1 = c1_res2 = 0
c2_res1 = c2_res2 = 0


while (line1 or line2):
	#several cases: 
	#both line1 and line2 are read anew and parsed now
	#one of the files has ended (reuse comp and set to 2)
	#both lines ended : break
	#cases for pointers go1-2: 
	# equal lines - both printed and we move on
	# one line is kept until it's combined with another or written out on its turn

	#both lines reached the end of line
	if (not (line1 or line2)):
		break
 	#otherwise, if pointers point to proceed, read next line
 	if (line1 and go1):
		line1 = infile1.readline()	
	if (line2 and go2):
		line2 = infile2.readline()
			
	if ((line1) and (list(line1)[0] not in [';', '#', '@']) ):
		#print line1[0:40],
		(c1_res1, c1_res2) = get_resids(line1)
	if ((line2) and (list(line2)[0] not in [';', '#', '@']) ):
		#print "\t",line2[0:40],
		(c2_res1, c2_res2) = get_resids(line2)

	comp = compare_contacts((c1_res1, c1_res2),(c2_res1, c2_res2))
	#check if one of the files has ended 
	if ((line2 and not line1) or (line1 and not line2)):
		comp=2
	#print " comp=",comp,
	#same contact, append frames		
	if (comp == 0):
		if (line1.strip() and line2.strip()): 		
		#new_row = line1[0:40]+" "+line2[0:40]+"\n"
			list1=str.split(line1)
			list2=str.split(line2)
			new_row = line1[0:50] + str.join(" ",list1[9::]) +" "+ str.join(" ",list2[9::])+"\n"
		#new_row = list1[9::]+list2[9::]
			outfile.write(new_row)
		go1=go2=1
	#first contact comes first, nothing to join. but keep line2 pointer 
	elif (comp == 1):
		if line1.strip():
			outfile.write(line1)
		go1=1
		go2=0		
	#second contact comes first, keep line1 pointer	
	elif (comp == -1):
		if line2.strip(): 
			outfile.write(line2)
		go1=0
		go2=1
	#one file ended
	elif (comp == 2): 
		if (not line1):
			if line2.strip():
				outfile.write(line2)	
			go1=0
			go2=1
		elif (not line2):
			if line1.strip():
				outfile.write(line1)
			go1=1
			go2=0
	#print " go1=",go1," go2=",go2
infile1.close()
infile2.close()
outfile.close()
