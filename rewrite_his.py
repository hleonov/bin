import sys, os

def finalize_his(hd, he, his, new):
   if he and hd:
       rep = "HIP"
   elif he:
       rep = "HIE"
   else:
       rep = "HID"
   his = his.replace("HIS",rep)
   new.append(his)
   hd = 0
   he = 0
   return hd, he, his, new

################## Main ########################

infile  = sys.argv[1]
outfile = sys.argv[2]
all_lines = open(infile).readlines()
new = []
his = ""
hd = 0
he = 0
still_his = 0

for line in all_lines:
    #any line but HIS is automatically written
    if ((line[:4] != "ATOM") or ("HIS" not in line)):
        if still_his:
            hd,he,his,new = finalize_his(hd, he, his, new)
            his = ""
            still_his = 0
        new.append(line)
    #HIS case, read entire residue and check for hydrogens
    else:
        still_his = 1
        his = his+line
        if "HE2" in line[13:16]:
            he = 1           
        if "HD2" in line[13:16]:
            hd = 1
            
#for HIS which are the last residue and line of the file
hd, he, his, new = finalize_his(hd, he, his, new)

out = open(outfile,'w')
for el in new:
   print >>out, el,
         
