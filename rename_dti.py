import sys, os
from glob import glob
import re
import pyutil

base = sys.argv[1]
new_name = sys.argv[2]

print base, new_name

dirs = os.listdir(".")
for d in dirs:
    match = re.search("%s(\d+\.\d+)" % base, d)
    if match:
        lmb = match.group(1)
        os.chdir(base+lmb)
        pyutil.run_command("mv fe_avg.txt fe_avg_%s.txt" % new_name)
        pyutil.run_command("mv fe_block_1000.txt fe_block_1000_%s.txt" % new_name)
        pyutil.run_command("mv fe_block_2000.txt fe_block_2000_%s.txt" % new_name)
        pyutil.run_command("mv fe_convergence.txt fe_convergence_%s.txt" % new_name)
        os.chdir('..')

files = ["lda_vs_dGdl.txt", "dG_vs_time_cum.txt", "dG_vs_time.txt", "fe_result.txt", "lda_vs_dGdl.png", "random_blocks.png" ]

for f in files:
   if os.path.exists(f):
      com="mv %s %s" %(f, f[:-4]+"."+new_name+f[-4:])
      pyutil.run_command(com)
