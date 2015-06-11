import sys, os
import optparse
import pyutil as pu
from glob import glob
import re
import string

# Define global variables
def init():
   global whatif, data, local_whatif, mybin, amberFF, job, mdp
   whatif        = "/usr/local/whatif/DO_WHATIF.COM"
   data          = "/home/hleonov/Data/"
   local_whatif  = "/home/hleonov/Data/whatif/"
   mybin         = "/home/hleonov/bin/"
   amberFF       = "/home/hleonov/Data/forcefields/amber99sb-ildn.ff"
   job['pr']     = data+"runjob_pr.sh"
   job['gpr']    = data+"runjob_gpr.sh"
   job['md']     = data+"runjob.sh"
   mdp['pr']     = data+"pr_conP.mdp"
   mdp['gpr']    = data+"gradual_pr.mdp"
   mdp['md']     = data+"md.mdp"
   
   
#########################################################################################
# Use whatif to protonate residues in the protein and then my script to change his name #
# to be suitable for GROMACS, Amber forcefield                                          #
#########################################################################################
def protonate_whatif(init_file):  
   #run whatif
   tmp_output = init_file[:-4]+"_h.pdb"
   output = init_file[:-4]+"_his_names.pdb"
   pu.run_command("cp "+local_whatif+"TOPOLOGY.* .")
   if (os.path.exists(tmp_output)):
      pu.run_command("rm "+tmp_output)
   commands = "getmol "+init_file+"\n\n\naddhyd\nall 0\nmakmol\n"
   commands +=          init_file+"\n"+tmp_output+"\nall 0\n\n\nexit\ny\n"
   pu.run_command("printf \""+commands+"\"| "+whatif)

   #cleanup
   pu.run_command("rm ALTERR.LOG PDBFILE* WHATIF.FIG fort.1 TOPOLOGY.*")
   #change his names for gromacs
   pu.run_command("python "+mybin+"rewrite_his.py "+tmp_output+" "+output)

   return output


#########################################################################################
# Create topology for the protein by using pdb2gmx. Then seperate it into top and itp   #
# files, and replace the forcefield directory with my local one (for Joung param.)      #
#########################################################################################  
def create_topology(infile, vsites, noignh):
   print "creating Topology for: "+infile
   output = "from_gmx.pdb"
   gmx_inp ="printf \"1\n5\n\" | " #amber
   #gmx_inp = "printf \"2\n1\n\" | " #charmm
   if (noignh is True):
      print "Note: Using pre-determined (probably whatif) protonation."
      pu.run_command(gmx_inp+"pdb2gmx -f "+infile+" -o "+output+" -p prot.top -i prot_pr.itp -vsite "+vsites)
   else: 
      print "Note: Using GROMACS (pdb2gmx) protonation"
      pu.run_command(gmx_inp+"pdb2gmx -f "+infile+" -o "+output+" -p prot.top -i prot_pr.itp -ignh -vsite "+vsites)
   
   #seperate to top and itp files, add Joung-corrected forcefield, and include protein topology
   flag=0
   lines = open("prot.top").readlines()
   top = open("system.top","w")
   itp = open("prot.itp","w")
   for l in lines:
      #For amber: include Joung param in forcefield, and include protein topology
      if(l.find('amber99sb-ildn.ff') != -1):
         l=l.replace("amber99sb-ildn.ff", amberFF )
      if (flag==0):  #before [atoms]
         if(l.find('[ moleculetype ]') != -1):         
            flag=1         
         else:
            print >>top, l,
      if (flag==1):  #in protein topology and before posres include
         if(l.find("Include water topology") != -1):
            flag=2
            print >>top, "#include \"prot.itp\""
         else:
            print >>itp, l,
      if (flag==2):  #after posres include
         print >>top, l,
   return output

#########################################################################################
# adds box,solvent, ions                                                                #
#########################################################################################

def add_box(infile):
   pu.run_command("editconf -f "+infile+" -o box.pdb -bt dodecahedron -d 1.2")
   pu.run_command("genbox -cp box.pdb -cs spc216.gro -o solvated.pdb -p "+opt.top_file)
   pu.run_command("grompp -f "+data+"min_steep.mdp -c solvated.pdb -p %s -maxwarn 1" %opt.top_file)
   pu.run_command("printf \"SOL\n\" | genion -o ionized.pdb -p %s -neutral -conc 0.154" %opt.top_file)
   
#########################################################################################
# Run EM on local machine                                                               #
#########################################################################################  
def run_EM(infile, mdp, top):
   name = mdp[:-4].split("_")[1]
   pu.run_command("grompp -f "+mdp+" -c "+infile+" -p "+top+" -maxwarn 1")
   pu.run_command("mdrun -s topol.tpr -deffnm em_"+name+" -v >& em_"+name+".log")
   return "em_"+name+".part0001.gro"

def getJobname(jobname, suffix=""):
   if (jobname is None):
      return os.getcwd().split("/")[-1]+suffix
   else:
      return jobname+suffix

#########################################################################################
# Submit a job to the owl cluster using a jobfile. Options: pr, gpr, md                 #
# If there is a need for a new mdp file or a new jobfile, pass the as input             #
#########################################################################################      
def submit_md(deffnm, infile, mdp, jobname, runjob_file, top, maxspan, days, pd):
   pu.run_command("grompp -maxwarn 1 -f "+mdp+" -c "+infile+" -p "+top)
#   pu.run_command("qsub  "+runjob_file+" "+jobname)
   if (pd):
      pu.run_command("gsub -d %s -m %s -N %s -P -r %s" %(days, maxspan, jobname, deffnm))
   else: 
      pu.run_command("gsub -d %s -m %s -N %s -r %s" %(days, maxspan, jobname, deffnm))
      
############################ Main ################################
usage = "Usage: %prog [options]"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-i", dest="input_pdb", help="Input structure to prep for simulations (pdb, gro)")
parser.add_option("--protonate", dest="whatif", action="store_true", help="Run WHATIF to protonate molecule")
parser.add_option("--gmx", dest="create_top", action="store_true", help="create Topology (pdb2gmx)")
parser.add_option("--noignh", dest="noignh", action="store_true", help="Do not use -ignh (pdb2gmx)")
parser.add_option("-v","--vsites", dest="vsites", default="hydrogens", help="vsite option for pdb2gmx")
parser.add_option("--box", dest="box", action="store_true", help="Add Box, Solvent  & Ions")
parser.add_option("--em", dest="em", action="store_true", help="Run Energy minimization")
parser.add_option("--top", dest="top_file", metavar="FILE", default="system.top", help="Use given top file")
parser.add_option("-a", dest="em_algorithm", default="steep", help="EM algorithm: steep (default), cg, lbfgs")
parser.add_option("--mdp", dest="mdp_file", metavar="FILE", help="Use given mdp file")
parser.add_option("--job", dest="job_file", metavar="FILE", help="Use given job file")
parser.add_option("-r","--run", dest="run", help="Run molecular dynamics (pr, gpr, md)")
parser.add_option("--name", dest="jobname", help="Job name for runjob script")
parser.add_option("--itp", dest="pr_itp", metavar="FILE", help="positional restraints itp file to add a B state to if there's only one")
parser.add_option("-d","--days", dest="days", default="2", help="days in queue")
parser.add_option("-m","--maxspan", dest="maxspan", default="48", help="hours per queue job")
parser.add_option("-P", dest="pd", action="store_true", help="run with particle decomposition")

(opt, args) = parser.parse_args()
print opt, args

global data, whatif, local_whatif, mybin, job, mdp

job = {}
mdp = {}
init()

next_input = opt.input_pdb

if (opt.whatif):
   next_input = protonate_whatif(next_input)
if (opt.create_top):
   next_input = create_topology(next_input, opt.vsites, opt.noignh)
if (opt.box):
   next_input = add_box(next_input)
if (opt.em):
   if (opt.mdp_file is not None):
      next_input = run_EM(next_input, opt.mdp_file, opt.top_file)
   else:   
      mdp = data+"min_"+opt.em_algorithm+".mdp"
      next_input = run_EM(next_input, mdp, opt.top_file)
if (opt.run):
   if (opt.job_file is not None):
      job[opt.run] = opt.job_file
   jobName = getJobname(opt.jobname,suffix="_"+opt.run)
   if (opt.run is "gpr"):
      print "in GPR"
      if (opt.pr_itp is not None):
        pyutil.run_command("perl ~/bin/make_pr_Bstate.pl %s" %opt.pr_itp)
      else:         
         print "Warning: Make sure prot_pr.itp is set with proper B states"   
   if (opt.mdp_file is not None):
      next_input = submit_md(opt.run, next_input, opt.mdp_file, jobName, job[opt.run], opt.top_file, opt.maxspan, opt.days, opt.pd)
   else:
      next_input = submit_md(opt.run, next_input, mdp[opt.run], jobName, job[opt.run], opt.top_file, opt.maxspan, opt.days, opt.pd)
