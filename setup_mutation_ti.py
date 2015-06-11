import sys, os
import optparse
import pyutil as pu
from glob import glob
import re
import string

# Define global variables
def init():

   #os.environ['GMXLIB'] = str("\"/home/hleonov/Programs/gromacs407/share/gromacs/top_amber99sb_gmx40:/home/hleonov/Programs/gromacs407/share/gromacs/top/\"")   
   
   print "\n",pu.run_command("echo $GMXLIB")
   print "\n",pu.run_command("which grompp")

   global whatif, data, mdpL, local_whatif, mybin, amberFF, job, days, span
   whatif        = "/usr/local/whatif/DO_WHATIF.COM"
   data          = "/home/hleonov/Data/"
   local_whatif  = "/home/hleonov/Data/whatif/"
   mybin         = "/home/hleonov/bin/"
   amberFF       = "/home/hleonov/Data/forcefields/amber99sb-ildn.ff"
   next['pr']     = "confout.part0001.gro"
   next['eq']     = "pr.part0001.gro"
   next['crooks'] = "eq.part0001.gro"
   mdpL['em']     = data+"/mti_mdp/em_lambda.mdp"  
   mdpL['pr']     = data+"/mti_mdp/equil_lambda.mdp"
   mdpL['eq']     = data+"/mti_mdp/crooks_GPU_equilibration_state.mdp"
#   mdpL['eq']     = data+"/mti_mdp/crooks_equilibration_state.mdp"
   mdpL['crooks'] = data+"/mti_mdp/crooks_non_eq_lambda.mdp"
   days['pr']     = 2
   days['eq']     = 16
   span['pr']     = 48
   span['eq']     = 48

   
def grompp(mdp, infile, top, out, useOld=False):
   if (useOld == True):
      grompp = "$gmx407/grompp"
   else: 
      grompp = "grompp"
   out, err = pu.run_command(grompp+ " -f "+mdp+" -c "+infile+" -p "+top+" -o "+out+" -maxwarn 2 >> setup.log")
   exit("Failed grompp command") if (err == -1) else ''   

def prep_mutation(pdb, script=None, rid=None, residue=None):
   err=0
   pu.run_command("rm setup.log")   
   
   ## Prep files
   pu.run_command("echo \"setup.sh log\" > setup.log")
   pu.run_command("echo \"Call: mutate %s %s\" >> setup.log" %(rid, residue))
   pu.run_command("tar -xzf /home/hleonov/Programs/gromacs407/share/gromacs/top/ffamber99sb_mut.tgz")
   pu.run_command("cp ~/Data/topology/residuetypes.dat .")

   #### perhaps add later a short function to protonate. for now, assuming protonated strucutres.
   
   ## Convert to amber format
   pu.run_command("echo \"Converting structure to amber format.\" >> setup.log")
   out, err = pu.run_command("python ~/bin/make_amber.py -f %s -o amber.gro &>> setup.log" %pdb)
   exit() if (err == -1) else ''   
   pu.run_command("mv amber.gro amber.pdb")

   ## Mutate
   pu.run_command("echo \"Mutating structure.\" >> setup.log")
   name = "mut_"
   if (rid is not None) and (residue is not None):
      name += rid+residue
      if not (os.path.exists(name)):
         pu.run_command("mkdir "+name)
      err=pu.run_command("printf \""+rid+"\n"+residue+"\nn\n\" | python ~/bin/mutate_beta.py -f amber.pdb -o "+name+"/"+name+".pdb &>> setup.log")
   else:
      name += "script"
      err=pu.run_command("python ~/bin/mutate_beta.py -f amber.pdb -script "+script+" -o "+name+"/"+name+".pdb &>> setup.log")
   exit() if (err == -1) else ''   
   
   pu.run_command("mv setup.log "+name+"/setup.log")
   dirfile = name+"/"+name

   ## pdb2gmx
   pu.run_command("echo \"Creating topology.\" >> "+name+"/setup.log")
   err=pu.run_command("$gmx407/pdb2gmx -f "+dirfile+".pdb -o "+dirfile+"_gmx.pdb -p "+name+"/topol.top -i "+name+"/posre.itp -ff amber99sb -water spce &>> "+name+"/setup.log")
   print "err code :"+str(err)

   exit() if (err == -1) else ''   

   ## Making B-state
   pu.run_command("echo \"Modifying topology for B-state.\" >> "+name+"/setup.log")
   pu.run_command("cp "+name+"/topol.top .")
   pu.run_command("python ~/bin/make_bstate.py  &>> "+name+"/setup.log")
   com = "sed 's/spce/\.\.\/spce/' newtop.top.ndx | sed 's/ions\.itp/\.\.\/ions\.itp/' | sed 's/mut_.*\/posre\.itp/\.\/posre\.itp/' > "+name+"/newtop.top"
   pu.run_command(com)

   ## Box + Solvate + Ions
   pu.run_command("echo \"defining Box.\" >> "+name+"/setup.log")
   pu.run_command("$gmx407/editconf -f "+dirfile+"_gmx.pdb -o "+dirfile+"_gmx.pdb -d 1.2 -bt dodecahedron  &>> "+name+"/setup.log")
   pu.run_command("echo \"Solvating in SPCE.\" >> "+name+"/setup.log")
   pu.run_command("$gmx407/genbox -cp "+dirfile+"_gmx.pdb -cs spc216.gro -o "+dirfile+"_water.pdb -p "+name+"/newtop.top &>> "+name+"/setup.log")
   pu.run_command("echo \"Adding ions.\" >> "+name+"/setup.log")
   grompp("~/Data/em.mdp", dirfile+"_water.pdb", name+"/newtop.top", name+"/water.tpr", useOld=True)
   #replaced by grompp   pu.run_command("$gmx407/grompp -f ~/Data/em.mdp -c "+dirfile+"_water.pdb -o "+name+"/water.tpr -p "+name+"/topol.top &>> "+name+"/setup.log")
   pu.run_command("echo SOL | $gmx407/genion -s "+name+"/water.tpr -o "+dirfile+"_ions.pdb -p "+name+"/newtop.top -pname NaJ -nname ClJ -conc 0.15 -neutral &>> "+name+"/setup.log")

#   pu.run_command("mv newtop.top.ndx "+name+"/newtop.top")
   if not (os.path.exists(name+"/lambda0") and os.path.exists(name+"/lambda1")):     
      pu.run_command("mkdir "+name+"/lambda0 "+name+"/lambda1")
   pu.run_command("rm "+name+"/#*")
   return name
      
def minimize(next_input, dir, lmb, mdp=None):
   if (mdp is None):
      mdp=mdpL['em']
   mdp=re.sub("\.mdp\Z", str(lmb)+".mdp", mdp) 
   print mdp  
   pu.run_command("echo \"Running EM\" >> "+dir+"/setup.log")
   #grompp(mdp, infile, top, out)
   grompp(mdp, dir+"/"+next_input, dir+"/newtop.top", dir+"/lambda"+str(lmb)+"/em.tpr")
   os.chdir(dir+"/lambda"+str(lmb))
   pu.run_command("mdrun_threads -nt 4 -v -s em.tpr &> em.log")
   os.chdir("../../")

def submit_md(next, mdp, run, d, sp, dir, lmb):
   if (d is None):
      d = str(days[run])
   if (sp is None):
      sp = str(span[run])

   mdp=re.sub("\.mdp\Z", str(lmb)+".mdp", mdp)
   lmb_dir = dir+"/lambda"+str(lmb)
   grompp(mdp,  lmb_dir+"/"+next, dir+"/newtop.top", lmb_dir+"/"+run+".tpr")

   os.chdir(lmb_dir)
   pu.run_command("rm runjob.out*")
   name = re.sub("/", "_", re.sub("/home/hleonov/Research/[a-z]*/?", "", os.getcwd()))
   name = re.sub("mutations_","", name)
   name = re.sub("lambda","", name)
   print "Make sure what cutoff scheme is used (Verlet, Group)? GPUs? CPUs?\n"
   if (run == "pr"):
      pu.run_command("gsub -f "+run+".tpr -N "+run+name+" -r "+run+" -d "+d+" -m "+sp)
   os.chdir("../../")

def rerun_job(next, run, d, sp, dir, lmb):
   if (d is None):
      d = str(days[run])
   if (sp is None):
      sp = str(span[run])
   lmb_dir = dir+"/lambda"+str(lmb)
   os.chdir(lmb_dir)
   pu.run_command("rm runjob.out*")
   #pu.run_command("gsub -f "+run+".tpr -N lmb"+str(lmb)+run+dir+" -r "+run+" -c "+run +".cpt -d "+d+" -m "+sp)
   name = re.sub("/", "_", re.sub("/home/hleonov/Research/[a-z]*/?", "", os.getcwd()))
   name = re.sub("mutations_","", name)
   name = re.sub("lambda","", name)

   #due to Error: Count mismatch for state entry FE-lambda, code count is 7, file count is 1, using another g_submit(.orig)
   if (os.path.exists(run+".cpt")):
      pu.run_command("/home/hleonov/Research/abc/mutations/gsub -f "+run+".tpr -N "+run+name+" -c "+run+".cpt -r "+run+" -d "+d+" -m "+sp)
   else:
      pu.run_command("/home/hleonov/Research/abc/mutations/gsub -f "+run+".tpr -N "+run+name+" -c state.cpt -r "+run+" -d "+d+" -m "+sp)
   os.chdir("../../")

def check_eq(dir, lmb):
   command = "trjcat -f eq.part*.xtc -o eq.xtc" 
   os.chdir(dir+"/lambda"+str(lmb))
   #concat trajectory parts of eq.part*.xtc
   if (not os.path.exists("eq.xtc")):
      pu.run_command(command)
   command = "printf \"0\n\" | trjconv -f eq.xtc -pbc mol -ur compact -o eq_pbc.xtc -s eq.tpr"
   if (not os.path.exists("eq_pbc.xtc")):
      pu.run_command(command)
   #RMSD
   if (not os.path.exists("rmsd.xvg")):
      command = "printf \"3\n3\n\" | g_rms -f eq_pbc.xtc -s pr.part0001.gro"
      pu.run_command(command)   
   #index file of core, core extraction and 2D projection
   pu.run_command("csh ~/Research/abc/mutations/2d_proj.csh ../"+dir+"_ions.pdb eq_pbc.xtc")
   os.chdir("../../")
##### Main ######
usage = "Usage: %prog [options]"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-i", dest="input_pdb", help="Input structure to prep for simulations (pdb, gro)")
parser.add_option("--dir", dest="workdir", help="Input workdir to work at")
parser.add_option("--mut", dest="mutate", action="store_true", help="Setup mutation")
parser.add_option("--resid", dest="rid", metavar="int", help="resid to mutate (numbering should start from 1)")
parser.add_option("--res", dest="res", help="one-letter residue code to mutate to")
parser.add_option("--script", dest="mscript", metavar="FILE", help="mutations script (for more than one)")
parser.add_option("--top", dest="top_file", metavar="FILE", default="system.top", help="Use given top file")
parser.add_option("--em", dest="min", action="store_true" , help="Run minimization for both lambdas")
parser.add_option("--noem", dest="no_em", action="store_true" , help="Do not run minimization at mutation function")
#parser.add_option("--equil", dest="equil", action="store_true" , help="Use given top file")
parser.add_option("--mdp", dest="mdp_file", metavar="FILE", help="Use given mdp file")
parser.add_option("--job", dest="job_file", metavar="FILE", help="Use given job file")
parser.add_option("-r","--run", dest="run", help="Run md or transitions (pr, eq, crooks, eqchk, analyze)")
parser.add_option("--name", dest="jobname", help="Job name for runjob script")
parser.add_option("-s","--skip", dest="skip", default="50", help="Frames skip in trajectory when creating initial states")
parser.add_option("--beg", dest="crooks_beg", default="10000", help="begin crooks at this time (ps)")
parser.add_option("--sw", dest="switch", default="200", help="transition time (ps)")
parser.add_option("-d","--days", dest="days", help="days in queue")
parser.add_option("-m","--maxspan", dest="maxspan", help="hours per queue job")
parser.add_option("-P", dest="pd", action="store_true", help="run with particle decomposition")
parser.add_option("--rerun", dest="rerun_lmb", help="rerun job for the given lambda (either 0 or 1), must be combined with --run option")
parser.add_option("--lmb", dest="lmb", metavar="int", help="work on specific lambda - 0 or 1")
parser.add_option("--old", dest="old", metavar="int", help="equil. from old mutation, change this resid name to standard")
parser.add_option("--new", dest="new", metavar="int", help="equil. from old mutation, change this resid name to a hybrid residue")
parser.add_option("-T", dest="T", metavar="int", default=310, help="Crooks Temperature")
parser.add_option("--bins", dest="bins", metavar="int", default=10, help="Number of bins for crooks analysis")
parser.add_option("--cleanup", dest="cleanup", action="store_true",help="Clean unncesseary files")
parser.add_option("--dhdl", dest="dhdlName", default="md.part0001.xvg",help="dhdl filename in morphes")


(opt, args) = parser.parse_args()
print opt, args
        
global data, whatif, local_whatif, mybin, mdpL, next, days, span

next = {}
mdpL = {}
days = {}
span = {}
init()

next_input = opt.input_pdb

if (opt.mutate):
   if (opt.mscript is not None) or ((opt.rid is not None) and (opt.res is not None)):
      if (next_input[-3:] != "gro"):
         pu.run_command("cat %s | sed 's/OC\([0-9]\)/O\\1 /' > tmp.pdb" %next_input )
         next_input = next_input.replace("pdb","gro")
         pu.run_command("editconf -f tmp.pdb -o %s" %(next_input))      
      name = prep_mutation(next_input, script=opt.mscript, rid=opt.rid, residue=opt.res)
      if (not opt.no_em):
         minimize(name+"_ions.pdb",dir=name, lmb=0)
         minimize(name+"_ions.pdb",dir=name, lmb=1)
      exit(0)      
   else:
      print "Mutation needs either a mutations script file, or a resid + one-letter code residue to mutate to"
      exit(-1)

if (opt.workdir is None):
   print "Needs a directory to be specified for further action. \n",
   print "For minimization, also needs an input file, usually *_ions.pdb. \n"
   exit(-1)
   
if (opt.min):
   if (next_input is None):
     print "For minimization, also needs an input file, usually *_ions.pdb. \n"
     exit(-1)
   else:
      if (opt.lmb is None):
         minimize(next_input, dir=opt.workdir, lmb=0)
         minimize(next_input, dir=opt.workdir, lmb=1)
      else:
         minimize(next_input, dir=opt.workdir, lmb=opt.lmb)

if (opt.run == "pr" or opt.run == "eq"):    
   ### PR or Equilibration ensemble generation (eq)     
   if (next_input is None):
      next_input = next[opt.run]
   if (opt.mdp_file is None):
      mdp = mdpL[opt.run]      
      
   if (opt.rerun_lmb is not None):
      print "Continuation of "+opt.run+" phase (rerun option)"
      rerun_job(next_input, opt.run, opt.days, opt.maxspan, dir=opt.workdir, lmb=opt.rerun_lmb)
   else:   #normal start of either pr or eq phase
      print "echo \"Setting up "+opt.run+" phase"
      if (opt.lmb is None):    #run both lambda 0 and 1
         submit_md(next_input, mdp, opt.run, opt.days, opt.maxspan, dir=opt.workdir, lmb=0)
         submit_md(next_input, mdp, opt.run, opt.days, opt.maxspan, dir=opt.workdir, lmb=1)
      else: #run only one of the lambdas a chosen, good when optimizing large amounts of independent mutations
	      submit_md(next_input, mdp, opt.run, opt.days, opt.maxspan, dir=opt.workdir, lmb=opt.lmb)
elif (opt.run == "crooks"):   
   #Run crooks for both lambdas, according to switch time and skip frames values provided.
   os.chdir(opt.workdir)
   for lmb in (range(2)):    
      mdp = mdpL[opt.run]     
      mdp=re.sub("\.mdp\Z", str(lmb)+".mdp", mdp) 
      if (opt.old is not None) and (opt.new is not None): #only lambda1 is eqiulibrated, 0 is taken from other simulations
         print "replacing residues old %s, new %s" % (opt.old, opt.new) 
         if (lmb == 1):
            command = "printf \"3\n3\n\" | g_rms -f lambda"+str(lmb)+"/eq.xtc -s lambda"+str(lmb)+"/pr.part0001.gro"
            pu.run_command(command)
            pu.run_command("python ~/bin/prepare_crooks_runs.py -d lambda"+str(lmb)+" -top ./newtop.top -mdp "+mdp+" -sw_time "+opt.switch+" -skip "+opt.skip)
         elif (lmb == 0):
            pu.run_command("python ~/bin/prepare_crooks_runs.py -d lambda"+str(lmb)+" -top ./newtop.top -mdp "+mdp+" -sw_time "+opt.switch+" -skip "+opt.skip+" -old "+opt.old+" -new "+opt.new)
      else: #lambda0 and lambda1 are both equilibrated  
         print "Same protocol for lambda0 and lambda1, no need to change residue names"
         command = "printf \"1\n1\n\" | g_rms -f lambda"+str(lmb)+"/eq.xtc -s lambda"+str(lmb)+"/pr.part0001.gro"
         if (not os.path.exists("rmsd.xvg")):
            pu.run_command(command)
         pu.run_command("python ~/bin/prepare_crooks_runs.py -d lambda"+str(lmb)+" -top ./newtop.top -mdp "+mdp+" -sw_time "+opt.switch+" -skip "+opt.skip+" -beg "+opt.crooks_beg )
   #pu.run_command("csh ../../crooks_batch.csh")
   os.chdir("../")   
	
#   lmb=1
#   mdp=re.sub("\.mdp\Z", str(lmb)+".mdp", mdp)
#   pu.run_command("python ~/bin/prepare_crooks_runs.py -d lambda1 -top ./newtop.top -mdp "+mdp+" -sw_time "+opt.switch+" -skip "+opt.skip)

elif (opt.run == "eqchk"):
   if (opt.lmb is None):
      check_eq(opt.workdir, 0)
      check_eq(opt.workdir, 1)
   else:
      check_eq(opt.workdir, opt.lmb)
   #plot pca   
   os.chdir(opt.workdir)
   command = "python ~/Research/abc/scatterPlots_on_main.py -o eq_pca.png --main ~/Research/abc/abc_structures/2d_1-2_all_all.xvg --cfile ~/Research/abc/abc_structures/color --fig_w 12 --fig_h 10 -p ~/Research/abc/mutations/eq_pca.param -s 0.5 -t \" %s \"" % opt.workdir
   pu.run_command(command)
   
   os.chdir("../")
   
elif (opt.run == "analyze"):
   os.chdir(opt.workdir)
   pu.run_command("../../fix_filenames.csh")
   if (not os.path.exists("integ0.dat")):
      pu.run_command("python ~/bin/analyze_crooks_overlap.py -pa lambda0/morphes/frame*/%s -pb lambda1/morphes/frame*/%s -T %s -nbins %s" %(opt.dhdlName, opt.dhdlName, opt.T, opt.bins))
   else:
      pu.run_command("python ~/bin/analyze_crooks_overlap.py -i0 integ0.dat -i1 integ1.dat -T %s -nbins %s" %(opt.T, opt.bins))
   os.chdir("../")

if (opt.cleanup):
   os.chdir(opt.workdir)
   pu.run_command("rm #* lambda*/morphes/#*")
   pu.run_command("clean_run_files.csh lambda0")
   pu.run_command("clean_run_files.csh lambda1")
   pu.run_command("rm lambda*/morphes/frame*/step*")
   pu.run_command("rm lambda*/morphes/frame*/core*")
   pu.run_command("rm lambda*/morphes/frame*/bench*")
   pu.run_command("rm lambda*/morphes/frame*/err*npme*")
   pu.run_command("rm lambda*/morphes/frame*/bench*")
   pu.run_command("rm lambda*/morphes/frame*/#*")
   
