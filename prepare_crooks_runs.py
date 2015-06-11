#!/usr/bin/env python
import sys, os, shutil
import commands
from glob import glob
from pymacs_beta.forcefield import MDP
from pymacs_beta import *
import pymacs as pm
import pyutil as pu
import re

def init():
   pu.run_command("export GMXLIB=\"/home/hleonov/Programs/gromacs407/share/gromacs/top_amber99sb_gmx40\"")

def run_command( func, string ):
    s = func.__name__+'(): '+ string
    err_file = func.__name__+'_ERROR.log'
    status, out = commands.getstatusoutput( string )
    if status != 0:
        print >>sys.stderr, 'ERROR in %s ' % func.__name__
        print >>sys.stderr, 'Output written to %s'  % err_file
        fp = open(err_file,'w')
        print >>fp, s
        print >>fp, out
        fp.close()
        sys.exit(1)
    else:
        print "%-90s" % s, ': ok'



def make_dir_tree(traj_file, tpr_file, skip = 4, begin=10001, aa_mut=None, newMut=None, oldMut=None):
    # catch all run?.? dirs
    # dirs = glob('run?.?')
    # for d in dirs:
    #    d = "run_"+str(i) 
    #print '\tPreparing run %s' % d 
    #if (not os.path.exists(d)):
    #    os.mkdir(d)
    #os.chdir(d)
    print "Preparing runs..."
    if (not os.path.exists("morphes")):
        os.mkdir('morphes')  
    
    cmd = 'echo 0| trjconv -f %s -s %s -skip %d -b %d -sep -o morphes/frame.pdb' % (traj_file, tpr_file, skip, begin)
    run_command(make_dir_tree, cmd )
    os.chdir('morphes')
        
    # put each file into single directory 
    # swap residue names if needed for one-for-all equilibration of lambda0    
    file_list = glob('frame*.pdb')
    
    #if the transition is mutating an amino acid
    if ((aa_mut is not None) and (newMut > -1) and (newMut > -1)): 
       cmd = "grep \"2"+aa_mut+" \" ../../mut_"+ str(newMut) +aa_mut+"_gmx.pdb | grep CA | awk '{print $4}'"
       newResName = pu.run_command(cmd)[0].rstrip()
       cmd =  "grep \" "+str(oldMut)+" \" ../../mut_"+ str(newMut) +aa_mut+"_gmx.pdb | grep CA | awk '{print $4}'"
       oldResName = pu.run_command(cmd)[0].rstrip()

    for f in file_list:
       print f
       name = f.split('.')[0]
       if (not os.path.exists(name)):
          os.mkdir(name)
          shutil.copy(f,name+"/"+f)
       elif (os.path.isdir(name)):
          shutil.copy(f, name+"/"+f)
       else: 
          shutil.copy(f, name+"/"+f)
       if ( (oldMut > -1) and (newMut > -1) ): 
       	  print "changing residue names! old="+str(oldMut)+" , new="+str(newMut)
          os.chdir(name)
          pm.model = pm.Model(f)		#load frame 
          pm.model.residues[oldMut-1].set_resname(str(oldResName))

          #pm.model.residues[0].set_resname("NGLU")     
          ##new residue will be re-mutated with Daniel's script due to topology disorder
          #pm.model.residues[newMut-1].set_resname(str(newResName))          
          #pm.model.residues[523].set_resname("CTHR")          
          #  del pm.model.residues[140]['1DHB']      
              
       	  pm.model.writePDB("tmp.pdb")
       	  #add the box dimensions line, extract only the protein and add solvent from the original frame, due to buggy renumbering by pymacs writer, over residue 9999 
       	  pu.run_command("grep CRYST ../"+f+" > cryst")
          pu.run_command("cat cryst tmp.pdb | grep -v \" SOL \" | grep -v \" Na\" | grep -v \" Cl\" > tmp2.pdb")
	  pu.run_command("grep \"SOL\\|Na\\|Cl\\|TER\\|ENDMDL\" ../"+f+" > solvent.pdb")
          pu.run_command("cat tmp2.pdb solvent.pdb > wrong_"+f)
 	  	  
###rm?	  pu.run_command("printf \"19\n\" | trjconv -f wrong_"+f+" -o wrong.gro -n ../../../system.ndx")
	  os.chdir('..')
          pu.run_command("tar -xzf /home/hleonov/Programs/gromacs407/share/gromacs/top/ffamber99sb_mut.tgz")

	  pu.run_command("printf \""+str(newMut)+"\n"+aa_mut+"\nn\n\" | python ~/bin/mutate_beta.py -f "+name+"/wrong_"+f+" -o "+name+"/mut.pdb")
          os.chdir(name)
	  start, err = pu.run_command("cat wrong_"+f+" | grep -v \"SOL\" | awk '{if ($5=="+str(newMut)+") print NR}' | head -1")
	  end, err = pu.run_command("cat wrong_"+f+" | grep -v \"SOL\" | awk '{if ($5=="+str(newMut)+") print NR}' | tail -1")
	  start = int(start)-1
	  end = int(end)+1
	  pu.run_command("head -"+str(start)+" wrong_"+f+" > head.pdb")
	  pu.run_command("cat mut.pdb | awk '{if ($5=="+str(newMut)+") print $0}' | grep -v \"SOL\" > mid.pdb")
	  pu.run_command("tail -n +"+str(end)+" wrong_"+f+" > end.pdb")
	  pu.run_command("cat head.pdb mid.pdb end.pdb > "+f)
 	  pu.run_command("rm tmp*.pdb cryst solvent.pdb head.pdb end.pdb")
          os.chdir('..')
    pu.run_command("rm frame*.pdb")
    pu.run_command("rm ffamber99sb.*")
    os.chdir('..')
    os.chdir('..')
    
def prepare_mdp_files( template_mdp, sw_time, sc_alpha, sc_sigma ):

    mdp = MDP().read( template_mdp )
    
    # make 0->1 file
    mdp['free-energy'] = 'yes'
    mdp['init-lambda'] = 0
    # time is in ps -> convert to steps
    nsteps = sw_time*500
    mdp['nsteps'] = int(nsteps)
    # delta lambda = 1./nsteps
    delta_lambda = 1./float(nsteps)
    mdp['delta-lambda'] = delta_lambda
    mdp['sc-alpha'] = sc_alpha
    mdp['sc-sigma'] = sc_sigma
    
    fp = open('crooks_TI_runA.mdp','w')
    print >>fp, "continuation = no"
    print >>fp, mdp
    fp.close()
    
    # make 1->0 file
    mdp['init-lambda'] = 1
    mdp['delta-lambda'] = -delta_lambda
    
    fp = open('crooks_TI_runB.mdp','w')
    print >>fp, "continuation = no"
    print >>fp, mdp
    fp.close()
    

def make_run_input_files(d, top_file):

#    dirs = glob('run?.?')
#    for d in dirs:
     #here we are at mut_XXX directory
     print '\n\tPreparing run input files %s' % d
     mdp_file = None
     if (d[-1] == 'A' or d[-1] == '0'):          #d.split('.')[0][-1] == 'A': # 0->1
         mdp_file = d+'/crooks_TI_runA.mdp'
     elif (d[-1] == 'B' or d[-1] == '1'):        #d.split('.')[0][-1] == 'B': # 1->0
         mdp_file = d+'/crooks_TI_runB.mdp'
     print "current directory: "+os.getcwd()
     path_to_search = os.path.join(d,'morphes')+'/frame*'
     print "searching for dirs: "+path_to_search
     match = re.search( "(lambda\d)" , os.getcwd())
     if (match):
        os.chdir("../")

     dir_list = glob(path_to_search)
     print dir_list
     for f in dir_list:
        print f
        fname = os.path.basename(f)+'.pdb'
        gro_file = os.path.join(f,fname)
        #cmd = '$gmx407/grompp -f %s -c %s -p %s -o %s/topol.tpr' % (mdp_file, gro_file, top_file, f)
        cmd = 'grompp -f %s -c %s -p %s -o %s/topol.tpr -maxwarn 2' % (mdp_file, gro_file, top_file, f)
        print cmd
        run_command( make_run_input_files, cmd )
        os.system( 'rm mdout.mdp ')
            
        
def main(argv):

    version = "1.0"

    options = [
        Option( "-sw_time", "real", 50., "switching time to use for FGTI runs"),
        #sc-alpha default: 0.3, sigma: 0.25? 
        Option( "-sc_alpha", "real", 0.5, "soft-core alpha"),
        Option( "-sc_sigma", "real", 0.3, "soft-core sigma"),
        Option( "-skip", "int", 10, "skip # frames"),
        Option( "-beg", "int", 10001, "begin at this time(ps)"),
        Option( "-old", "int", -1, "equilibration done & taken from this mutation number (diff directory)"),
        Option( "-new", "int", -1, "crooks should be done with a mutation in this residue number"),
        Option( "-resMut", "string", "A", "Performing an amino acid mutation, options are amino acid letter (one code) or 'no'"),
##         Option( "-box_size", "float", 1.2, "distance from solute to box"),
##         Option( "-conc", "float", 0.15, "ion concentration"),
##         Option( "-vsite", "bool", False, "use virtual sites"),
##         Option( "-fe.dti", "bool", False, "do discrete TI setup"),
##         Option( "-dti_run_time", "float", 10., "Simulation time [ns] at each lambda point"),
##         Option( "-fe.crooks", "bool", False, "do Crooks setup"),
##         Option( "-n_crooks_runs", "int", 1, "setup # crooks runs for each mutation"),
##         Option( "-crooks_run_time", "float", 50., "Simulation time [ns] for crooks equilibrium runs"),
##         Option( "-skip_md_setup", "bool", False, "skip md setup and use -f as starting configuration for free energy runs")
        
        ]
    
    files = [
        FileOption("-d", "r",["dir"],"", "directory with equlibrated states"),
        FileOption("-mdp", "r",["mdp"],"TI_template.mdp", "template TI mdp file"),
        FileOption("-trj", "r",["trr","xtc"], "eq.xtc", "trajectory file name of equilibrated state"),
        FileOption("-top", "r",["top"], "newtop.top", "topology file for the system with a transition"),
        FileOption("-tpr", "r",["tpr"], "eq.tpr", "tpr file needed for trjconv to extract frames from the eq. trajectory"),
##         FileOption("-crooks_mdp", "r/o",["mdp"],"template", "template run mdp file for crooks equilibrium runs"),
##         FileOption("-dti_mdp", "r/o",["mdp"],"template", "template run mdp file for discrete TI calculations"),
##         FileOption("-min_mdp", "r/o",["mdp"],"em", "template minimization mdp file ( for TI or Crooks )"),
##         FileOption("-lambda_steps", "r/o",["txt"],"lambda_steps", "text file with lambda steps for DTI runs"),
        ]
    
    
    
    help_text = ("Script for setting up plain FGTI runs",
                 )

    cmdl = Commandline( argv, options = options,
                        fileoptions = files,
                        program_desc = help_text,
                        check_for_existing_files = False, version = version)


    here = os.getcwd()

    run_dir = cmdl['-d']
    mdp_file = cmdl['-mdp']
    trj_file = cmdl['-trj']
    print trj_file, 
    top_file = cmdl['-top']
    tpr_file = cmdl['-tpr']
    old_mut = cmdl['-old']
    new_mut = cmdl['-new']
    sw_time = cmdl['-sw_time']
    sc_alpha = cmdl['-sc_alpha']
    sc_sigma = cmdl['-sc_sigma']
    res_mut = cmdl['-resMut']
    
    init()
    print '\n\t Preparing FGTI runs in directory..: %s' % run_dir
    print '\t Template mdp file to use............: %s' % mdp_file
    print '\t Switching time to use...............: %8d ps' % int( sw_time )
    print '\t Soft-core alpha to use..............: %8.3f' % sc_alpha
    print '\t Soft-core sigma to use..............: %8.3f' % sc_sigma
    print '\t Equil. lambda0, prev mutated residue to change to standard ...: %s' % old_mut
    print '\t Equil. lambda0, replace standard with hybrid name for residue : %s' % new_mut
    print '\n'
    sys.stdout.flush()

    os.chdir( run_dir )  
    print os.getcwd()
    print '\t Preparing mdp input files........... '
   
    prepare_mdp_files( mdp_file, sw_time, sc_alpha, sc_sigma )
    
    
    print '\t Preparing directory tree............ '
    make_dir_tree(trj_file, tpr_file,  skip = cmdl['-skip'], begin = cmdl['-beg'], aa_mut=res_mut , newMut=new_mut, oldMut=old_mut)
    make_run_input_files(run_dir, top_file)
    os.chdir( here )
    print '\n\t............... DONE .................\n'


if __name__=='__main__':
    
    main( sys.argv ) 
