#!/usr/bin/env python
import sys, os, shutil
import commands
from glob import glob
from pymacs_beta.forcefield import MDP
from pymacs_beta import *
import pyutil as pu
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



def make_dir_tree(skip = 4):

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
    
    cmd = 'echo 0| trjconv -f ./eq.xtc -s ./eq.tpr -skip %d -b 5001 -sep -o morphes/frame.gro' % skip
    run_command(make_dir_tree, cmd )
    os.chdir('morphes')
        
    # put each file into single directory
    file_list = glob('frame*.gro')
    for f in file_list:
       name = f.split('.')[0]
       if (not os.path.exists(name)):
          os.mkdir(name)
          shutil.move(f,name)
       elif (os.path.isdir(name)):
          os.rename(f, name+"/"+f)
       else: 
         os.rename(f, name)
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
    print >>fp, mdp
    fp.close()
    
    # make 1->0 file
    mdp['init-lambda'] = 1
    mdp['delta-lambda'] = -delta_lambda
    
    fp = open('crooks_TI_runB.mdp','w')
    print >>fp, mdp
    fp.close()
    

def make_run_input_files(d, top_file):
    
#    dirs = glob('run?.?')
#    for d in dirs:
#        print '\n\tPreparing run input files %s' % d
     mdp_file = None
     if (d[-1] == 'A' or d[-1] == '0'):          #d.split('.')[0][-1] == 'A': # 0->1
         mdp_file = d+'/crooks_TI_runA.mdp'
     elif (d[-1] == 'B' or d[-1] == '1'):        #d.split('.')[0][-1] == 'B': # 1->0
         mdp_file = d+'/crooks_TI_runB.mdp'
     dir_list = glob(os.path.join(d,'morphes')+'/frame*')
     for f in dir_list:
        print f
        fname = os.path.basename(f)+'.gro'
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
        Option( "-sc_alpha", "real", 0.3, "soft-core alpha"),
        Option( "-sc_sigma", "real", 0.25, "soft-core sigma"),
        Option( "-skip", "int", 10, "skip # frames"),
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
        FileOption("-top", "r", ["top"], "newtop.top", "topology file for the system with a transition"),
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
    sw_time = cmdl['-sw_time']
    sc_alpha = cmdl['-sc_alpha']
    sc_sigma = cmdl['-sc_sigma']
    top_file = cmdl['-top']
    init()
    print '\n\t Preparing FGTI runs in directory..: %s' % run_dir
    print '\t Template mdp file to use............: %s' % mdp_file
    print '\t Switching time to use...............: %8d ps' % int( sw_time )
    print '\t Soft-core alpha to use..............: %8.3f' % sc_alpha
    print '\t Soft-core sigma to use..............: %8.3f' % sc_sigma
    print '\n'
    sys.stdout.flush()

    os.chdir( run_dir )  

    print '\t Preparing mdp input files........... '
   
    prepare_mdp_files( mdp_file, sw_time, sc_alpha, sc_sigma )
    
    
    print '\t Preparing directory tree............ '
    make_dir_tree( skip = cmdl['-skip'])
    make_run_input_files(run_dir, top_file)
    os.chdir( here )
    print '\n\t............... DONE .................\n'


if __name__=='__main__':
    
    main( sys.argv ) 
