#!/usr/bin/perl -w
use strict;
use Cwd;
use Math::Complex;
use POSIX qw(ceil floor);
my $usage =
"Usage: $0 {[tprfile1] ... [tprfileN]}\nStarts GROMACS jobs for all .tpr files according to the options defined in the head of the script.
If specific .tpr files are specified, start GROMACS jobs for these only.
If no commandline arguments are given, start jobs for all .tpr files found in the current directory.\n";

# run GROMACS on the Hydra Cluster of the Compute Center Garching
# RTU, Nov 2013
#   adapted from RZG example scripts / Andrea's and Carstens's scripts
#
# for each .tpr file found in the current directory:
#   Create job directory containing a link to the .tpr.
#   Write a job script for submiting to the IBM LoadLeveler queuing system
#    The job script contains multiple steps each running for a certain
#     time (the maximum allowed walltime per step).
#    Each step calls a simple run script which runs the actual calculation
#     starting from the results of the previous step.
# NOTE that a checkpoint file from a previous job can be reused.
#      It has to be found in the job directory and has to be named
#      following the pattern "step.nnnn.cpt", where "nnnn" is a four-digit
#      number. If there are multiple such files, the one with the highest
#      number will be used.

################################

my %global = (
              jobprefix => "E1",
              qsub => "llsubmit",
#              gmx_module => "gromacs/4.6.5-gpu",
              cuda_module => "cuda/5.5",
              gmx_cshrc => "/u/cuk/gmx/465-gcc464-ibmmpi13-cuda55/bin/GMXRC",
#              gmx_cshrc => "/u/cuk/gmx/464-gcc464-ibmmpi13-cuda55-compel/bin/GMXRC",
              mpirun => "poe",
              progname => "mdrun_mpi",
#              tunepmename => "g_tune_pme_mpi",
              tunepmename => "",
              job_script_name => "submit.csh",
              run_script_name => "run.csh",
              steps_in_script => 1,
              threads_per_mpiproc => 10,
              mpiproc_per_node    => 2,
#              npme_per_node       => 0,
              nodes               => 1,
	      gpus_per_node       => 2,
              memory              => "4gb",
              benchmark           => "no",
              benchmark_min_nodes => 1,
              benchmark_max_nodes => 32,
#              edi_suffix => "_to_4NPQ",
             );
# GROMACS additional command line options,
# flags can be set as key with empty argument string,
# user defined "cpi, cpo and s will be ignored (set automatically below)
my %gmxopts = (
              maxh         => 24,
              dlb          => "yes",
              dds          => 0.55,
#	      ei           => "*.edi",
#              nsteps    => 2,
#              nsteps    => 10000,
#              noconfout    => "",
#              resethway    => "",
             );
# LoadLeveler options, user defined "shell" and "name" will be ignored (assigned automatically below) 
my %llopts = (
# days+hours:minutes:seconds is the format returned by llq -l, but strangely it is not accepted by llsubmit
#	      wall_clock_limit => sprintf("%d+%02d:%02d:%02d", floor($gmxopts{maxh}/24), floor(($gmxopts{maxh} % 24)), floor((($gmxopts{maxh} * 60) % 60)), ceil((($gmxopts{maxh} * 3600) % 60)) ),
	      wall_clock_limit => sprintf("%d:%02d:%02d", floor($gmxopts{maxh}), floor((($gmxopts{maxh} * 60) % 60)), ceil((($gmxopts{maxh} * 3600) % 60)) ),
              output => "\$(job_name).\$(jobid).out",
              error => "\$(job_name).\$(jobid).err",
              notification => "always",
              job_type => "parallel",
              notify_user => getpwuid($<)."\@rzg.mpg.de",
              "network.MPI" => "sn_all,not_shared,us",
              requirements => "(Feature==\"gpu\")",
              node_usage => "not_shared",
              tasks_per_node => $global{mpiproc_per_node},
              task_affinity => "core($global{threads_per_mpiproc})",
              resources => "ConsumableCpus($global{threads_per_mpiproc}) ConsumableMemory($global{memory})",
              node => "$global{nodes}",
             );

################################

if ($global{steps_in_script} < 1){
 die ("Exiting because no steps are requested through \$global{steps_in_script} = $global{steps_in_script} (has to be > 0)\n");
}
if (defined $global{npme_per_node} && defined $global{mpiproc_per_node} && $global{npme_per_node} >= $global{mpiproc_per_node}){
 die ("Exiting because \$global{npme_per_node} >= \$global{mpiproc_per_node} -- no MPI ranks for calculating particle-particle interactions left.\n");
}
elsif (defined $global{npme_per_node} && $global{npme_per_node} / ($global{mpiproc_per_node} - $global{npme_per_node}) > 1){
 printf STDERR ("WARNING: ratio of PME and PP MPI ranks is unusually high:\n");
 printf STDERR ("\$global{npme_per_node} / (\$global{mpiproc_per_node} - \$global{npme_per_node}) = %.2f\n", $global{npme_per_node} / ($global{mpiproc_per_node} - $global{npme_per_node}));
}

if ( ( !defined $global{gmx_module} || $global{gmx_module} eq "" ) &&
     ( !defined $global{gmx_cshrc}  || $global{gmx_cshrc}  eq "" ) ) {
 printf STDERR ("ERROR: Can not determine location of your GROMACS installation, define either gmx_module or gmx_cshrc in %global\n");
 exit(0);
} elsif (defined $global{gmx_module} && defined $global{gmx_cshrc}) {
 printf STDERR ("WARNING: found \$global{gmx_module} and \$global{gmx_cshrc}, will use the module!\n");
 $global{gmx_cshrc} = "";
}

my @tprfiles = glob("*.tpr");
if (defined $ARGV[0]){
 @tprfiles = @ARGV;
 $_ = $tprfiles[0];
 if (! /^.+\.tpr/ || $ARGV[0] eq "-h"){
  print $usage;
  exit(0);
 }
}

if ($global{benchmark} eq "no"){
 foreach my $tpr (@tprfiles) {
  $_ = $tpr;
  s/\.tpr//g;
  my $name = $_;
  if (defined $global{"edi_suffix"}) {
   my $this_edifile = $name.$global{"edi_suffix"}.".edi";
   $_ = $this_edifile;
   s/^(equi|prod|EQUI|PROD)\.//g;
   s/\.tip\d{1}p//g;
   $this_edifile = $_;
   if ( -e $this_edifile ) {
    $gmxopts{"ei"} = $this_edifile;
    printf( "Found .edi file \"%s\" matching pattern tpr_name\$global{\"edi_suffix\"}.edi\n", $this_edifile);
   }
  }
  printf("Starting GROMACS job for %s ...\n", $name);
  master_sr($name, "", \%global, \%gmxopts, \%llopts);
 }
}
else {
 foreach my $tpr (@tprfiles) {
  $_ = $tpr;
  s/\.tpr//g;
  my $name = $_;
  printf("Starting GROMACS benchmark for %s ...\n", $name);
  foreach my $this_nodes ( 1, 2, 4, 6, 8, 10, 12, 16, 32 ) {
   if ($this_nodes < $global{benchmark_min_nodes} || $this_nodes > $global{benchmark_max_nodes}){
    next;
   }
   printf("Starting benchmark job on %04d nodes ...\n", $this_nodes);
   my %global2 = %global;
   $global2{nodes} = $this_nodes;
   my %llopts2 = %llopts;
   $llopts2{node} = $this_nodes;
   master_sr($name, sprintf("%04d_nodes", $this_nodes), \%global2, \%gmxopts, \%llopts2);
  }
 }
}

exit(0);

#-----------------------------------------------------------------------------------------
# subroutines
#-----------------------------------------------------------------------------------------

sub master_sr {
 my ($name, $dir_suffix, $g_ref, $gmxo_ref, $llo_ref) = @_;
 my %g = %$g_ref; my %gmxo = %$gmxo_ref; my %llo = %$llo_ref;
 my $rootdir = cwd();
 my $jobname = "GMX.".$name;
 $jobname = $g{jobprefix}.".".$name unless ($g{jobprefix} eq "");
 my $dirname = $jobname;
 if ($dir_suffix ne "") {
  $dirname = $dirname.".".$dir_suffix;
 }
 if ( ! -d $dirname ) {
  mkdir ("$dirname");
 }
 # write the job script
 my $job_script = "$dirname/$g{job_script_name}";
 open(JOB, ">$job_script") || die ("Cannot Open File $job_script");
 print JOB "#!/bin/tcsh\n";
 print JOB "# written by $0\n";
 print JOB "#@ shell = /bin/tcsh\n";
 print JOB "#@ job_name = $name\n";
 # dont' forget to add additional options from %llopts
 foreach my $thisopt (sort keys %llo) {
  print JOB "#@ $thisopt = $llo{$thisopt}\n" unless ( $thisopt eq "shell" || $thisopt eq "job_name" || $thisopt eq "network.MPI" );
 }
 print JOB "\n";
 my $last_step = 0;
 # check for preexisting checkpoint files from previous runs
 my @prev_cpts = glob("$dirname/step.????.cpt");
 foreach my $prev_cpt (@prev_cpts) {
  $_ = $prev_cpt;
  s/^.*\/step\.(\d{4})\.cpt/$1/g;
  my $tmpi = sprintf("%d", $_);
  $last_step = $tmpi unless ($tmpi < $last_step);
 }
 for (my $i = 0; $i < $g{steps_in_script}; $i++) {
  print JOB "\n################################################################################\n";
  printf JOB ("#@ step_name = step.%04d\n", ($i + 1));
  printf JOB ("#@ executable = ./%s\n", $g{run_script_name});
  printf JOB ("#@ arguments = %04d %04d\n", $i + $last_step, ($i + $last_step + 1));
  printf JOB ("#@ dependency = (step.%04d == 0)||(step.%04d == CC_REMOVED)\n", $i, $i) unless ($i == 0);
  # In the RZG setup of LoadLeveler, specifying these options is mandatory for each separate task,
  # even though we already specified them globally before.
  foreach my $thisopt ( "node", "tasks_per_node", "task_affinity", "network.MPI" ) {
   print JOB "#@ $thisopt = $llo{$thisopt}\n";
  }
  printf JOB ("#@ queue\n");
  print JOB "################################################################################\n";
 }
 close(JOB);
 chmod(0755, $job_script);
 # write the run script
 my $run_script = "$dirname/$g{run_script_name}";
 open(RUN, ">$run_script") || die ("Cannot Open File $run_script");
 print RUN "#!/bin/tcsh\n";
 print RUN "# written by $0\n";
 print RUN "# adapted from Andrea Vaiana's script, arguments: ID-number of previous check point file, # ID-number of new check point file\n\n";
 print RUN "module load $g{cuda_module}\n" unless ( $g{cuda_module} eq "" );
 print RUN "module load $g{gmx_module}\n\n" unless ( ! defined $g{gmx_module} || $g{gmx_module} eq "" );
 print RUN "source $g{gmx_cshrc}\n\n"       unless ( ! defined $g{gmx_cshrc}  || $g{gmx_cshrc}  eq "" );
 print RUN "setenv PROCS `echo \"$global{nodes} * $global{mpiproc_per_node}\" | bc`\n";
 print RUN "cd $rootdir/$dirname\n";
 print RUN "if ( -e ../$name.tpr && ! -e ./$name.tpr ) ln -sf ../$name.tpr .\n";
 print RUN "set edifiles=`ls .. | grep -E \".edi\" | head -n 1` && if ( \${edifiles} != \"\" ) ln -sf ../*.edi .\n\n";
 print RUN "echo \"Running step \${LOADL_STEP_ID} on master node \${HOST} ...\"\n";
 if (defined $g{tunepmename} && $g{tunepmename} ne ""){
  my $nsteps_tune_pme = 5000;
  if (defined $gmxo{nsteps}){
   $nsteps_tune_pme = $gmxo{nsteps};
   if ($nsteps_tune_pme > 50000) {
    $nsteps_tune_pme = 50000;
   }
  }
  print RUN "if ( ! -e perf.out ) then\n";
  print RUN " $g{mpirun} $g{tunepmename} -s $name.tpr -steps $nsteps_tune_pme -resethway\n";
  print RUN "endif\n";
  print RUN "npme=`cat perf.out | grep \"Best performance was achieved with\" | perl -pe \'s/.*Best performance was achieved with (\\d+) PME nodes.*/\$1/g\'`\n";
  # FIXME: g_tune_pme should also determine the opptimum topology for the mapping of GPUs to PP MPI ranks
  #        the topology would have to be grep'ed from the output and inserted into the call to $g{progname}
  print RUN "$g{mpirun} $g{progname} -cpi step.\$1.cpt -cpo step.\$2.cpt -s $name.tpr -deffnm $name.long -npme \${npme}";
 } else {
  print RUN "$g{mpirun} $g{progname} -cpi step.\$1.cpt -cpo step.\$2.cpt -s $name.tpr -deffnm $name.long";
 }
 # check the partitioning of MPI ranks between PP and PME
 if (defined $g{mpiproc_per_node}) {
  if    (defined $g{npme_per_node}) { $g{npp_per_node} = $g{mpiproc_per_node} - $g{npme_per_node}; }
  elsif (defined $g{npp_per_node}) { $g{npme_per_node} = $g{mpiproc_per_node} - $g{npme_per_node}; }
  else {$g{npme_per_node} = 0; $g{npp_per_node} = $g{mpiproc_per_node};}
  $gmxo{npme} = $g{npme_per_node} * $g{nodes};
 }
 # tell GROMACS how to assign mutiple GPUs to multiple MPI ranks
 if (defined $g{gpus_per_node} && defined $g{npp_per_node} && $g{gpus_per_node} > 1 && $g{npp_per_node} > 1){
  my $mpiproc_per_gpu = floor($g{npp_per_node} / $g{gpus_per_node});
  my $gpu_to_mpiproc_topology = "-gpu_id "; 
  for (my $i = 0; $i < $g{gpus_per_node}; $i++){
   for (my $j = 0; $j < $mpiproc_per_gpu; $j++) {
    $gpu_to_mpiproc_topology = $gpu_to_mpiproc_topology."".$i;
   }
  }
  # all mpiproc need to be assigned a GPU even if the number of MPI ranks is a non-integer multiple of the number of GPUs
  for (my $j = 0; $j < ($g{npp_per_node} % $g{gpus_per_node}); $j++) {
   $gpu_to_mpiproc_topology = $gpu_to_mpiproc_topology."".($g{gpus_per_node} - 1);
  } 
  print RUN " ".$gpu_to_mpiproc_topology unless ($g{tunepmename} ne ""); 
 } 
 # dont' forget to add additional options from %gmxopts
 foreach my $thisopt (sort keys %gmxo) {
  print RUN " -$thisopt $gmxo{$thisopt}" unless ($thisopt eq "cpi" || $thisopt eq "cpo" || $thisopt eq "s" || ($thisopt eq "npme" && $g{tunepmename} ne "") || ((($g{gpus_per_node} > 1 && $g{npp_per_node} > 1) || $g{tunepmename} ne "") && $thisopt eq "gpu_id"));
 }
 print RUN "\n";
 print RUN "# if the previous .cpt already has the full number of simulation steps, mdrun will exit without writing a checkpoint file\n";
 print RUN "if ( ! -e step.\$2.cpt && -e step.\$1.cpt ) then\n";
 print RUN " ln -s step.\$1.cpt step.\$2.cpt\n";
 print RUN "endif\n";
 print RUN "\n";
 close(RUN);
 chmod(0755, $run_script);
 system("cd $dirname && $g{qsub} $g{job_script_name}");
 printf STDOUT ("Started simulation of \"%s\" in <= %d steps of <= %d h duration each, beginning at step %d.\n", $name, $g{steps_in_script}, $gmxo{maxh}, $last_step + 1);
}

sub in_array {
 my ($scalar, $arrayref) = @_;
 my @array = @$arrayref; my $in = 0;
 CHECK: foreach my $element (@array){
  if ($element eq $scalar){$in = 1; last CHECK;}
 }
 return ($in);
}

