#!/usr/bin/perl -w

# gradual_pr.pl
# ==============
# Perform gradual PR - perform a short period of PR for decending force constants until zero is reached.
# Usage example 1: perl gradual_pr.pl pr_1.mdp system_from_EM-2.gro dmpc_pr.itp 1000 10 8
# Usage example 2 (with crash recovery): perl gradual_pr.pl pr_1.mdp system_from_EM-2.gro dmpc_pr.itp 1000 10 8 1
# Notes:
# 1. Must have a valid pr.mdp file 
# 2. Must have lam-mpi active
# 3. Usually perform EM before this step unless resuming from mid-PR
# 4. Must have all itp and top files in the directory. 
 
$pr_file = shift || &print_help;	#mdp file
$em_file = shift || &print_help;	#Initial structure
$itp_list = shift || &print_help;	#	X_pr.itp file of the restrained molecule (dmpc_pr.itp) as written in X.itp file.
$max_f 	 =  shift || &print_help;	#force constant to start with
$min_f   =  shift || &print_help;	#force constant to finish with
$df 	 = shift || &print_help;	#force steps to decend by
$dt	      = shift || &print_help;	#timestep. 'f' for function, constant otherwise     
$np 	 = shift || &print_help;	#number of processors to use
$recover = shift;

$f=$max_f;
undef(@list);
read_itp_list($itp_list);
unless (-e "backup.$pr_file") {
	`cp $pr_file backup.$pr_file`;
}

while ($f>=$min_f) {
	write_itp($f);
	#system("rm \#*");
	$nsteps = compute_nsteps($f);
	change_mdp($nsteps);
	 
	#if it's the first PR and not recovering from a crash
	if ($f == $max_f && !$recover) {
		#for gmx3 add  "-np $np"
		system("grompp -f $pr_file -c $em_file -n system.ndx -p system.top -o system_for_PR.tpr") == 0
			|| die "Cannot execute grompp\n";

	} else {
		$i=$f+$df;
		#for gmx3 add  "-np $np"
		system("grompp -f $pr_file -c system_from_PR$i.gro -n system.ndx -p system.top -o system_for_PR.tpr") == 0
			|| die "Cannot execute grompp\n";
	}
	#gmx3
	#system("mpirun -np $np mdrun -v -s system_for_PR.tpr -o system_from_PR$f.trr -c system_from_PR$f.gro -e system_from_PR$f.edr -x system_from_PR$f.xtc -np $np >& logfile.pr1") == 0
		#|| die "Cannot execute mpirun\n";
		
	system("mpirun -np $np /usr/local/gromacs/bin/mdrun_mpi -v -s system_for_PR.tpr -deffnm system_from_PR$f >& logfile.pr1") == 0
			|| die "Cannot execute mpirun\n";

	$f = $f - $df;
}

sub change_mdp {
	my $num =  $_[0];
	#if (system("sed \'s/^nsteps\(.*\)= \([0-9]*\)/nsteps = $num/\' backup.$pr_file > out") == 0) {
	if (system("sed 's/^nsteps.*/nsteps = $num/' backup.$pr_file > $pr_file") == 0) {
		print "Success\n";
	}
}

# Using function: y = -40x +50000
# to start from F=1000, 20ps and end at F=0, 100ps
sub compute_nsteps {
	if ($dt =~ m/f/) {
		return (-40*$_[0] + 50000);	
		#return (-232*$_[0] + 100000); 
	} else {
		return ($dt/0.002);
	}
}

sub read_itp_list {
	foreach $file (`cat $_[0]`) {
		chomp $file;
		unless (-e "backup.$file") {
			`cp $file backup.$file`;
		}
		push(@list, $file); 
	}
}

sub write_itp {
	my $f = $_[0];
	foreach $pr_itp (@list) {
		next if (!defined $pr_itp);
		open(IN, "backup.$pr_itp");
		open(OUT, ">$pr_itp");
		while ($line = <IN>) {
			if ($line =~ m/(\s+\d+\s+\d+)/) {
				my $begin = $1;
				print OUT "$begin\t$f\t$f\t$f\n";
			} else {
				print OUT $line;
			}
		}
		close(OUT);
		close(IN);
	}
	
}

sub print_help {
	print "Usage: perl $0 <pr-file> <initial-struct> <itp-list> <maxF> <minF> <df> <dt> <np>\n";
	print "Where:\tpr-file : pr.mdp parameters file.\n";
	print "\t\tem-file : Starting structure, preferably after EM\n";
	print "\t\titp-list: A list containing the itp of the restrained molecules (X_pr.itp). One line each\n";
	print "\t\tmaxF : Initial force constant to use in PR\n";
	print "\t\tminF : Final force constant to use (minus 1)\n";
	print "\t\tdf : Force step to decend by\n";
	print "\t\tdt : Timestep (ps). 'f' for pre-programmed function. Constant otherwise.\n";
	print "\t\tnp : Number of processors to use\n";
#	print "\t\trecover: Optional. use to continue from a gradual PR which crashed.\n";
	exit;
}
