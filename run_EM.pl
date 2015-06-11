#!/usr/bin/perl -w

$minZ = shift;
$maxZ = shift;

foreach $file (`ls -1 center*.pdb`) {
	$file =~ s/\s//g;
	$base = $file;
	$base =~ s/\.pdb//;
	$num = $1 if ($base =~ m/_(-?\d+)/);
	next if (defined $minZ && defined $maxZ && ($num>$maxZ || $num<$minZ));
	#$tpr = $base."_for_EM.tpr";
	$tpr = "for_EM.tpr";
	print "$file\t$tpr\n";
	$output = "from_EM_".$base;
	if (system("grompp -f ../em.mdp -c $file -p ../system.top -o $tpr") == 0) {
		if (system("mdrun -v -s $tpr -o $output.trr -c $output.pdb -e $output.edr -x $output.xtc >& $base.logfile.em") == 0) {
			print "EM completed\n";
		} else {
			die "EM mdrun failed\n";
		}

	} else {
		die "Error: EM grompp failed\n";
	}
}
#
