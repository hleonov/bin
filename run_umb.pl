#!/usr/bin/perl -w
undef(%coords);
&read_pull("pull.ppa");
&read_min_coord;
$K = shift || die "Usage: perl $0 <force> <np>"; #spring force
$np = shift || die "Usage: perl $0 <force> <np>\n";
foreach $file (`ls -1 from_EM_center*.pdb`) {
        $file =~ s/\s//g;
        $base = $file;
        $base =~ s/\.pdb//;
	#	$num = $1 if ($base =~ m/(-?\d+\.\d)/);
		$num = $1 if ($base =~ m/(-?\d+)/);
		print "base=$base\tnum=$num\n";
        $pdb = $base.".pdb";
        print "$file\t$pdb\n";
        $output = $base."_pull_k$K";
		&create_pull_config($num);
		unless (-e "$output.pdo") {
          if (system("grompp -f pmf.mdp -c $pdb -n pull.ndx -p ../system.top -o pull_run.tpr -np $np") == 0) {
                if (system("mpirun -np $np mdrun -pi pull_K$K.ppa -po $output.ppa -pd $output.pdo -s pull_run.tpr -o $output.trr -c results.$output.pdb -e $output.edr -x $output.xtc -v -np $np >& $base.logfile.pull") == 0) {
                        print "pull-md completed\n";
                } else {
                        die "pull mdrun failed\n";
                }
          } else {
                die "Error: pull-md grompp failed\n";
          }
		}
}

sub create_pull_config {
	open(OUT, ">pull_K$K.ppa") or die "Cannot open pull_K$K.ppa\n$!";
	print OUT $pull_str;
	print OUT "K1	= $K\n";
	print OUT "Pos1	= $coords{$_[0]}\n";
	close(OUT);
}

sub read_min_coord {
	open(IN,"min_coords.dat");
	while (<IN>) {
		if (m/(\S+)\s+(\S+\s+\S+\s+\S+)/) {
			$coords{$1} = $2;
		}
	}
	close(IN);
}

sub read_pull {
	$pull_str = `grep -v '^Pos1' $_[0] | grep -v '^K1'`;
}
