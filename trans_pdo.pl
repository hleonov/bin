#!/usr/bin/perl -w
# Reads a .pdo output from umbrella sampling and turns it into WHAM input.
# Output: A file with two columns. (1) timestep (2) system coordinates along the reaction coordinate

$pdo = shift || die "Usage:  perl $0 <pdo-file> <ppa-file>\n";
$ppa = shift || die "Usage:  perl $0 <pdo-file> <ppa-file>\n";

$npdo = $1 if ($pdo =~ m/(\S+)\.pdo/);
#Pos1                     = 3.29664268494 3.15019283295 2.39795703888
@pos = split(/\s+/, `grep Pos1 $ppa`);
@k = split(/\s+/,`grep K1 $ppa`);
open(META, ">>meta.dat") || die "cannot open meta.dat\n $!";
print META "$npdo.npdo\t$pos[4]\t$k[2]\n"; 

open(IN, $pdo) || die "Cannot open file $pdo\n $!";
open(OUT, ">$npdo.npdo");
while (<IN>) {
	if (m/(^\d+\.\d+)\s+(-?\d+\.\d+)/) {
		my $time = $1;
		my $relative = $2;
		my $abs = $pos[4]+$relative;
		print OUT "$time\t$abs\n";
				
	} else { 
		print OUT $_;
	}
}
close (IN);
close(OUT);
