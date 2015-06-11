#!/usr/bin/perl -w 

# remove first X ps from the simulations. 

$in = shift || die "Usage: perl $0 <infile> <ps_num>\n";
$ps = shift || die "Usage: perl $0 <infile> <ps_num>\n";

open(IN, $in) || die "Cannot open $in\n$!\n";
$in =~ s/\.pdo//;
open(OUT, ">$in\_wo$ps.pdo") || die "Cannot open output file for writing $in\_wo$ps.pdo\n$!\n";
while (<IN>) {
	if (m/^(\w\S+)\s+\S+/) {
		$cur_ps = $1;
		print OUT $_ if ($cur_ps > $ps);
	} else {
		print OUT $_;
	}
}
close(OUT);
close(IN);
