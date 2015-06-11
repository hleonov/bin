#!/usr/bin/perl -w

$in = shift || die "Usage: perl $0 <posres.itp>\n";
`mv $in $in.backup`;
open(IN, $in.".backup")|| die "Cannot open posres input file $in.backup\n";
open(OUT, ">$in") || die "Cannot open posres.itp input file $in\n";
while($line = <IN>) {
	if ($line =~ m/\;\s+atom\s+type\s+fx\s+fy\s+fz/) {
		chomp $line;
		$line = $line."    fxB    fyB    fzB\n";
	} elsif ($line =~ m/^\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s*$/) {
		chomp $line;
		$line = $line."   0   0   0\n";
	}
	print OUT $line;
} 

close(IN);
close(OUT);
