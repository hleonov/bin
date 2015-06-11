#!/usr/bin/perl -w

# This script will seperate a pdb file with many models into single model files.

$infile = shift || die "Usage: perl $0 <pdb-file>";

$out = $infile;
$out =~ s/\.pdb//;
$curr_model = "";
$n = 1;

open(IN,"$infile") || die "Cannot open pdb input $infile\n";
while (<IN>) {
	if (m/ENDMDL/) {
		write_pdb($curr_model, "$out\_$n.pdb");
		$curr_model = "";
		$n++;
	} else {
		$curr_model .= $_;
	}	 
}
close(IN);

sub write_pdb {
	open(OUT, ">$_[1]");
	printf OUT $_[0];
	close(OUT);
}
