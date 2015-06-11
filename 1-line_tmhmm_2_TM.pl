#!/usr/bin/perl -w
# 
# Reads TMHMM short-format prediction
# Extract TMs in Fasta format, while keeping full header.

$file = shift || die "Usage: perl $0 <short-tmhmm> <fasta> <output>\n";
$seq_file = shift || die "Usage: perl $0 <short-tmhmm> <fasta> <output>\n";
$outfile = shift || die "Usage: perl $0 <short-tmhmm> <fasta> <output>\n";

#$file = "sp_subcellular_secreted+periplasm.txt.30.0.6.nr.phobius"; 
#$seq_file = "sp_subcellular_secreted+periplasm.txt.30.0.6.nr";

undef(%seq_hash);
undef(%header_hash);

read_fasta($seq_file);
open(INFO, $file) || die "cannot open file for reading: $file\n $!"; 
open(SEQ, ">$outfile") || die "Cannot write seq file\n $!";

while (<INFO>) {
	# match bitopic with a signal sequence:
	#gi|222148130|ref|YP_002549087.1|        len=223 ExpAA=84.06     First60=23.87   PredHel=4       Topology=o20-42i63-85o144-166i187-209o
	my @parts = split(/\s+/); 
	$parts[5] =~ s/Topology=//;
	my @TMs = split(/[io]/,$parts[5]);
	foreach $tm (@TMs) {		
		next if ($tm eq "");
		($start, $end) = split(/-/,$tm);
		$seq = substr($seq_hash{$parts[0]}, $start-1, $end-$start+1);
		print SEQ ">$header_hash{$parts[0]}\n$seq\n";
	}
}	

close(INFO); 
close(SEQ);

sub read_fasta {
	open(FAS, $_[0]) || die "cannot open fasta file for reading: $file\n $!"; 
	$header = "";
	$seq = "";
	while (<FAS>) {
		chomp;
		if (m/^>/) {
			if ($header ne "") {
				$seq_hash{$header} = "$seq";
			}
			$seq = "";
			$header = $1 if ($_ =~ m/^>(\S+)/);	#keep sequence ID
			$header_hash{$header} = $_; 	#keep all annotation
		} else {
			$seq .= $_;
		}
	}
	$seq_hash{$header} = "$seq";
	close(FAS);
}
