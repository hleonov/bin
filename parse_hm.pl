#!/usr/bin/perl -w
# 
# Reads Phobius prediction (TMs and signal peptides). 
# Extract sequences without signal peptide (cut it when exists)

$file = shift || die "Usage: perl $0 <short-phobius> <fasta> <output>\n";
$seq_file = shift || die "Usage: perl $0 <short-phobius> <fasta> <output>\n";
$outfile = shift || die "Usage: perl $0 <short-phobius> <fasta> <output>\n";

#$file = "sp_subcellular_secreted+periplasm.txt.30.0.6.nr.phobius"; 
#$seq_file = "sp_subcellular_secreted+periplasm.txt.30.0.6.nr";

undef(%seq_hash);
undef(%header_hash);

read_fasta($seq_file);
open(INFO, $file) || die "cannot open file for reading: $file\n $!"; 
open(SEQ, ">$outfile") || die "Cannot write seq file\n $!";

<INFO>;
while (<INFO>) {
	# match bitopic with a signal sequence:
	# Q14761|PTCA_HUMAN               1  Y n3-15c20/21o36-54i
    # match bitopic without a signal sequence:
	# P18031|PTN1_HUMAN               1  0 o409-430i
	if (m/(\S+)\s+\d+\s+(Y)\s+.*\/(\d+)/ ||
		m/(\S+)\s+\d+\s+(0)\s+(.)/){
		my $id = $1;
		my $signal = $2;
		my $start = $3;	
		my $seq = $seq_hash{$id};
		if ($signal eq 'Y') {	
			$seq = substr($seq_hash{$id}, $start-1);
		} 
		print SEQ "$header_hash{$id}\n$seq\n";	
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
