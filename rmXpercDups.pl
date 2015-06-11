#!/usr/bin/perl -w

use Bio::SeqIO;
$ID_THRESH = shift or die "Usage: perl $0 idFracThresh infile\n";

my $in = Bio::SeqIO->new(-fh => \*ARGV, -format => 'fasta');
my $out = Bio::SeqIO->new(-format => 'fasta');
my @sequs;
undef(@sequs);
while (my $seq = $in->next_seq) {
	my $i=0;
	my $idflag=1;
	while ($idflag && $i<=$#sequs) {
		my $id = getIden($sequs[$i]->seq,$seq->seq);
		$idflag=0 if ($id > $ID_THRESH);
		$i++;
	}
	if ($idflag) {
		$sequs[$#sequs+1] = $seq;
      $out->write_seq($seq);
	}
}

sub getIden {
	return 0 if (length($_[0]) != length($_[1]));
	my @s1 = split(//,$_[0]);
	my @s2 = split(//,$_[1]);
	my $iden = 0;
	for my $i (0 .. $#s1) {
		if ($s1[$i] eq $s2[$i]) {
			$iden++;
		}
	}
	return $iden/($#s1+1);	
}
