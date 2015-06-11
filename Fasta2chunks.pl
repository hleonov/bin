#!/usr/bin/perl -w

#input fasta files in regular (not one line) form.
#outputs files as the number of sequences, with the name.

my $head = undef;
my $body = "";
$infile = shift || die "Usage: perl $0 <infile> <num_in_chunk>\n";
$chunk_size = shift || die "Usage: perl $0 <infile> <num_in_chunk>\n";

$to_print = "";
$n = 0;
$chunk = 1;

open (IN, "$infile");
while (<IN>) {
  if (m/^>/) {
    $n++;
	if ($n % $chunk_size == 0) {
		print_out($to_print);
		$to_print = "";
	}
	$to_print .= $_;
	
  } else {
    $to_print .= $_;
  }
}
print_out($to_print);

sub print_out {
  my $file = "$infile.part.$chunk.fasta";
  open(OUT,">$file");
  print OUT "$to_print";
  close(OUT);
  $chunk++;
}
