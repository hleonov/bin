#!/usr/bin/perl -w

$cfile = shift || die "Usage: perl $0 <coord-file> <freefile> <outfile>\n";
$freefile = shift || die "Usage: perl $0 <coord-file> <freefile> <outfile>\n";
$out = shift || die "Usage: perl $0 <coord-file> <freefile> <outfile>\n";

my $cdata = read_freefile($freefile);
read_coords($cfile);


open (OUT, ">$out") || die "Cannot open $out\n$!";

print OUT "\@    xaxis  ticklabel on\n";
print OUT "\@    xaxis  ticklabel format general\n";
print OUT "\@    xaxis  ticklabel prec 5\n";
print OUT "\@    xaxis  ticklabel append \"\"\n";
print OUT "\@    xaxis  ticklabel start type auto\n";
print OUT "\@    xaxis  ticklabel stop type auto\n";
print OUT "\@    xaxis  tick place both\n";
print OUT "\@    xaxis  tick spec type both\n";
print OUT "\@    xaxis  tick spec 30\n";
print OUT "\@    xaxis  ticklabel layout vertical\n";

for (my $i=0; $i<=$#aa; $i++) {
	print OUT "\@    xaxis  tick major $i, $z[$i]\n";
	print OUT "\@    xaxis  ticklabel $i, \"$z[$i]  $aa[$i]\"\n";
}
print OUT "\@target G0.S0\n\@type xy\n";
print OUT "$cdata\n";
close(OUT);

sub read_freefile {
	my $data = "";
	open(IN, $_[0]) || die "Cannot open $_[0]\n$!";
	while (<IN>) {
		if (m/^(\d+\.\d+)\s+(\d+\.\d+)/) {
			my $ang = $1*10;
			$data .= "$ang $2\n";
		}
	}
	close(IN);
	return $data;
}

sub read_coords {
	my $i=0;
	open(IN, $_[0]) || die "Cannot open $_[0]\n$!";
	while (<IN>) {
		if (m/\d+\s+(\w{3})\s+(\d+\.\d{3})/) {
			$aa[$i] = $1;
			$z[$i] = $2;
			$i++;
		}
	}
	close(IN);
}
