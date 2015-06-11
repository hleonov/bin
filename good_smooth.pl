#!/usr/bin/perl -w

$in		= shift || die "Usage: perl $0 <in> <out> <odd_avg_num>\n";
$out 	= shift || die "Usage: perl $0 <in> <out> <odd_avg_num>\n";
$param 	= shift || die "Usage: perl $0 <in> <out> <odd_avg_num>\n";

undef(@x);
undef(@y);
undef(@y_out);

read_data();

for my $k ( 1 .. $#y) {
	smooth($k);
}

open(OUT, ">$out") || die "Cannot open file for writing: $out\n";
for my $i (1 .. $#x) {
	print OUT "$x[$i]\t";
	for my $k (1 .. $#y) {
		printf OUT "%3.3f\t", $y_out[$k][$i];
	}
	print OUT "\n";
}
close(OUT);

# smooth(ycol)
sub smooth {
	my $col   = $_[0];
	my $half = int($param/2);
	my $N =$#x;
	# beginning
	for my $i ( 1 .. $half ) {
		$y_out[$col][$i] = average(\@{$y[$col]}, 1, $i+$half);
	}
	# middle
	for my $i ( $half+1 ..$N-$half-1 ) {
		$y_out[$col][$i] = average(\@{$y[$col]}, $i-$half, $i+$half);
	}
	# end
	for my $i ($N-$half .. $N) {
		$y_out[$col][$i] = average(\@{$y[$col]}, $i-$half, $N);
	}
}

# average(array_ref, k1, k2)
# Returns the sum of arrays k1 ... k2 in the given array
sub average {
	my $sum = 0;
	my @arr = @{$_[0]};
	for my $i ($_[1] .. $_[2]) {
		$sum += $arr[$i];
	}
	return $sum / ($_[2] - $_[1] +1); 
}

# read_data(filename)
# for both cols and rows: stored from the 1st place in the array, not 0th !
sub read_data {
	open (IN, $in) || die "Cannot read input file $in\n";
	my $r=1;
	while (<IN>) {
		#skip for comments and non-data lines
		next unless (m/^(\d+)/);

		#Account for multi-columns (one x, many y's)
		my @parts = split(/\s+/, $_);
		$x[$r] = $parts[0];
		for my $k ( 1 .. $#parts) {		
			$y[$k][$r] = $parts[$k];	
		}
		$r++;
	}
	print "number of y cols: $#y \n";
	print "number of data points: $#x \n";
	close(IN);
}
