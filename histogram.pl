#!/usr/bin/perl -w
use POSIX qw (ceil floor);

undef(@data);
#$min = 42;
#$max = 47;
$min = 9999999;    #really big int
$max = -9999999;   #really small int
$norm = 0;
$infile 	= shift || die "Usage: perl $0 <infile> <bin_num> [norm 1]\n";
$num_bins 	= shift || die "Usage: perl $0 <infile> <bin_num> [norm 1]\n";
$norm = shift;

#read the data file into an array and find min/max
open(IN, $infile) || die "Cannot open $infile\n$!";
$i=0;
while (<IN>) {
	if (m/^\s*(\S+\.?\S*)\s*$/) {
		$x = $1;
		$min = $x if ($x < $min);
		$max = $x if ($x > $max);
		$data[$i] = $x;
		$i++;
	} else {
		die "Error: Data file contain more than one column.\n";
	}
}
close(IN);
$max+=0.000001;	# so that large edge won't get a new seperate bin.
$range = $max - $min;

# Prepare arrays: 
#	1. Construct array of bin limits
# 	2. Initialize the histogram
for $i (0..$num_bins-1) {
    $bins[$i]= $min + ( $range * ($i / $num_bins));
	$histogram[$i] = 0;
}

# Divide data into bins and count
# get statistics
$sum=0;
foreach $x (@data) {
   $bin_ind = int(1.0 * ($x-$min) * ($num_bins / $range) );
   $histogram[$bin_ind]++;
   $sum += $x;
}
$mean = $sum/(1+$#data);

$sqsum=0;
foreach $x (@data) {
   $sqsum += ($x-$mean)**2
}
$stddev = sqrt($sqsum/(1+$#data));
printf STDERR "Mean = %5.6f\nStddev = %5.6f\n",$mean, $stddev;
# Print histogram
for $i (0..$#histogram) {
   if ($norm) {
   		$histogram[$i] = ($histogram[$i] / ($#data+1) ); 
	}
	#printf: put g instead of f will align float to the right
   printf ("%5.3f\t%12.7f\n",$bins[$i], $histogram[$i]);
}

# --------------------------------------------------------------

sub round {
	return int($_[0] + .5 * ($_[0] <=> 0));
}
