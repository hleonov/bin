#!/usr/bin/perl -w

$base = shift || die "perl $0 <base-name> <movie-name> <last-frame>\n";
$out  = shift || die "perl $0 <base-name> <movie-name> <last-frame>\n";
$last = shift || die "perl $0 <base-name> <movie-name> <last-frame>\n";

$parfile = "movie.par";
open(PAR, ">$parfile") || die "Cannot open parameter file\n";
#print PAR "PATTERN    IBBPBBPBBPBBPBB\n";
print PAR "PATTERN    I\n";
print PAR "FORCE_ENCODE_LAST_FRAME\n";		# force anim loopable
print PAR "OUTPUT	  $out.mpg\n";
print PAR "INPUT_DIR  .\n";
print PAR "INPUT\n";

printf PAR "%s.*.ppm \[%.4d-%.4d\]\n", $base, "0000", $last; 
#print PAR "$base.*.ppm \[0000-$last\]"

print PAR "END_INPUT\n";
print PAR "BASE_FILE_FORMAT PPM\n";
print PAR "INPUT_CONVERT *\n";
print PAR "GOP_SIZE 1\n";
print PAR "SLICES_PER_FRAME 1\n";
print PAR "PIXEL FULL\n";

#print PAR "PIXEL HALF\n";
print PAR "RANGE 32\n";
print PAR "PSEARCH_ALG LOGARITHMIC\n";
print PAR "BSEARCH_ALG CROSS2\n";
print PAR "IQSCALE 1\n";
print PAR "PQSCALE 1\n"; #5
print PAR "BQSCALE 1\n";	#12
print PAR "REFERENCE_FRAME DECODED";
close(PAR);


`ppmtompeg movie.par`;
