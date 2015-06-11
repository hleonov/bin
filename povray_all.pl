#!/usr/bin/env perl
# supply a size

# a script for pov ray

$size = $ARGV[0];
@files = split(/\s+/,`ls *.pov`);

foreach $file (@files){
  print "$file\n";
  $file =~/(.+)\.pov/;
  $prefix = $1;
  $out = $prefix.".ppm";
  # for transpereceny
  #system ("/sw/bin/povray +H$size +W$size -I$file -O$out +D +X +A +FT +UA");
  system ("/opt/local/bin/povray +H$size +W$size -I$file -O$out +D +X +A +FP Display=off");
  #system ("/sw/bin/povray +H$size +W$size -I$file -O$out -D +X +A +FT");
}

