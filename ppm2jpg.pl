#!/usr/bin/env perl 
@files = split(/\s+/,`ls *.ppm`);
foreach $file (@files){
  $file =~ /(.+)\.ppm/;
  $start = $1;
  $jpg = $start.".jpg";
  print "$file => $jpg\n";
  system("/opt/local/bin/ppmtojpeg $file > $jpg");
}
