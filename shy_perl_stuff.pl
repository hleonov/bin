#!/usr/bin/env perl
# a bunch of subroutines that contain only local vairables (i.e. my is used)
# the subroutined should not depened on anything exepct themselves
# usage requires only putting the folling line in your script:
# require "/Users/shy/Bin/shy_perl_stuff.pl";

use List::Util qw[min max];
use File::Glob ':glob'; # this overrides the glob to take accoutn of white
                        # spaces
#use strict;

sub help {
  # Takes 2 arugments: 
  # 1 is a message
  # 2 is a flag
  # if an flag is supplied, print the message and exit
  if ($_[1]) {
    print $_[0];
    exit;
  }
}

sub in_file {
  # check to see if a files exists.
  # if not pester untill you get one that does
  my $file = $_[0];
  while(!(-e $file)){
    print STDERR "File \"$file\" doesn't exist. Please supply new file name\n";
    chop($file = <STDIN>);
  }
  return $file;
}
sub out_file {
  # retruns a file name but checks to see if a files exists.
  # if yes ask if you can overwrite, or ask for a new file name and recheck.
  # if a 2nd element is supplied (using the -rf switch) then it doesn't check.
  my $file = $_[0];
  my $rf   = $_[1];
  if ($rf){
    return $file;
  }
  else{
    while(-e $file){
      print STDERR "File \"$file\" exists, overwrite? [y]\n";
      chop($response = <STDIN>);
      if($response eq "y"){
    	last;
      }
      else{
    	print "Please suplly output file name.\n";
    	chop($file = <STDIN>);
      }
    }
    return $file;
  }
}
sub line_message {
  # print a message whenever the counter is devisible by the periodicity evenly
  my $counter = $_[0];
  my $periodicity = $_[1];
  if($counter/$periodicity == int($counter/$periodicity)){
    print "Analysing line $counter\n";
  }
}
sub trim {
  # returns a string without any white charcaters
  my $out = $_[0];
  $out =~ s/\s//g;
  return $out;
}
sub round {
  # returns a rounded number
  return sprintf("%.0f", $_[0]);
}
sub distance {
  # distanec between two points
  my $xa = $_[0];
  my $ya = $_[1];
  my $za = $_[2];
  my $xb = $_[3];
  my $yb = $_[4];
  my $zb = $_[5];
  return sqrt(($xa-$xb)*($xa-$xb)+($ya-$yb)*($ya-$yb)+($za-$zb)*($za-$zb))
}
sub return_color_white_to_black {
  # this retruns the rgb value from 0 to $max_color (defualt is 255)
  # when given the following (in this order):
  # min max value
  # this a gray scale from white to black
  my $max_color = 255;
  my $min   = $_[0];
  my $max   = $_[1];
  my $value = $_[2];
  my $c;
  if($value < $min){
    print STDERR "Value $value is smaller than minimum value $min\n";
  }
  elsif($value > $max){
    print STDERR "Value $value is bigger than maximum value $max\n";
  }
  else{
    $c = ($max_color * ($value - $min) / ($max - $min));
    $c = $max_color - $c;
    $c = sprintf("%.0f", $c);
    return ($c, $c, $c);
  }
}
sub return_color_red_to_blue {
  # this retruns the rgb value from 0 to $max_color (defualt is 255)
  # when given the following (in this order):
  # min max value
  # the scale here is from red to blue without any green.
  # the numbers are intergers
  my $max_color = 255;
  my $min   = $_[0];
  my $max   = $_[1];
  my $value = $_[2];
  my $r;
  my $g;
  my $b;
  if($value < $min){
    print STDERR "Value $value is smaller than minimum value $min\n";
  }
  elsif($value > $max){
    print STDERR "Value $value is bigger than maximum value $max\n";
  }
  else{
    $b = ($max_color * ($value - $min) / ($max - $min));
    # now round it
    $b = sprintf("%.0f", $b);
    $r = $max_color - $b;
    $g = 0;
    return ($r, $g, $b);
  }
}
sub return_color_blue_to_red {
  # this retruns the rgb value from 0 to $max_color (defualt is 255)
  # when given the following (in this order):
  # min max value
  # the scale here is from blue to red without any green.
  # the numbers are intergers
  my $max_color = 255;
  my $min   = $_[0];
  my $max   = $_[1];
  my $value = $_[2];
  my $r;
  my $g;
  my $b;
  if($value < $min){
    print STDERR "Value $value is smaller than minimum value $min\n";
  }
  elsif($value > $max){
    print STDERR "Value $value is bigger than maximum value $max\n";
  }
  else{
    $r = ($max_color * ($value - $min) / ($max - $min));
    # now round it
    $r = sprintf("%.0f", $r);
    $b = $max_color - $r;
    $g = 0;
    return ($r, $g, $b);
  }
}
sub return_color_red_to_green_to_blue {
  # this retruns the rgb value from 0 to $max_color (defualt is 255)
  # when given the following (in this order):
  # min max value
  # the scale here is from red to green to blue.
  # the numbers are intergers
  my $max_color = 255;
  my $min   = $_[0];
  my $max   = $_[1];
  my $value = $_[2];
  my $r;
  my $g;
  my $b;
  my $x;
  if($value < $min){
    print STDERR "Value $value is smaller than minimum value $min\n";
  }
  elsif($value > $max){
    print STDERR "Value $value is bigger than maximum value $max\n";
  }
  else{
    $x = ($max_color * ($value - $min) / ($max - $min));
    if($x <= $max_color / 2){
      # only red and green vary and blue is set to zero
      $g = $x * 2;
      $g = sprintf("%.0f", $g);
      $r = $max_color - $g;
      $b = 0;
    }
    else{
      # only blue and green vary and red is set to zero
      $r = 0;
      $g = ($max_color - $x) * 2;
      $g = sprintf("%.0f", $g);
      $b = $max_color - $g;
    }
    return ($r, $g, $b);
  }
}
sub histogram {
  # this returns a histogram of an array
  # in 1 array that has both x and y
  # the inputs are the number of bins, and then the actual data
  my $bin_number = $_[0];
  my @list       = splice(@_, 1);
  my $min        = &min(@list);
  my $max        = &max(@list);
  my $range      = $max - $min;
  my $bin_size   = $range / $bin_number;
  my @out, @xy, $x, $y;
  foreach my $i (@list){
    my $location = int(($i - $min) / $bin_size);
    $out[$location]++;  
  }
  # now assign int0 an x,y which can then be split
  for (my $j = 0; $j < $bin_number; $j++){
    $x = $min + ($j * $bin_size);
    if($out[$j]){
      $y = $out[$j];
    }
    else {
      $y = 0;
    }
    $xy[$j] = $x." ".$y;
  }
  return @xy;
}
sub recurse($) {
  # this yeilds all files in a directory
  # and in its subdirectories recusrivly.
  # it uses the glob from File::Glob that takes
  # into account white spaces in files names
  my($path)  = @_;
  # append a trailing / if it's not there
  $path .= '/' if($path !~ /\/$/);
  # loop through the files contained in the directory
  for my $eachFile (glob($path.'*')) {
    # if the file is a directory
    if( -d $eachFile) {
      # pass the directory to the routine ( recursion )
      recurse($eachFile,@result);
    } else {
      # print the file
      # print "$eachFile\n";
      # add it to the array
      push(@an_array_that_i_will_never_use_elsewhere, $eachFile);
    }
  }
  return @an_array_that_i_will_never_use_elsewhere;
}












# this is needed to return true.
1;

















