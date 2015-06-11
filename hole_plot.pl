#!/usr/bin/perl
if($#ARGV != 2){
  print "Usage: hole_plot.pl file time resolution\n"; 
  print "      \"file\" is the name of the hole output file to read\n"; 
  print "      \"time\" is the length of the simulation in ns\n"; 
  print "      \"resolution\" how often to sampel a time point\n"; 
  print "      Other options are in the code its self (e.g. number of ticks, initial time etc.)\n"; 
  exit;
}

$file = $ARGV[0];
$time = $ARGV[1];   # the total time in ns
                    # its needs to be calculted
$step = $ARGV[2];   # this determines how many time point are selected out of the entire list
                    # i.e. it will select every nth time point.

# optional parameters
$number_of_ticks_time = 4; # this is the number of ticks to plot on the time axis
$number_of_ticks_depth = 4; # this is the number of ticks to plot on the depth axis
$number_of_ticks_pore = 4; # this is the number of ticks to plot on the pore axis
$initial_time = 3;          # this is the start of the time
$LabelAxes = 1; # if it is anything else then the axes are not labled
# ranges for the plots
$LowerPoreRange  = 0; # fixed here    
$UpperPoreRange  = 8; # fixed here
$LowerTimeRange  = 0; # changed later on  
$UpperTimeRange  = 0; # changed later on   
$LowerDepthRange = 0; # changed later on   
$UpperDepthRange = 0; # changed later on  


# initialzation
$i = 0;
$j = 0;
$radius_min = 1000;
$radius_max = 0;
$min = 10000000000000;
open(IN,$file) || die "File \"$file\" cannot be found or opened\n";
while(<IN>){
  # cetch the resolution of the hole run
  if(/sample (.{10}) .+! distance between planes/){
    $resolution = $1;
    $resolution =~ s/\s//g; # remove white spaces
  }
  if(/.{51}\(sampled\)/ || /.{49}\(mid-point\)/){
    @fields = split;
    $position[$i][$j] = $fields[0];
    $radius[$i][$j]   = $fields[1];
    $j++;
  }
  if(/ cenxyz.cvec      radius  cen_.+/){ # a new structure
    if($min > $j && $i > 0){# Remember the smallest list but not for the first one
      $min = $j;
      $start = $position[$i][0];
      $end = $position[$i][$j];
      $length = $start - $end;
    }
    $i++;
    $j = 0;
  }
}
if($min > $j && $i > 0){# Remember the smallest list. Perhpahs its the last one
  $min = $j;
  $start = $position[$i][0];
  $end = $position[$i][$j];
  $length = $start - $end;
}


# all the data is in the arrays
print "labelsize = 12;\n";
print "data = {\n";
# run thru all the lines/strcutures
for($k = 1; $k <= $i; $k += $step){
  print "{";
  # run thru all the points
  for($l = 0; $l < $min - 1; $l++){
    print "$radius[$k][$l],";
    # figure out the maximu and minimu radii
    if($radius_min > $radius[$k][$l]){
      $radius_min = $radius[$k][$l];
    }
    if($radius_max < $radius[$k][$l]){
      $radius_max = $radius[$k][$l];
    }
  }
  print "$radius[$k][$min - 1]"; # the last one
  if($k == $i){
    print "}\n";
  }
  else{
    print "},\n";
  }
}
print "};\n";

# calcualte the depth of the smallest line
# the number of points is $min
# the rosolution of hole is 

$depth = $min * $resolution;

$LowerTimeRange  = 0;
$UpperTimeRange  = $i / $step;
$LowerDepthRange = 0;
$UpperDepthRange = $min;


print "LowerPoreRange  = $LowerPoreRange;\n";
print "UpperPoreRange  = $UpperPoreRange;\n";
print "LowerTimeRange  = $LowerTimeRange;\n";
print "UpperTimeRange  = $UpperTimeRange;\n";
print "LowerDepthRange = $LowerDepthRange;\n";
print "UpperDepthRange = $UpperDepthRange;\n\n";

print "LabelAxes = $LabelAxes; (* If this is anything else then 1 then the axes will not be labled *)\n\n";

print "If[LabelAxes == 1, PoreLabel  = StyleForm[\"Pore/\\[Angstrom]\",FontSize->labelsize], PoreLabel  = \"\"];\n";
print "If[LabelAxes == 1, TimeLabel  = StyleForm[\"Time/ns\",FontSize->labelsize], TimeLabel  = \"\"];\n";
print "If[LabelAxes == 1, DepthLabel  = StyleForm[\"Depth/\\[Angstrom]\",FontSize->labelsize], DepthLabel  = \"\"];\n";



# this does the colored contour plot
print "ListContourPlot[data,ColorFunction->(Hue[#/1.25,1,1]&)\n";
print ",PlotRange -> {{LowerDepthRange,UpperDepthRange},{LowerTimeRange,UpperTimeRange},{LowerPoreRange,UpperPoreRange}}\n";  
print ",FrameLabel->{DepthLabel,TimeLabel}\n";
#print ",FrameLabel->{StyleForm[\"Depth/\\[Angstrom]\",FontSize->labelsize],StyleForm[\"Time/ns\",FontSize->labelsize]}\n";
print ",RotateLabel -> False\n";
print  ",FrameTicks -> {{";
# now the ticks
# printing x ticks (Depth)
# $min is the distances
print  "{1,\"0\"},"; # this is the first tick
for($n = 1; $n < $number_of_ticks_depth; $n++){
  $tick = ($n * $depth)/$number_of_ticks_depth;
  $a = ($n * $min)/($number_of_ticks_depth);
  print  "{$a,\"$tick\"},";
}
print  "{$min,\"$depth\"}},\n{";# this is the last tick
# done printing the x ticks

# printing y ticks (time)
print  "{$step,\"$initial_time\"},"; # this is the first tick
for($n = 1; $n < $number_of_ticks_time; $n++){
  $tick = ($n * $time)/$number_of_ticks_time;
  $a = ($n * $i)/($number_of_ticks_time * $step);
  $tick += $initial_time;
  print  "{$a,\"$tick\"},";
}
# this is for the last tick
$a = $i/$step;
$total_time = $initial_time + $time;
print  "{$a,\"$total_time\"}";
# done printing the y ticks
print  "},{},{}}\n]\;\n\n\n";
# done printing the ticks


# this does the colored 3D plot
print "ListPlot3D[data,ColorFunction->(Hue[#/1.25,1,1]&),Mesh->False";
print ",PlotRange -> {{LowerDepthRange,UpperDepthRange},{LowerTimeRange,UpperTimeRange},{LowerPoreRange,UpperPoreRange}}\n";  
print ",AxesLabel->{DepthLabel,TimeLabel,PoreLabel}\n";
#print ",AxesLabel->{StyleForm[\"Depth/\\[Angstrom]\",FontSize->labelsize],StyleForm[\"Time/ns\",FontSize->labelsize],StyleForm[\"Pore/\\[Angstrom]\",FontSize->labelsize]}\n";
print  ",Ticks -> {{";
# now the ticks
# printing x ticks (Depth)
# $min is the distances
print  "{1,\"0\"},"; # this is the first tick
for($n = 1; $n < $number_of_ticks_depth; $n++){
  $tick = ($n * $depth)/$number_of_ticks_depth;
  $a = ($n * $min)/($number_of_ticks_depth);
  print  "{$a,\"$tick\"},";
}
print  "{$min,\"$depth\"}},\n{";# this is the last tick
# done printing the x ticks

# printing y ticks (time)
print  "{$step,\"$initial_time\"},"; # this is the first tick
for($n = 1; $n < $number_of_ticks_time; $n++){
  $tick = ($n * $time)/$number_of_ticks_time;
  $a = ($n * $i)/($number_of_ticks_time * $step);
  $tick += $initial_time;
  print  "{$a,\"$tick\"},";
}
# this is for the last tick
$a = $i/$step;
$total_time = $initial_time + $time;
print  "{$a,\"$total_time\"}},";
# done printing the y ticks
# now the Z ticks
print  "{{$radius_min,\"$radius_min\"},"; # this is the first tick
# this is for the last tick
print  "{$radius_max,\"$radius_max\"}}";
print  "}\n]\;\n\n\n";
# done printing the ticks




# this does the BW contour plot
print "ListContourPlot[data,ColorFunction->(GrayLevel[#]&)\n";
print ",PlotRange -> {{LowerDepthRange,UpperDepthRange},{LowerTimeRange,UpperTimeRange},{LowerPoreRange,UpperPoreRange}}\n";  
print ",FrameLabel->{DepthLabel,TimeLabel}\n";
#print ",FrameLabel->{StyleForm[\"Depth/\\[Angstrom]\",FontSize->labelsize],StyleForm[\"Time/ns\",FontSize->labelsize]}\n";
print ",RotateLabel -> False\n";
print  ",FrameTicks -> {{";
# now the ticks
# printing x ticks (Depth)
# $min is the distances
print  "{1,\"0\"},"; # this is the first tick
for($n = 1; $n < $number_of_ticks_depth; $n++){
  $tick = ($n * $depth)/$number_of_ticks_depth;
  $a = ($n * $min)/($number_of_ticks_depth);
  print  "{$a,\"$tick\"},";
}
print  "{$min,\"$depth\"}},\n{";# this is the last tick
# done printing the x ticks

# printing y ticks (time)
print  "{$step,\"$initial_time\"},"; # this is the first tick
for($n = 1; $n < $number_of_ticks_time; $n++){
  $tick = ($n * $time)/$number_of_ticks_time;
  $a = ($n * $i)/($number_of_ticks_time * $step);
  $tick += $initial_time;
  print  "{$a,\"$tick\"},";
}
# this is for the last tick
$a = $i/$step;
$total_time = $initial_time + $time;
print  "{$a,\"$total_time\"}";
# done printing the y ticks
print  "},{},{}}\n]\;\n\n\n";
# done printing the ticks



# this does the BW 3D plot
print "ListPlot3D[data,ColorFunction->(GrayLevel[#]&),Mesh->False";
print ",PlotRange -> {{LowerDepthRange,UpperDepthRange},{LowerTimeRange,UpperTimeRange},{LowerPoreRange,UpperPoreRange}}\n";  
print ",AxesLabel->{DepthLabel,TimeLabel,PoreLabel}\n";
#print ",AxesLabel->{StyleForm[\"Depth/\\[Angstrom]\",FontSize->labelsize],StyleForm[\"Time/ns\",FontSize->labelsize],StyleForm[\"Pore/\\[Angstrom]\",FontSize->labelsize]}\n";
print  ",Ticks -> {{";
# now the ticks
# printing x ticks (Depth)
# $min is the distances
print  "{1,\"0\"},"; # this is the first tick
for($n = 1; $n < $number_of_ticks_depth; $n++){
  $tick = ($n * $depth)/$number_of_ticks_depth;
  $a = ($n * $min)/($number_of_ticks_depth);
  print  "{$a,\"$tick\"},";
}
print  "{$min,\"$depth\"}},\n{";# this is the last tick
# done printing the x ticks

# printing y ticks (time)
print  "{$step,\"$initial_time\"},"; # this is the first tick
for($n = 1; $n < $number_of_ticks_time; $n++){
  $tick = ($n * $time)/$number_of_ticks_time;
  $a = ($n * $i)/($number_of_ticks_time * $step);
  $tick += $initial_time;
  print  "{$a,\"$tick\"},";
}
# this is for the last tick
$a = $i/$step;
$total_time = $initial_time + $time;
print  "{$a,\"$total_time\"}},";
# done printing the y ticks
# now the Z ticks
print  "{{$radius_min,\"$radius_min\"},"; # this is the first tick
# this is for the last tick
print  "{$radius_max,\"$radius_max\"}}";
print  "}\n]\;\n";




#print "(*  Depth = $depth  *)\n";
















