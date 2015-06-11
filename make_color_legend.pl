#!/usr/bin/perl
# this program make a color legend table for
# mathematica, whereby you input the color ranges,
# and the number of objects.
unless($#ARGV == 4){
  print "Usage: make_color_legend.pl low high number_of_color number_of_ticks landscape|protrait\n";
  exit;
}
$start      = $ARGV[0];
$end        = $ARGV[1];
$number     = $ARGV[2];
$ticks      = $ARGV[3];
$landscape  = $ARGV[4];
print "ColorLegend = {\n";
for($i = $start; $i < $end; $i += ($end - $start) / $number){
  print "{$i, $i},\n";
}
print "{$end, $end}\n";
print "}\;\n";
# figure out whether to print out it in landscape or portrait
if($landscape eq "landscape"){
  print "ColorLegend = Transpose[ColorLegend]\;";
}
print "ListDensityPlot[ColorLegend, \n";
print "ColorFunction -> (Hue[\#/1.25,1,1]\&),Mesh->False,\n";
if($landscape eq "landscape"){
  print "AspectRatio -> 0.1,\n";
}
else{
  print "AspectRatio -> 10,\n";
}
#taking care of the frame ticks
$ticks = $ticks - 1;
# if it is in landscape
if($landscape eq "landscape"){
  print "FrameTicks -> {{\n";
  for($i = 0; $i < $ticks; $i++){
    $tmp   = $start + $i * ($end - $start) / $ticks;
    $place = $i * ($number / $ticks);
    print "{$place, \"$tmp\"},";
  }
  $tmp   = $start + $ticks * ($end - $start) /$ticks;
  $place = $ticks * ($number / $ticks);
  print "{$place, \"$tmp\"}";
  print "},{},{},{}}\n";
}
# if it is in protrait
else{
  print "FrameTicks -> {{},{\n";
  for($i = 0; $i < $ticks; $i++){
    $tmp   = $start + $i * ($end - $start) / $ticks;
    $place = $i * ($number / $ticks);
    print "{$place, \"$tmp\"},";
  }
  $tmp   = $start + $ticks * ($end - $start) /$ticks;
  $place = $ticks * ($number / $ticks);
  print "{$place, \"$tmp\"}";
  print "},{},{}}\n";
}

# done with the frameticks
print "]\;\n";
