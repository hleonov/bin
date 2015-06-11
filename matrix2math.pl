# make a 2D matrix for mathematica
# from a simple 2D matrix in text
$file = shift || die "perl $0 <infile> [rev?]\n";
$reverse = shift;

open (IN, $file) || die "Error: cannot open $file\n";
$min = 1e38;
$max = -1e38;
$text = "data = {";
while(<IN>){
  $row++;
  # trim leading spaces
  s/^\s+//;
  # trim lagging spaces
  s/\s+$//;
  @elements = split(/\s+/,$_);
  $text .= "{";
  for ($i = 0; $i < $#elements; $i++) {
    if ($elements[$i] > $max){
      $max = $elements[$i];
      $row_max = $row;
      $column_max = $i + 1;
    }
    if ($elements[$i] < $min){
      $min = $elements[$i];
      $row_min = $row;
      $column_min= $i + 1;
    }
    $text .= "$elements[$i],";  
  }
  $text .= "$elements[$#elements]},\n";
}
# to get rid of the annoying last comma
$text = substr($text,0,length($text)-2);
print $text;
print "};\n";
if ($reverse) {
	print "revdata = Map[Reverse, data];\n";
	print "ListPlot3D[revdata, ColorFunction->(Hue[0.9*#]&), Mesh->False, PlotRange->{{0,60},{0,200},{0,15}}];\n";
	#print "ListDensityPlot[revdata, ColorFunction->(Hue[0.9*#]&), Mesh->False, PlotRange->{{0,60},{0,200},{0,15}}];\n";
} else {
	print "ListPlot3D[data, ColorFunction->(Hue[0.9*#]&), Mesh->False, PlotRange->{{0,60},{0,200},{0,15}}];\n";
	#print "ListDensityPlot[data, ColorFunction->(Hue[0.9*#]&), Mesh->False, PlotRange->{{0,60},{0,200},{0,15}}];\n";
	
}
print STDERR "Max of $max at Row $row_max, Column $column_max\n";
print STDERR "Min of $min at Row $row_min, Column $column_min\n";

close(IN);











