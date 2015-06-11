#!/usr/bin/perl

#Below is a subroutine that produces ListDensityPlots from a 2D associative
#array of data.

#The inputs to the subroutine are:

#1, (Reference to) Associative array of data - gives the point values
#2, (Reference to) A list of keys for the X-axis
#3, (Reference to) A list of keys for the Y-axis
#4, Output file name

#Note that the axis keys function both to extract the array value and to
#make the labels for the figure.

#The colour is set by the (Hue[h,s,b]&) entry. Within this field "#" refers
#to the normalised value of the function (or array in this case), it will
#be between 0 and 1.

#Thus, in the below examples as "#" goes from 0 to 1:
#The hue varies from 0.5 to 1
#The saturation from 0   to 1
#The brightness from 1   to 0.6


#-----------------------------------------------------------------------------
%array_of_data = undef;
@x_keys = undef;
@y_keys = undef;
$in = shift;
$outfile = $in.'.nb';
$dataName = $in;
$dataName =~ s/_//g;   $dataName =~ s/stat\.//g;
&readData($in);
#&plot(\%array_of_data,\@x_keys,\@y_keys,$outfile);
&plot(\%array_of_data,$outfile);
sub plot {
  my ($hash,$outfile)=@_;
  my $str ="\{";
  open (Q, ">$outfile");
  print Q "$dataName=";
  foreach $x (sort keys(%array_of_data)) {
    $str .="\{";
    foreach $y (sort keys(%{$array_of_data{$x}})) {
  
#  @xx=@$lx;
#  @yy=@$ly;
 # for $x (@xx) {
 #     next if !defined($x);

#    for $y (@yy) {
#	next if !defined($y);
      my $s = $array_of_data{$x}{$y};
	  print "hash $x $y = $s\n";
      if (! $s) {
        $s="0.0";
      }
      $str .= "$s,";
    }
    $str =~ s/,$/},\n/;
  }
  $str =~ s/,\n$/};\n/;
  print Q "$str\nListDensityPlot[1-$dataName,Mesh -> False,ColorFunctionScaling -> False, FrameTicks -> {{";
#ColorFunction -> (Hue[0.5*#+0.5,#,1-#*0.4]&),
  $n=0;
  for $y (@yy) {
      next if !defined($y);
   my $m=$n+0.5;
   print Q "{$m,\"$y\"}";
    $n++;
    unless ($y eq $yy[$#yy]) {
      print Q ",";
    }
  }
  print Q "},{";
  $n=0;
  for $x (@xx) {
      next if !defined($x);
    my $m=$n+0.5;
    $name=$x;
    $name =~ /([a-z])/;
    my $l=$1;
    $l =~ tr/[a-z]/[A-Z]/;
    $name =~ s/([a-z])_/$l\./;
    print Q "{$m,\"$name\"}";
    $n++;
    unless ($x eq $xx[$#xx]) {
      print Q ",";
    }
  }
  print Q "}}]\n\n";
  close Q;
}

sub readData {
    open(DATA,$_[0]);
    my %xkeys = undef; my %ykeys = undef;
    while (<DATA>) {
	if (m/(\S+)\s(\S+)\s(\S+)/) {
	    (my $x, my $y, my $val) = ($1,$3,$2);
		chomp $val;
	    $xkeys{$x}++;
	    $ykeys{$y}++;
	   # push (@x_keys, $x) if ($xkeys{$x} == 1);
	   # push (@y_keys, $y) if ($ykeys{$y} == 1);
	    $array_of_data{$x}{$y} = $val;
		#print "hash $x $y = $array_of_data{$x}{$y}\n";
	}
    }
	
    close(DATA);
}
