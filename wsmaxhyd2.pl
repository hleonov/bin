#!/usr/bin/perl

# a script to  retrieve the hydrphobicity from a 
# file containing a list of sequences


# energy tables (0= GES; 1=KD; 2=Eisenberg; 3= Moti)
%{$T[0]} = (A => 1.6, C => 2.0, D => -9.2, E => -8.2, F => 3.7,
	         G => 1.0, H => -3.0, I => 3.1, K => -8.8, L => 2.8,    
	         M => 3.4, N => -4.8, P => -0.2, Q => -4.1, R => -12.3,
	   	   S => 0.6, T => 1.2, V => 2.6, W => 1.9, Y => -0.7);

#%{$T[0]} = (A => 0, C => -7.6, D => -43.6, E => -8.2, F => 10.1,
#	         G => 1.0, H => -6.2, I => 6.3, K => -48.8, L => 12.4,    
#	         M => 7.4, N => -30.4, P => -25, Q => -9.7, R => -65.9,
#	   	   S => -7.4, T => -9.2, V => 10.6, W => 3.6, Y => 4.1);


$inFile = shift;				
print "Enter Minimal Win size: \n";
$MIN_WIN = <STDIN>;
chomp $MIN_WIN;
print "Enter Maximal Win size: \n";
$MAX_WIN = <STDIN>;
chomp $MAX_WIN;
print "Enter intervals: \n";
$DW = <STDIN>;
print "Enter start position: \n";
$start= <STDIN>;
print "Enter end position: \n";
$end = <STDIN>;
open(OUT,">$inFile.out");
open(IN,"$inFile");
while(<IN>){
  chop;
  if(/^>/){   # it is a new one
    if($seq){ # it is not the first one
		energy($seq);
    }
    $name = $_;
    $seq = "";
  }
  else{
    s/\s//g;
    $seq .=$_ ;
  }
}
energy($seq) ;
close(OUT);
close(IN);

sub energy {
	$maxVal = -200;
	$maxWin = $MIN_WIN;
	$maxPos = -1;
	undef(@s);
   @s = split(//,$_[0]);
	for (my $ws=$MIN_WIN; $ws <=$MAX_WIN; $ws+=$DW) {
		next if ($ws>$#s+1);
 	 	my @energy;
	 	undef(@energy);
    	for my $i ($start .. $end-$ws) {
			next if ($end<$start);
		for my $t (0 .. $#T) {
		  	for $j ( 0 .. $ws-1) { #sum on window size
				$energy[$t][$i] += $T[$t]{$s[$i+$j]};					
			}
			setMax($energy[$t][$i],$i,$ws);
			
		}	 
    }
   }
    #printing results
	 print OUT "$name";
 	 print OUT "\nWin: $maxWin\n";
    print OUT "Position: $maxPos\n";
	 printf OUT "Energy: %10.1f\n", $maxVal;
 	print OUT "---------------------------------------------------------\n";
}

#setMax(newVal,pos,win)
sub setMax {
	if ($maxVal < $_[0]) {
		$maxVal = $_[0];
		$maxPos = $_[1];
		$maxWin = $_[2];
	}
}

















