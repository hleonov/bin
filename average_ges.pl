#!/usr/bin/perl

# Engelman = ges scale				
%G = (A => 1.6, C => 2.0, D => -9.2, E => -8.2, F => 3.7,
	   G => 1.0, H => -3.0, I => 3.1, K => -8.8, L => 2.8,    
	   M => 3.4, N => -4.8, P => -0.2, Q => -4.1, R => -12.3,
	   S => 0.6, T => 1.2, V => 2.6, W => 1.9, Y => -0.7);

# new Von Heine scale	   
%VH = (A => -0.11, C => 0.13, D => -3.49, E => -2.68, F => 0.32,
	   G => -0.74, H => -2.06, I => 0.6, K => -2.71, L => 0.55,    
	   M => 0.1, N => -2.05, P => -2.23, Q => -2.36, R => -2.58,
	   S => -0.84, T => -0.52, V => 0.31, W => -0.3, Y => -0.68);
	   	   
				
$in = shift;
open(IN, $in);
$seq = "";
$i=0;
while (<IN>) {
	 chomp;
  if(m/^>(.*)/){   # it is a new one
    if($seq){ # it is not the first one
      if ($seq ne "") {
		@res=energy($seq) ;
		for $k (0 .. $#res) {
			$e[$k]+=$res[$k];
	#		$e2+=$res[1];
		}
	  }
    }
    $name = $_;
    $seq = "";
	 $i++;
  }
  else{
    s/\s//g;
    $seq .=$_ ;
  }
}
if ($seq ne "") {
	@res=energy($seq) ;
	for $k (0 .. $#res) {
		$e[$k]+=$res[$k];
	}
	$i++;
}
#print "Average GES: ".$e[0]/$i."\nAverage VH: ".$e[1]/$i."\n";
print "amino acid Average GES: ".$e[2]/$i."\namino acid Average VH: ".$e[3]/$i."\n";
sub energy {
	my @sequ = split(//,uc($_[0]));
	my $sum1 = 0;
	my $sum2 = 0;
	for my $j (0 .. $#sequ) {
		$sum1 += $G{$sequ[$j]};
		$sum2 += $VH{$sequ[$j]};
	}	
	return ($sum1,$sum2,$sum1/($#sequ+1),$sum2/($#sequ+1));
	#return ($sum1/($#sequ+1),$sum2/($#sequ+1)) ; (return average of amino acid ges)
}

