#!/usr/bin/perl

# Input: Fasta file of protein sequences.
# Output: a 2D histogram. (file: infile.2D_prot.histo)
#			The X-axis is the energy threshold. 
#			The Y-axis is the percentage of proteins which have a region 
#			(possible size MIN_WIN..MAX_WIN with a GES value > thresh

# energy tables (0= Engelman; 1=KD; 2=Eisenberg; 3= Moti)
%{$T[0]} = (A => 1.6, C => 2.0, D => -9.2, E => -8.2, F => 3.7,
	         G => 1.0, H => -3.0, I => 3.1, K => -8.8, L => 2.8,    
	         M => 3.4, N => -4.8, P => -0.2, Q => -4.1, R => -12.3,
	   	   S => 0.6, T => 1.2, V => 2.6, W => 1.9, Y => -0.7);

$PERC = 0.75;
$MIN_THRESH = 10;
$MAX_THRESH = 40;
$DW = 1;
$MIN_WIN=10;
$MAX_WIN=30;

$inFile = shift or die "Usage: perl $0 fasta-seq-files\n";				

undef(@proteins);   #$proteins[i] = X if there are X proteins s.t there exists a hyd region with GES>i

open(IN,"$inFile") or die "Cannot find $inFile: !$ \n";
#$total = 3458235;
$total = `grep '^>' $inFile | wc -l`;
#$total =~ s/\s//g;
while(<IN>){
  chop;
  if(/^>/){   # it is a new one
    if($seq){ # it is not the first one
      
		energy($seq);
    }
    $name = substr($_,4,4);	
    $seq = "";
  }
  else{
    s/\s//g;
    $seq .=$_ ;
  }
}
energy($seq) ;
open(PR,">$inFile.2D_prot.histo");

for my $thresh ($MIN_THRESH .. $MAX_THRESH) {
	
			if (!defined $proteins[$thresh]) {
				$proteins[$thresh] = 0;
			}
			$proteins[$thresh] = 100*$proteins[$thresh] /$total;
			print PR "$thresh $proteins[$thresh]\n";
		
}
close(PR);
close(IN);

sub energy {
	my @s = split(//,$_[0]);
	my $best_hyd = -500;		#just a really low value	
	for (my $ws=$MIN_WIN; $ws <=$MAX_WIN; $ws+=$DW) {
	
		my @energy;
		undef(@energy);

		#computes energy foreach starting position
		for $j ( 0 .. $ws-1) {							
			$energy[0] += $T[0]{$s[$j]};					
		}
		$best_hyd = $energy[0] if ($energy[0] > $best_hyd);
#		print "energy[0] = $energy[0]\n";
	   for my $i (1 .. $#s-$ws+1) {	#start pos
			$energy[$i] = $energy[$i-1]-$T[0]{$s[$i-1]}+$T[0]{$s[$i+$ws-1]};
			$best_hyd = $energy[$i] if ($energy[$i] > $best_hyd);
#			print "energy[$i] = $energy[$i]\n";		
		}
	} #win size loop
	
	for my $thresh ($MIN_THRESH .. $MAX_THRESH) {
		#count this protein if its best hydrophobic segment passes $thresh
		if ($best_hyd >= $thresh) {
			$proteins[$thresh]++;
		}  	
	}#thresh loop	

}
