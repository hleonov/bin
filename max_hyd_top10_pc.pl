#!/usr/bin/perl

# a script to  retrieve the hydrphobicity from a 
# file containing a list of sequences


# energy tables (0= GES; 1=KD; 2=Eisenberg; 3= Moti)
# GES
%{$T[1]} = (A => 1.6, C => 2.0, D => -9.2, E => -8.2, F => 3.7,
	         G => 1.0, H => -3.0, I => 3.1, K => -8.8, L => 2.8,    
	         M => 3.4, N => -4.8, P => -0.2, Q => -4.1, R => -12.3,
	   	   S => 0.6, T => 1.2, V => 2.6, W => 1.9, Y => -0.7);

# new Von Heine scale      
%{$T[2]} = (A => -0.11, C => 0.13, D => -3.49, E => -2.68, F => 0.32,
           G => -0.74, H => -2.06, I => 0.6, K => -2.71, L => 0.55,    
           M => 0.1, N => -2.05, P => -2.23, Q => -2.36, R => -2.58,
           S => -0.84, T => -0.52, V => 0.31, W => -0.3, Y => -0.68);

#Moti scale??
#%{$T[0]} = (A => 0, C => -7.6, D => -43.6, E => -8.2, F => 10.1,
#	         G => 1.0, H => -6.2, I => 6.3, K => -48.8, L => 12.4,    
#	         M => 7.4, N => -30.4, P => -25, Q => -9.7, R => -65.9,
#	   	   S => -7.4, T => -9.2, V => 10.6, W => 3.6, Y => 4.1);

@table_name = ("","GES", "VH");

print "Enter filename of FASTA input:\n";
$inFile = <STDIN>;
chomp $inFile;
print "Enter Minimal Win size:\n";
$MIN_WIN = <STDIN>;
chomp $MIN_WIN;
print "Enter Maximal Win size:\n";
$MAX_WIN = <STDIN>;
chomp $MAX_WIN;
print "Enter desired number of hits:\n";
$Nbest = <STDIN>;
chomp $Nbest;
print "Enter the length cutoff to merge sequences:\n";
$merge_range = <STDIN>;
chomp $merge_range;

$table = 1; #1=GES, 2=VH
$DW =1; 

my %energy;
my $in = $inFile;
$in =~ s/\..{3}$//; 
my $outFile = "$in\_$MIN_WIN\-$MAX_WIN\_top$Nbest\_$table_name[$table]\_m$merge_range.txt";
print "Wrote: $outFile\n";

&read_fasta($inFile,$outFile); #will run "calc_energy" and "print_results"

sub print_results { 
	print OUT "\n*** Table ". ($t+1) .": $table_name[$table] ***\n";
	print OUT "\n$_[0]\n\nWin-size\tPosition\tHydrophobicy (Descending)\n==================================================\n";
	
	undef(%printed);	
	$topi=0;
	foreach $key (sort hash_value_by_descending (keys(%energy))) {
 	 	 
  		 last if ($topi > $Nbest) ;
		 @temp = split(/_/, $key);
		 $flag=0;
		 
		 #Check if we already printed out something similar
		 for ($r=0; $r<=$merge_range; $r++) {
		 	$flag=1 if (defined $printed{$temp[1]+$r});
			$flag=1 if (defined $printed{$temp[1]-$r});
		 }
		 
		 #overlapping segment was not printed yet, so print out this result and store.
		 if ($flag==0) {
			 $topi++;
			 $printed{$temp[1]} = $temp[0];
			 my $sub_seq = substr($_[1], $temp[1], $temp[0]);
			 printf OUT "%d\t\t\t%d\t\t\t%3.1f\t\t\t%s\n",$temp[0], $temp[1], $energy{$key}, $sub_seq; 
#	  	 	print OUT "$temp[0]\t\t\t$temp[1]\t\t\t$energy{$key}\t\t\t$sub_seq\n";
		}
	}
}

sub calc_energy {
	undef(%energy);
	#split sequence into letters
	undef(@s);
    @s = split(//,$_[0]);
	$start = 0;
	$end = $#s;
	
	#iterate on different win size values
	for (my $ws=$MIN_WIN; $ws <=$MAX_WIN; $ws+=$DW) {
		next if ($ws>$#s+1);
	 	
		#iterate on sequence
    	for my $i ($start .. $end-$ws) {
			#next if ($end<$start);
			
				#iterate from 0 to window size and sum on each amino-acid
			  	for $j ( 0 .. $ws-1) { #sum on window size
					$energy{"$ws\_$i"} += $T[$table]{$s[$i+$j]};	
					#$energy[$ws][$i] += $T[$t]{$s[$i+$j]};					
				}
		}	 
    }
   	
}

sub read_fasta {
	open(OUT,">$_[1]");
	open(IN,"$_[0]");
	while(<IN>){
 		chop;
 	 	if(/^>/){   # new sequence
	  		if($seq){ #  not the first one in the file
				print "calc for $name\n";
				&calc_energy($seq);
				&print_results($name, $seq);
    		}
   			$name = $_;
	    	$seq = "";
  		}
  		else{
    		s/\s//g;
    		$seq .=$_ ;
  		}
	}
	calc_energy(uc($seq)) ;
	&print_results($name, $seq);
	close(IN);
	close(OUT);
}

sub hash_value_by_descending {
   $energy{$b} <=> $energy{$a};
}

sub hash_value_by_ascending {
   $energy{$a} <=> $energy{$b};
}












