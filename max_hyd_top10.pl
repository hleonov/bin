#!/usr/bin/perl

# a script to  retrieve the hydrphobicity from a 
# file containing a list of sequences


# energy tables (1= GES; 2=KD; 3=Whimley-White) - deleted tables 2-3.
%{$T[1]} = (A => 1.6, C => 2.0, D => -9.2, E => -8.2, F => 3.7,
	         G => 1.0, H => -3.0, I => 3.1, K => -8.8, L => 2.8,    
	         M => 3.4, N => -4.8, P => -0.2, Q => -4.1, R => -12.3,
	   	   S => 0.6, T => 1.2, V => 2.6, W => 1.9, Y => -0.7);


@table_name = ("","GES", "KD","WW");
$inFile = shift || die "Usage: perl $0 <infile> <minwin> <maxwin> <#best> <merge_range>\n";
$MIN_WIN = shift || die "Usage: perl $0 <infile> <minwin> <maxwin> <#best> <merge_range>\n";
$MAX_WIN = shift || die "Usage: perl $0 <infile> <minwin> <maxwin> <#best> <merge_range>\n";
$Nbest = shift || die "Usage: perl $0 <infile> <minwin> <maxwin> <#best> <merge_range>\n";
$merge_range = shift || die "Usage: perl $0 <infile> <minwin> <maxwin> <#best> <merge_range>\n";

$table = 1; 
$DW =1; 

my %energy;
my $in = $inFile;
$in =~ s/\..{3}$//;
my $outFile = "$in\_$MIN_WIN\-$MAX_WIN\_top$Nbest\_$table_name[$table]\_m$merge_range.txt";
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
				print "Computing hydrophobicity for $name\n";
				&calc_energy(uc($seq));
				&special_chars(uc($seq));
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
	#for last
	print "Computing hydrophobicity for $name\n";
	calc_energy(uc($seq));
	&special_chars(uc($seq));
	&print_results($name, $seq);
	close(IN);
	close(OUT);
}

sub special_chars {
	my $sequ = $_[0];
	if ($sequ =~ m/[BJOUXZ]/ || $sequ =~ m/[^A-Z]/) {
		print "Warning: special non- amino acids characters were found in the sequence\n";
		while ($sequ =~ m/[BJOUXZ]/g) {
			print pos($sequ)."(".substr($sequ,pos($sequ)-1,1).")\n";
		}
		while ($sequ =~ m/[^A-Z]/g) {
			print pos($sequ)."(".substr($sequ,pos($sequ)-1,1).")\n";
		}
	}
	
}

sub hash_value_by_descending {
   $energy{$b} <=> $energy{$a};
}

sub hash_value_by_ascending {
   $energy{$a} <=> $energy{$b};
}












