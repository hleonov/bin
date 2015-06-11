#!/usr/bin/perl -w 

# Reads a hole output (log) with multiple trajectories to analyze 
# and creates a 3 column dat file <frame> <rad> <z> for gnuplot.

$file = shift || die "Usage: perl $0 <hole log> <out> <min> <max> <dz>\n";
$out = shift || die "Usage: perl $0 <hole log> <out> <min> <max> <dz>\n";
#$min_limit  = shift || die "Usage: perl $0 <hole log> <out> <min> <max> <dz>\n";
#$max_limit = shift || die "Usage: perl $0 <hole log> <out> <min> <max> <dz>\n";
$dz = shift || die "Usage: perl $0 <hole log> <out> <min> <max> <dz>\n";
undef(@rad);
$SM_CUTOFF = 2;
$z_min = 10000;
$z_max = -10000;
$offset = 0;
$endrad=14;
&read_file;
#$minL = getZbin($min_limit, $dz);
#$maxL = getZbin($max_limit, $dz);

&fix_matrix;
&fixed_out;
#$type eq "gnu" ? &gnu_out : &matrix_out;



sub matrix_out {
	open(OUT, ">$out");
	for ($frame=0; $frame<=$#hole; $frame++) {
		for ($z=$z_min; $z<=$z_max; $z++) {
		#print "min: $min[$frame]\t max: $max[$frame]\n";
		#for ($z=$z_min; $z<=$z_max; $z++) {
			#foreach $z (sort {$a<=>$b} keys( %{$hole[$frame]} )) {
			if (defined $hole[$frame][$z]) {
			    $to_print = $hole[$frame][$z] / $num[$frame][$z];
			} 
#			else {
#				$to_print = 7;
#			}	
			
			elsif ($z<$min[$frame]) {
				#print "(Z) $z < $min[$frame] (min) for frame $frame\n";
				my $tz = $min[$frame];
				$to_print = $hole[$frame][$tz] / $num[$frame][$tz];
			} elsif ($z>$max[$frame]) {
				
				#print "(Z) $z > $max[$frame] (max) for frame $frame\n";
				my $tz = $max[$frame];
				$to_print = $hole[$frame][$tz] / $num[$frame][$tz];
			}		
			
	
#			print "frame : $frame\t Z: $z\t cul Radius: $hole[$frame][$z]  Radius:  $to_print\t Number: $num[$frame][$z]\n";
			

			#print OUT "$hole[$frame]{$z}\t";
			printf OUT "%3.3f\t",$to_print;
		}
		print OUT "\n";
	}
	close(OUT);
}

sub fixed_out {
open(OUT, ">$out");
	for ($frame=0; $frame<=$#hole; $frame++) {
		#$tmp = $hole[$frame][$min[$frame]] / $num[$frame][$min[$frame]];
		#print "min: $frame\t $min[$frame]\t $tmp\t $rad[$frame][$min[$frame]] \n";
		for ($z=$z_min; $z<=$z_max; $z++) {
			if (defined $rad[$frame][$z]) {
			    $to_print = $rad[$frame][$z]; 
			} 
			printf OUT "%3.3f\t",$to_print;
		}
		print OUT "\n";
	}
	close(OUT);
}
sub fix_matrix {
	
	for ($frame=0; $frame<=$#hole; $frame++) {
		#print "frame $frame \t min: $min[$frame]\t max: $max[$frame]\n";
		
		# Go over all range of Z and put the last calculated value where missing
		for ($z=$z_min; $z<=$z_max; $z++) {
			if (defined $hole[$frame][$z]) {
			    $rad[$frame][$z] = $hole[$frame][$z] / $num[$frame][$z];
		#		print "frame $frame: defined $z: $rad[$frame][$z] \n";
			}	elsif ($z<$min[$frame]) {
				my $tz = $min[$frame];
				$rad[$frame][$z] = $hole[$frame][$tz] / $num[$frame][$tz];
		#		print "frame $frame: $z<$min[$frame]: $rad[$frame][$z] \n";
			} elsif ($z>$max[$frame]) {
				my $tz = $max[$frame];
				$rad[$frame][$z] = $hole[$frame][$tz] / $num[$frame][$tz];
		#		print "frame $frame: $z>$max[$frame]: $rad[$frame][$z] \n";
			}		 
		}
	
		my $tz1 = $min[$frame];
		if ($tz1 > $z_min) {
			#print "frame: $frame\t tz: $tz1\t$rad[$frame][$tz1]\t$rad[$frame][$tz1+1]\n";
			#my $diff1 = abs($rad[$frame][$tz1] - $rad[$frame][$tz1+1]);
			my $diff1 = abs($endrad- $rad[$frame][$tz1+1]);
			if  ($diff1 > $SM_CUTOFF ) {
				$ds = $diff1/($tz1-$z_min);
				for ($z=$tz1; $z>=$z_min; $z--) {
					$rad[$frame][$z] = $rad[$frame][$z+1] + $ds;
				}
			}
		}
			
		my $tz2 = $max[$frame];
		if ($tz2 < $z_max) {
			#my $diff2 = abs($rad[$frame][$tz2] - $rad[$frame][$tz2-1]);
			my $diff2 = abs($endrad - $rad[$frame][$tz2-1]);
			if  ($diff2 > $SM_CUTOFF ) {
				$ds = $diff2/($z_max-$tz2);
				for ($z=$tz2; $z<=$z_max; $z++) {
					$rad[$frame][$z] = $rad[$frame][$z-1] + $ds;
				}
			}
		}
	}
}


sub getZbin {
	return int($_[0]/$_[1]) + abs($offset);
}
sub read_file {

	my $frame = "";
	my $z = "";
	my $radius = "";
	$z_min = 10000;
	$z_max = -10000;
	open(IN,$file) || die "Cannot open $file\n$!\n";
	while (<IN>) {
		if (m/From file:\s\S+_(\d+)\.pdb/) {
			$frame = $1;
			$max[$frame] =  -10000;
			$min[$frame] = 10000;
		}
		elsif (m/cenxyz.cvec\s+radius\scen.line.di\sinteg/) {
			$flag = 1;
		}
		elsif ($flag && m/^\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/) {
			$z = $1;
			$radius = $2;
			#$zbin = int($z);
			if ($z<0) {
				$offset = $z<$offset ? int($z/$dz) : $offset;
			} 
			$zbin = getZbin($z, $dz);
			
			#update min and max (in discrete bin steps)
			$z_min = $z_min > $zbin ? $zbin : $z_min;
			$z_max = $z_max < $zbin ? $zbin : $z_max;
			
			#update frame min-max (Zbin values)
			$min[$frame] = $min[$frame] > $zbin ? $zbin : $min[$frame];
			$max[$frame] = $max[$frame] < $zbin ? $zbin : $max[$frame];
			
			$hole[$frame][$zbin] += $radius;
			$num[$frame][$zbin]++;
			
		}
		elsif (m/Minimum radius found/) {
			$flag = 0;
		}
	}
	close(IN);
	print "min: $z_min\nmax: $z_max\n";
}

sub gnu_out {
	open(OUT, ">$out");
	for ($frame=0; $frame<=$#hole; $frame++) {
		for ($z=$z_min; $z<=$z_max; $z++) {
#		foreach $z (sort {$a<=>$b} keys( %{$hole[$frame]} )) {
			if (defined $hole[$frame][$z]) {
				$to_print = $hole[$frame][$z] / $num[$frame][$z];
			} else {
				$to_print = 7;
			}
			#print OUT "$frame $hole[$frame]{$z} $z\n";
			printf OUT "%d %3.3f %d\n", $frame, $to_print, $z;
		}
		print OUT "\n";
	}
	close(OUT);
}
