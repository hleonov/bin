#!/usr/bin/perl -w

# Engelman = ges scale				
%T = (A => 1.6, C => 2.0, D => -9.2, E => -8.2, F => 3.7,
	   G => 1.0, H => -3.0, I => 3.1, K => -8.8, L => 2.8,    
	   M => 3.4, N => -4.8, P => -0.2, Q => -4.1, R => -12.3,
	   S => 0.6, T => 1.2, V => 2.6, W => 1.9, Y => -0.7);
				
$in = shift;
open(IN, $in);
while (<IN>) {
	chomp;
	if (m/(\w)\s(.*)/) {
		$aa  = $1;
		$p = $2;
		$e += $p*$T{$aa};
	}
}
close(IN);
print "energy: $e\n";
