#!/usr/bin/perl -w

$profile = shift || die "Usage: perl $0 <profile> \n";
open(FREE, $profile) || die "Error: Cannot open $profile\n $!\n";
while (<FREE>) {
	if (m/^(\S+)\s+(\S+)/) {
		($first, $second) = ($1, $2); 
#		print "$first\t$second\n";
		my @first = split(/e/, $first);
		my @second = split(/e/, $second);
		my $z = $first[0]*(10*($first[1]+1));	#switch to Angstrom and get rid of exponential writing
		my $free = $second[0]*(10**($second[1]));	
		print "$z \t $free\n";
		#$rc{$z} = $free;
	}
}
close (FREE);
exit;
