#!/usr/bin/perl -w

$infile = shift || die "Usage: perl $0 <itp>\n";

open(IN,$infile) || die "Cannot open input file $infile\n $!\n";
$out=$infile;
$out =~ s/\.itp/_topB/;

open(ITP, ">$out.itp") 	 || die "Cannot open output file $out\n $!\n";
open(FF, ">ff$out.itp")  || die "Cannot open output file ff$out\n $!\n";
print FF "[ atomtypes ]\n";

$flag=0;
undef(%at_names);

while (<IN>) {
	if (m/^\;/) {
		print ITP $_;			
	} elsif (m/\[\s+(\S+)\s+\]/) {
		print ITP $_;
		my $directive = $1;
		$flag = ($directive eq 'atoms') ? 1 : 0;
	} else {
		if ($flag) {
			chomp;
		    #1     O2		    1    atp    O1G      1   -0.95260  15.999400
			if (($type, $mass) = m/^\s*\S+\s+(\S+).*\S+\s+(\S+)\s*$/) {
				print "$type\t$mass\n";
				print ITP $_."\tDUM_$type\t0.0\t$mass\n";
				if (!defined $at_names{"DUM_$type"}) {
					print FF "DUM_$type\t 0.0\t0.0\tA\t0.0\t0.0\n"
				}
				$at_names{"DUM_$type"} = 1;
			}
		} else {
			print ITP $_;
		}
	}
}

close(FF);
close(ITP);
close(IN);
