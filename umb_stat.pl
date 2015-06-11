#!/usr/bin/perl -w

# This script outputs some statistics about the bias from the reference
# point during the umbrella sampling results.
# for individual windows, and overall. 

#init to opposite very small numbers
$all_max = -100;
$all_average = 0;
undef(@all_list);
foreach my $file (`ls -1 *.pdo`) {
	chomp $file;
	$one_max = -100;
	$one_avg = 0;
	undef(@one_list);
	open (PDO,"$file") || die "Cannot open $file\n";
	while (<PDO>) {
		if (m/(\d+.\d+)\s+(-?\d+.\d+)/) {
			$shift = abs($2);
			$one_avg += $shift;
			$all_average += $shift;
			push (@all_list, $shift);
			push (@one_list, $shift); 
			if ($one_max < $shift) {
				$one_max = $shift;
			}
		}
	}
	close(PDO);
	$one_max = $one_max*10;
	$one_avg = $one_avg*10 / ($#one_list+1);
	$one_std = 10*std(\@one_list, $one_avg/10);
	printf "File: %s\t Max: %2.3f\t Average: %2.3f\t Std: %2.3f\n", $file, $one_max, $one_avg, $one_std;
	if ($all_max < $one_max) {
		$all_max = $one_max;
	}
}

#multiply all by 10 so units are in Angstrom
$all_max = $all_max;
$all_average = $all_average*10 / ($#all_list+1);
$std_dev = 10*std(\@all_list, $all_average/10);
printf "\n\nAll max: %2.3f\nAll Average: %2.3f\nAll Std: %2.3f\n", $all_max, $all_average, $std_dev;

sub std {
	my @arr = @{$_[0]};
	my $avg = $_[1];
	my $std = 0;
	foreach $elem (@arr) {
		$std += ($elem-$avg)**2;
	}
	$std /= ($#arr+1);
	return sqrt($std);
}
