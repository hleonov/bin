#!/usr/bin/perl -w

# join_trj.pl
# unite gradual PR files continuously. 
# Usage example: perl ../join_trj.pl 1000 10 20 system_from_PR
# Note: Must change the function in "compute_nsteps" to match

$max = shift || die "Usage: perl $0 <max> <min> <df> <dt> <base-file> <out-base>\n";
$min = shift || die "Usage: perl $0 <max> <min> <df> <dt> <base-file> <out-base>\n";
$df  = shift || die "Usage: perl $0 <max> <min> <df> <dt> <base-file> <out-base>\n";
$dt  = shift || die "Usage: perl $0 <max> <min< <df> <dt> <base-file> <out-base>\n";
$base = shift || die "Usage: perl $0 <max> <min> <df> <dt> <base-file> <out-base>\n";
$out = shift || die "Usage: perl $0 <max> <min> <df> <dt> <base-file> <out-base>\n";
 
$f=$max;
$xtc_str = "";
$trr_str = "";
$edr_str = "";
$time_str = "printf \"";
$time=0;

while ($f>=$min) {
	$time_str .= "$time\\n";

	$xtc_str .= "$base$f.xtc ";
	$trr_str .= "$base$f.trr ";
	$edr_str .= "$base$f.edr ";

	if ($dt =~ m/f/) {	#linear function of time
		$time += compute_time($f);
	} else {	#constant time
		$time += $dt;
	}
	$f = $f - $df;
}

$ene_str = $time_str . "\" | eneconv -settime ";
$time_str .= "\" | trjcat -settime ";

$xtc_com = $time_str . "-o $out.xtc -f " . $xtc_str . "\n";
$trr_com = $time_str . "-o $out.trr -f " . $trr_str . "\n";
$edr_com = $ene_str . "-o $out.edr -f " . $edr_str . "\n";

sub compute_time {
	return (-40*$_[0] + 50000)*0.002; #return time in ps
		
}

unless (system ($xtc_com) == 0) {
	print "Error in creating $out.xtc\n";
}
unless (system ($trr_com) == 0) {
	print "Error in creating $out.trr\n";
}
unless (system ($edr_com) == 0) {
	print "Error in creating $out.edr\n";
}

