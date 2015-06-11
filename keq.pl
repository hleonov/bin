# in kJ
$r = 8.314;
$t = 310;
$k = exp(-($ARGV[0]*1000)/($r*$t));
print "$ARGV[0] => $k\n";

