#!/usr/bin/perl -w

$K = shift || die "perl $0 <K>\n";
$min = -14.5;
$max = 14.5;
$half = 1;
`mv meta.dat meta.dat.backup`;
for ($i=$min; $i<=$max; $i+=0.5) {
	if ($half) {
		$file = "center_$i"."_pull_k$K";
		$half = 0;
	} else {
		$file = "center_$i".".0_pull_k$K";
		$half = 1;
	}
	#print $file."\n";
	`perl ../trans_pdo.pl $file.pdo $file.ppa`;

}
