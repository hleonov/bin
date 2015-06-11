#!/usr/bin/perl
# originally by Peter Tieleman/
# modif by Graham Smith 

# this script removes all lipids with specified atom (usually P) 
# within at a certain distance from the centre
#ATOM      1  C1  POP     1      10.362  47.025  13.902  0.00  0.00

# usage make_hole.pl -f 128popw_eq.pdb -o edited.pdb -cx 2.7 -cy 3.5 -r 2.0 -lipid POP -lipat P8 

# set up a parsable command line 
%argh = @ARGV ; 

# parse 
if (defined($argh{"-f"})){
    $infile = $argh{"-f"} ; 
} else {
    print " input file not set (-f) - using bilayer.pdb\n" ; 
    $infile = "bilayer.pdb" ; 
}
if (defined($argh{"-o"})){
    $outfile = $argh{"-o"} ; 
} else {
    print " output file not set (-o) - using edited.pdb\n" ; 
    $outfile = "edited.pdb" ; 
}
if (defined($argh{"-lipid"})){
    $lipid = $argh{"-lipid"} ; 
} else {
    print " lipid residue not set (-lipid) - using POP \n" ; 
    $lipid = "POP"; 
}
if (defined($argh{"-lipat"})){
    $lipat = $argh{"-lipat"} ; 
} else {
    print " lipid atom not set (-lipat) - using P8 \n" ; 
    $lipat = "P8"; 
}
if (defined($argh{"-r"})){
    $r = $argh{"-r"} ; 
} else {
    print " radius not set (-r) - using 1 nm \n" ; 
    $r = 1.0 ; 
}
$r  *= 10; 

if (defined($argh{"-cx"})){
    $cx = $argh{"-cx"} ; 
    $cx *= 10; 
}
if (defined($argh{"-cy"})){
    $cy = $argh{"-cy"} ; 
    $cy *= 10; 
}

$r2 = $r*$r;

open(IN, "$infile") or die "cant open $infile (-f argument)" ; 
@lines = (<IN>);
$nr = scalar @lines;
close(IN); 

# just get av z 
$split = 0 ; 
$avx = 0 ; 
$avy = 0 ; 
$ct = 0.001 ; 
for ($i=0;$i<$nr;$i++)
{
    $_ = $lines[$i]; 
    next unless /^ATOM/ or /^HETATM/; 
    ($atom,$anr,$name,$res,$resnr,$x,$y,$z,$dum1,$dum2) = split(' ',$lines[$i]);
    next unless $name eq $lipat;

    $split += $z ; 
    $avx += $x ; 
    $avy += $y ; 
    $ct += 1 ; 
}
$split /= $ct ; 
$avx /= $ct ; 
$avy /= $ct ; 
if (not(defined($argh{"-cx"}))){
    $cx = $avx ; 
    $cxn = $cx/10; 
    print " centre x crd not set (-cy) - using average = $cxn nm \n" ; 
}
if (not(defined($argh{"-cy"}))){
    $cy = $avy ; 
    $cyn = $cy/10; 
    print " centre y crd not set (-cy) - using average = $cyn nm \n" ; 
}

# now do it properly 

for ($i=0;$i<$nr;$i++)
{
    $_ = $lines[$i]; 
    next unless /^ATOM/ or /^HETATM/; 
    ($atom,$anr,$name,$res,$resnr,$x,$y,$z,$dum1,$dum2) = split(' ',$lines[$i]);
    next unless $name eq $lipat;

    $dx = $x - $cx; $dy = $y - $cy;
    $dx2 = $dx*$dx; $dy2 = $dy*$dy;
    if ($dx2 < $r2 && $dy2 < $r2 && ($dx2+$dy2 < $r2))
    {
	$keep[$resnr] = 1;
	if ($z < $split) {print STDERR "removing: side1: $resnr\n";}
	else {print STDERR "removing: side2: $resnr\n";}
    }
}

open(OUT, ">$outfile") or die "cant open $outfile (-o argument)" ; 

for ($j=0;$j<$nr;$j++)
{
    $_ = $lines[$j]; 
    next unless /^ATOM/ or /^HETATM/; 

    ($atom,$anr,$name,$res,$resnr,$x,$y,$z,$dum1,$dum2) = split(' ',$lines[$j]);
 #     print $keep[$resnr];
    if ( $keep[$resnr] != 1  || $res ne $lipid  )
    {
	print OUT $lines[$j];
    }
}
close (OUT); 
exit(0); 













