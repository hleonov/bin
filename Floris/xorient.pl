#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my (@in, @out);
my $ligandsel;
my $veclength = 35;

print "$0: translate center of selection to origin,\nand orient so that a vector of specified length in direction +x points away from the protein\nAlso rotate around x axis to maximise nr of atoms in cube with xyz axis sides\n";
if (! GetOptions('i=s' => \@in, 'o=s' => \@out, 'l=s' => \$ligandsel, 'd=f' => \$veclength)) {
  die("Usage: \n  -i [input files] coordinates.[gro|pdb|rst7] and (optional) topology.[parm7|top]\n  -o [output file] coordinates.pdb\n  -l VMD selection defining ligand / 'interesting' region\n  -d length of vector in angstrom to align along x axis for maximal protein distance (default 35)");
}

if (!(scalar(@in) && scalar(@out))) {
  die("Usage: \n  -i [input files] coordinates.[gro|pdb|rst7] and (optional) topology.[parm7|top]\n  -o [output file] coordinates.pdb\n  -l VMD selection defining ligand / 'interesting' region\n  -d length of vector in angstrom to align along x axis for maximal protein distance (default 35)");
}

my $topin;
my $coorin;
my $coorout;

foreach my $infile (@in) {
  if ($infile =~ /\.top$/) {
    die("only .parm7 accepted as topology file (VMD's fault)\n");
  }
  if ($infile =~ /\.parm7$/) {
    if (defined($topin)) {die("double definition of input topology\n"); }
    $topin = $infile;
  }
  if ($infile =~ /\.gro$/) {
    if (defined($topin)) {die("double defition of input coordinates\n");}
    $coorin = $infile;
  }
  elsif ($infile =~ /\.pdb$/) {
    if (defined($coorin)) {die("double defition of input coordinates\n");}
    $coorin = $infile;
  }
  elsif ($infile =~ /\.rst7$/) {
    if (defined($coorin)) {die("double defition of input coordinates\n");}
    $coorin = $infile;
  }
}
if (($coorin =~ /\.rst7$/) && (!defined($topin))) {
  die("rst7 input coordinates require a parm7 topology file\n");
}

foreach my $outfile (@out) {
  if ($outfile =~ /\.pdb$/) {
    if (defined($coorout)) {die("double defition of output coordinates\n");}
    $coorout = $outfile;
  }
}


if (!defined($coorout)) {die("no output file specified for coordinates\n");}
if (!defined($coorin)) {die("no input file specified for coordinates\n");}
if (!defined($ligandsel)) {die("no VMD selection defining ligand specified (-l)\n");}



my $vmdcmd = "";
if (defined($topin)) {
  $vmdcmd .= "mol new $topin\n";
  $vmdcmd .= "mol addfile $coorin\n";
}
else {
  $vmdcmd .= "mol new $coorin\n";
}
$vmdcmd .= "expr srand(1)\n";
$vmdcmd .= "set sel [atomselect top \"$ligandsel\"]\n";
$vmdcmd .= "set distance $veclength\n";
$vmdcmd .= "set bestvector {}\n";
$vmdcmd .= "set lowest 999999999\n";
#get geometric centre (wrongly labelled COM)\
$vmdcmd .= "set com [measure center \$sel]\n";
#translate COM to origin
$vmdcmd .= "[atomselect top \"all\"] move [transoffset [vecscale \$com -1]]\n";
# send random vectors out from origin
$vmdcmd .= "for {set i 0} {\$i < 1000} { incr i} {\n";
$vmdcmd .= "  if { ! [expr (\$i+1) % 100]} {puts [expr \$i+1]}\n";
$vmdcmd .= "  set randvec \"[expr rand() - 0.5] [expr rand() - 0.5] [expr rand() - 0.5]\"\n";
$vmdcmd .= "  set randvec [vecscale \$randvec [expr \$distance / [veclength \$randvec]]]\n";
$vmdcmd .= "  set sel [atomselect top \"(protein or nucleic) and (((sqr(x - [lrange \$randvec 0 0]  )+sqr(y-[lrange \$randvec 1 1])+sqr(z- [lrange \$randvec 2 2]))^0.5) < \$distance)\"]\n";
$vmdcmd .= "  set penalty [\$sel num]\n";
# penalty of natoms
$vmdcmd .= "  if {\$penalty < \$lowest} {\n";
$vmdcmd .= "    puts \"  score \$penalty\"\n";
$vmdcmd .= "    set lowest \$penalty\n";
$vmdcmd .= "    set bestvector \$randvec\n";
$vmdcmd .= "  }\n";
$vmdcmd .= "}\n";

# send random vectors out from best so far
$vmdcmd .= "for {set i 0} {\$i < 1000} { incr i} {\n";
$vmdcmd .= "  set randvec \"[expr [lrange \$bestvector 0 0] + rand() - 0.5] [expr [lrange \$bestvector 1 1] + rand() - 0.5] [expr [lrange \$bestvector 2 2] + rand() - 0.5]\"\n";
$vmdcmd .= "  set randvec [vecscale \$randvec [expr \$distance / [veclength \$randvec]]]\n";
$vmdcmd .= "  set sel [atomselect top \"(protein or nucleic) and (((sqr(x - [lrange \$randvec 0 0]  )+sqr(y-[lrange \$randvec 1 1])+sqr(z- [lrange \$randvec 2 2]))^0.5) < \$distance)\"]\n";
$vmdcmd .= "  set penalty [\$sel num]\n";
# penalty of natoms
$vmdcmd .= "  if {\$penalty < \$lowest} {\n";
$vmdcmd .= "    puts \"  score \$penalty\"\n";
$vmdcmd .= "    set lowest \$penalty\n";
$vmdcmd .= "    set bestvector \$randvec\n";
$vmdcmd .= "  }\n";
$vmdcmd .= "}\n";

# optimise (slow!)
$vmdcmd .= "set penalty 0\n";
$vmdcmd .= "foreach coord [\$sel get {x y z}] {\n";
$vmdcmd .= "  set penalty [expr \$penalty + 1/( ([lrange \$coord 0 0]-[lrange \$randvec 0 0])**2 + ([lrange \$coord 1 1]-[lrange \$randvec 1 1])**2 + ([lrange \$coord 2 2]-[lrange \$randvec 2 2])**2)]\n";
$vmdcmd .= "}\n";

$vmdcmd .= "for {set i 0} {\$i < 50} { incr i} {\n";
$vmdcmd .= "  set randvec \"[expr [lrange \$bestvector 0 0] + rand() - 0.5] [expr [lrange \$bestvector 1 1] + rand() - 0.5] [expr [lrange \$bestvector 2 2] + rand() - 0.5]\"\n";
$vmdcmd .= "  set randvec [vecscale \$randvec [expr \$distance / [veclength \$randvec]]]\n";
$vmdcmd .= "  set sel [atomselect top \"(protein or nucleic) and (((sqr(x - [lrange \$randvec 0 0]  )+sqr(y-[lrange \$randvec 1 1])+sqr(z- [lrange \$randvec 2 2]))^0.5) < \$distance)\"]\n";
$vmdcmd .= "  set penalty 0\n";
$vmdcmd .= "  foreach coord [\$sel get {x y z}] {\n";
# penalty of distance^-2
$vmdcmd .= "    set penalty [expr \$penalty + 1/( ([lrange \$coord 0 0]-[lrange \$randvec 0 0])**2 + ([lrange \$coord 1 1]-[lrange \$randvec 1 1])**2 + ([lrange \$coord 2 2]-[lrange \$randvec 2 2])**2)]\n";
$vmdcmd .= "  }\n";
$vmdcmd .= "  if {\$penalty < \$lowest} {\n";
$vmdcmd .= "    puts \"  iteration \$i score \$penalty\"\n";
$vmdcmd .= "    set lowest \$penalty\n";
$vmdcmd .= "    set bestvector \$randvec\n";
$vmdcmd .= "  }\n";
$vmdcmd .= "}\n";
$vmdcmd .= "puts \"Rotation matrix to put vector away from protein along x axis [transvecinv \$bestvector]\"\n";
$vmdcmd .= "[atomselect top \"all\"] move [transvecinv \$bestvector]\n";

my $halflength = $veclength/2;
$vmdcmd .= "set maxat [[atomselect top \"((protein or nucleic) and (x < $halflength && x > (-$halflength) && y < $halflength && y > (-$halflength) && z < $halflength && z > (-$halflength))) and not (same residue as ((x > $halflength || x < (-$halflength) || y > $halflength || y < (-$halflength) || z > $halflength || z < (-$halflength))))\"] num]\n";
$vmdcmd .= "set maxrot 0 \n";
$vmdcmd .= "for {set i 1} {\$i <= 360} { incr i} {\n";
$vmdcmd .= "  [atomselect top \"all\"] move [transaxis x 1]\n";
$vmdcmd .= "  set nat [[atomselect top \"((protein or nucleic) and (x < $halflength && x > (-$halflength) && y < $halflength && y > (-$halflength) && z < $halflength && z > (-$halflength))) and not (same residue as ((x > $halflength || x < (-$halflength) || y > $halflength || y < (-$halflength) || z > $halflength || z < (-$halflength))))\"] num]\n";
$vmdcmd .= "  if {\$nat > \$maxat} {\n";
$vmdcmd .= "    set maxat \$nat\n";
$vmdcmd .= "    set maxrot \$i\n";
$vmdcmd .= "  }\n";
$vmdcmd .= "}\n";
$vmdcmd .= "puts \"  best rotation \$maxrot\"\n";
$vmdcmd .= "[atomselect top \"all\"] move [transaxis x \$maxrot]\n";

# rotate in 90 degree increments so atom 1 is as far 'down' (in the y direction) as possible, for consistency between runs
$vmdcmd .= "set miny 99999\n";
$vmdcmd .= "set maxrot 0\n";
$vmdcmd .= "for {set i 0} {\$i < 4} {incr i} {\n";
$vmdcmd .= "  if {[[atomselect top \"serial 1\"] get y] < \$miny} {\n";
$vmdcmd .= "    set miny [[atomselect top \"serial 1\"] get y]\n";
$vmdcmd .= "    set maxrot \$i\n";
$vmdcmd .= "  }\n";
$vmdcmd .= "  [atomselect top \"all\"] move [transaxis x 90]\n\n";
$vmdcmd .= "}\n";
$vmdcmd .= "[atomselect top \"all\"] move [transaxis x [expr \$maxrot * 90]]\n\n";
$vmdcmd .= "\n";
$vmdcmd .= "\n";
$vmdcmd .= "\n";
$vmdcmd .= "\n";

$vmdcmd .= "[atomselect top \"all\"] writepdb $coorout\n";


#print $vmdcmd;
open(VMDIN,">vmdin");
print VMDIN $vmdcmd;
close(VMDIN);
`cat vmdin | ~/vmdtext `;


