#!/usr/bin/perl
#Usage: perl $0 system_from_MD.trr  ../system.ndx system_for_MD.tpr

$trj   	= shift || die "Usage: perl $0 <trr/xtc> <index> <tpr>\n";
#$init  	= shift || die "Usage: perl $0 <trr/xtc> <initial pdb/gro> <index> <tpr>\n";
$index 	= shift || die "Usage: perl $0 <trr/xtc> <index> <tpr>\n";
$tpr	= shift || die "Usage: perl $0 <trr/xtc> <index> <tpr>\n";	

$nojump	= "nojump.xtc";
$whole	= "whole.xtc";
$final	= "final_fit.xtc";

$com1 = "printf \"1\\n0\\n\" | trjconv -f $trj -o $nojump -pbc nojump -center rect -s $tpr -n $index";
$com2 = "printf \"0\\n\" | trjconv -ur rect -pbc whole -f $nojump -o $whole -s $tpr";
$com3 = "printf \"1\\n0\\n\" | trjconv -fit rot+trans -f $whole -o $final -s $tpr";

system($com1) ==0 || die "Cannot exec nojump\n $!\n";
system($com2) ==0 || die "Cannot exec whole\n $!\n";
system($com3) ==0 || die "Cannot exec fit\n $!\n";

system("rm $nojump $whole") == 0 || die "Cannot remove temporary files\n";
