#!/usr/bin/env perl
# this script reads in a gmx file, a pdb file and places the charges of
# each atom in the b factor column
# in the GMX the atom sstart at 0 while in the pdb they start form 1.
# Warning: this only works when there is one molecule in the tpr file
# otherwise, the numbering starts over from 1 and messes up the charges

$gmx = $ARGV[0];
$pdb = $ARGV[1];

$pdb =~ /(.+)\.pdb/;
$out = $1."_with_charges.pdb";

open(GMX,$gmx);
open(PDB,$pdb);
open(OUT,">$out");

# 1st read the charged from the gmx dump file
while(<GMX>){
  if(/            atom\[(.{6})\].{54} q=(.{12}),.+/){
    $id = $1;
    $charge = $2;
    $id =~ s/\s+//g;
    $charge =~ s/\s+//g;
    $charge += 0;
    $charges[$id] = $charge;
    print "$id $charge\n";
    #$total_charge += $charges;
  }
}
close(GMX);
# now assign the charges to the PDB
while(<PDB>){
  if(/^ATOM(.{7})(.{50}).*/){
    $id = $1;
    $else = $2;
    $id =~ s/\s+//g;
    print  OUT "ATOM";
    printf OUT ("%7i",$id);
    print  OUT $else;
    printf OUT ("%7.4f\n",$charges[$id-1]);
    #print "$charges[$id-1]\n";
    $total_charge += $charges[$id-1];
  }
  else{
    print OUT;
  }
}
close(PDB);
close(OUT);
print "$total_charge\n";









