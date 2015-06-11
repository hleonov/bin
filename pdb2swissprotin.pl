#!/usr/bin/env perl

# this reads in the cath list and add another entry for the swiss prot ID
# the data are from:
# http://release.cathdb.info/v3.2.0/CathDomainPdb.S35.v3.2.0.tgz
# http://release.cathdb.info/v3.2.0/CathDomainPdb.S60.v3.2.0.tgz
# http://www.bioinf.org.uk/pdbsws/pdbsws_chain.txt

# first read in the mapping
open (IN,"data/pdb2swiss.txt");
while (<IN>) {
  ($pdb, $chain, $swissprot) = split;
  $map{$pdb} = $swissprot;
}
close(IN);
# now map them
open (IN,"data/cath_35_list.txt");
open (OUT,">data/cath_35_list_swiss_prot.txt");
while (<IN>) {
  chop;
  /(^....).+/;
  $pdb = $1;
  $entry = $_;
  if ($map{$pdb} && $map{$pdb} ne "?") {
    print OUT "$entry $map{$pdb}\n";
  } else {
    print OUT "$entry 000000\n";
    print STDERR "No match for $pdb\n";
  }
}
close(IN);
close(OUT);

open (IN,"data/cath_60_list.txt");
open (OUT,">data/cath_60_list_swiss_prot.txt");
while (<IN>) {
  chop;
  /(^....).+/;
  $pdb = $1;
  $entry = $_;
  if ($map{$pdb} && $map{$pdb} ne "?") {
    print OUT "$entry $map{$pdb}\n";
  } else {
    print OUT "$entry 000000\n";
    print STDERR "No match for $pdb\n";
  }
}
close(IN);
close(OUT);

