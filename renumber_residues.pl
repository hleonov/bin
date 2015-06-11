#!/usr/bin/env perl
require "/home/hleonov/bin/shy_perl_stuff.pl";
use Getopt::Long qw (:config default);
$mesaage = "Usage:\nrenumber_residues.pl -f infile [-o outfile] [-n number_to_add] [-h]
For Example:
renumber_residues.pl -f system.pdb -o system_fixed.pdb -n 10\n";
# defaults
$n = 0;
# parse options
GetOptions ("f=s" => \$infile,  # string 	 
            "o=s" => \$outfile,	# string 	 
            "n=i" => \$n,	# integer 	 
            "h"   => \$help)	# flag	 	 
	    || die $mesaage;	    	
if ($help) {
  print $mesaage;
  exit;
}
unless ($outfile){
  $infile =~ /(.+)\.pdb/;
  $outfile = $1."_fixed.pdb";
}


open (IN,  &in_file($infile));
if($outfile){
  open (OUT, ">".&out_file($outfile));
}
while(<IN>){
  $global_counter++;
  &line_message($global_counter,1000);
  if(/^(ATOM  |HETATM)(.{15})(.)(....)(.+)/){
    $a     = $1;
    $b     = $2;
    $chain = $3;
    $resid = $4;
    $d     = $5;
    if($old_chain ne $chain){
      $resid_modifier = $resid;
    }
    $resid = $resid + 1 - $resid_modifier;
    $resid += $n;
    print OUT $a;
    print OUT $b;
    print OUT $chain;
    printf OUT ("%4.0f",$resid);
    print OUT "$d\n";
    $old_chain = $chain;
  }
  else{
    print OUT;
  }
}
close(IN);
close(OUT);











