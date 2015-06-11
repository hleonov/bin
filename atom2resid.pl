#!/usr/bin/perl
undef(%resid_hash);
undef($chain);
$message = "Usage: perl $0 <contacts.dat> <pdb> <out> <low freq. cutoff> <high cutoff> [noChain]";

$in = shift || die $message;
$pdb = shift || die $message;
$out = shift || die $message;
$low = shift || die $message;
$high = shift || die $message;
$no_chain = shift;

$col = ($no_chain ? 5 : 6);
open(IN,$in) || die "cannot open input file $in\n";
open(OUT, ">$out") || die "cannot open output file $out\n";
while (<IN>) {
   @parts = split();
 #  print $parts[0]," ", $parts[1],"\n"; 
   $grep_res1 = `grep \"ATOM\\s\\s*$parts[0]\\s\" $pdb | awk '{print \$4, \$${col}}'`;
 #  print `grep \"ATOM\\s\\s*$parts[0]\\s\" $pdb`;
   $grep_res2 = `grep \"ATOM\\s\\s*$parts[1]\\s\" $pdb | awk '{print \$4, \$${col}}'`;
 #  print `grep \"ATOM\\s\\s*$parts[1]\\s\" $pdb`;
   $grep_res1 =~ s/\n//g;
   $grep_res2 =~ s/\n//g;
   
   $resid_hash{$grap_res1}{$grep_res2} += $parts[2];   
   if ($parts[2] <= $high && $parts[2] >= $low) {
      print OUT "$grep_res1 \t $grep_res2 \t $parts[2]\n";
   }
}
close(IN);
close(OUT);
