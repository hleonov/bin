#!/usr/bin/perl
foreach $file (@ARGV) {
  $file =~ /(.+)\.pdb/; 
  $name = $1;
  $out = $name."_fixed.pdb";
  open(IN,"$file");
  if(-e "$out"){
    unlink("$out");
  }
  open(OUT,">>$out");
  $first = "true";
  while ($_ = <IN>) {
    if(/^ATOM/ || /^HETATM/) {
      /^(.{20}).{3}(.{47})(.{3})/;
      $beg = $1;
      $end = $2;
      $seg = $3;
      $seg =~ s/\s//g;
      $chain = &change;
      if($old_chain ne $chain && $first eq "false") { # its a new chain but not the first one
        print OUT "TER\n";
      }
      print OUT "$beg$chain$end\n";
      $first = "false";
      $old_chain = $chain;
    }
    else {
      print OUT;
    }
  }
  close(IN);
  close(OUT);
}
sub change {
  if ($seg == 1) {
    $new_chain = " A ";
  }
  elsif ($seg == 2) {
    $new_chain = " B ";
  }
  elsif ($seg == 3) {
    $new_chain = " C ";
  }
  elsif ($seg == 4) {
    $new_chain = " D ";
  }
  elsif ($seg == 5) {
    $new_chain = " E ";
  }
  elsif ($seg == 6) {
    $new_chain = " F ";
  }
  elsif ($seg == 7) {
    $new_chain = " G ";
  }
  elsif ($seg == 8) {
    $new_chain = " H ";
  }
  elsif ($seg == 9) {
    $new_chain = " I ";
  }
  elsif ($seg == 10) {
    $new_chain = " J ";
  }
  elsif ($seg == 11) {
    $new_chain = " K ";
  }
  elsif ($seg == 12) {
    $new_chain = " L ";
  }
  elsif ($seg == 13) {
    $new_chain = " M ";
  }
  elsif ($seg == 14) {
    $new_chain = " N ";
  }
  elsif ($seg == 15) {
    $new_chain = " O ";
  }
  elsif ($seg == 16) {
    $new_chain = " P ";
  }
  elsif ($seg == 17) {
    $new_chain = " Q ";
  }
  elsif ($seg == 18) {
    $new_chain = " R ";
  }
  elsif ($seg == 19) {
    $new_chain = " S ";
  }
  elsif ($seg == 20) {
    $new_chain = " T ";
  }
  return $new_chain;
}
