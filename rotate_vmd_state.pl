#!/usr/local/perl/bin/perl

# This script will rotate a pdb according to a vmd state file.
# so rotate the pdb in vmd, save the state and then apply this script.
# If you want you can also supply the rotation matrix yourslef.



%argh = @ARGV ;
if(defined($argh{"-h"})){
  print "Usage -f my_protien.pdb -v myprotein.vmd [-o outfile.pdb] [-h] [-matrix a11_a12_a23_a21_a22_a23_a31_a32_a33]\n";
  exit;
}
if (defined($argh{"-f"})){
  $pdbfile = $argh{"-f"} ;
  open(INPDB,$pdbfile) || die "Input pdb file \"$infile\" does not exist\n";
  unless($pdbfile =~ /.+\.pdb/){
    print "Input pdb file $pdbfile doesn't contain an .pdb ending\n";
    print "Usage -f my_protien.pdb -v myprotein.vmd [-h] [-matrix a11 a12 a23 a21 a22 a23 a31 a32 a33]\n";
    exit;
  }
}
else {
  print "Input pdb file not set (-f)\n" ;
  print "Usage -f my_protien.pdb -v myprotein.vmd [-o outfile.pdb] [-h] [-matrix a11_a12_a23_a21_a22_a23_a31_a32_a33]\n";
  exit;
}
if (defined($argh{"-o"})){
  $outfile = $argh{"-o"} ;
}
else{
  # take care of opening the outputfile
  $pdbfile =~ /(.+)\.pdb/;
  $name = $1;
  $outfile = $name . "_rotated.pdb";
}
if (defined($argh{"-v"})){
  $vmdfile = $argh{"-v"} ;
  open(INVMD,$vmdfile) || die "Input vmd file \"$vmdfile\" does not exist\n";
}
else {
  print "Input vmd file not set (-v)\n" ;
  print "Usage -f my_protien.pdb -v myprotein.vmd [-o outfile.pdb] [-h] [-matrix a11_a12_a23_a21_a22_a23_a31_a32_a33]\n";
  exit;
}
if (defined($argh{"-matrix"})){
  $matrix = $argh{"-matrix"};
  ($a11, $a12, $a13, $a21, $a22, $a23, $a31, $a32, $a33) = split(/_/,$matrix);
}
else{
  # extarct matrix
  while(<INVMD>){
    if(/^set viewpoints\(\[molinfo top\]\) (.+)/){ # it is the matric info line
      $matrices = $1;
      $matrices =~ s/\{//g;      
      $matrices =~ s/\}//g;
      # now we are just left with numbers
      @n = split(/\s+/,$matrices);
      # the rotation matrix is
      ($a11, $a12, $a13, $a21, $a22, $a23, $a31, $a32, $a33) = ($n[16], $n[17], $n[18], $n[20], $n[21], $n[22], $n[24], $n[25], $n[26]); 
    }
  }
  close(INVMD);
}
# take care of opening the outputfile
$outfile = &overwrite($outfile);
open(OUT,">$outfile");


while(<INPDB>){
  if(/^ATOM/ || /^HETATM/){
    /^(.{30})(.{8})(.{8})(.{8})(.*)/;
    $before = $1;
    $x = $2;
    $y = $3;
    $z = $4;
    $after = $5;
    $x =~ s/\s+//g;    
    $y =~ s/\s+//g;    
    $z =~ s/\s+//g;
    $new_x = ($a11 * $x) + ($a12 * $y) + ($a13 * $z);
    $new_y = ($a21 * $x) + ($a22 * $y) + ($a23 * $z);
    $new_z = ($a31 * $x) + ($a32 * $y) + ($a33 * $z);
    print  OUT $before;
    printf OUT  (("%8.3f",$new_x));
    printf OUT  (("%8.3f",$new_y));
    printf OUT  (("%8.3f",$new_z));
    print  OUT  "$after\n";
  }
  else{
    print OUT;
  }
}

sub overwrite {
  my $file = $_[0];
  while(-e $file){
    print "File \"$file\" exists, overwrite? [y]\n";
    chop($response = <STDIN>);
    if($response eq "y"){
      last;
    }
    else{
      print "Please suplly output file name.\n";
      chop($file = <STDIN>);
    }
  }
  return $file;
}






















