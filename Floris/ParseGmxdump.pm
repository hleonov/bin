package ParseGmxdump;
use strict;

use lib "/home/fbuelen/fec";
use ToputilsShared;

use Exporter;  
our @ISA = qw(Exporter);
our @EXPORT = qw(ParseGmxdump);

my $degtorad = 57.2957795;

sub ParseGmxdump {
  my $gmxdumpfilein = $_[0];
  open (GMXD, "$gmxdumpfilein");

  my $gmxdumpfile;

  my @settles;
  my @bonds;
  my @angles;
  my @pdihs;
  my @rbdihs;
  my @lj14s;
  my @posres;
  my @posrestypes;
  while (my $line = <GMXD>) {
    $gmxdumpfile .= $line;
    if ($line =~ / \(SETTLE\) /) {push @settles, $line;}
    elsif ($line =~ / \(LJ1[45][^\s\)]*\) /) {push @lj14s, $line;}
    elsif ($line =~ / \(BONDS\) /) {push @bonds, $line;}
    elsif ($line =~ / \(POSRES\) /) {push @posres, $line;}
    elsif ($line =~ /functype\[\s*\d+\]=POSRES/) {push @posrestypes, $line;}
    elsif ($line =~ / \(ANGLES\) /) {push @angles, $line;}
    elsif ($line =~ / \(PDIHS\) /) {push @pdihs, $line;}
    elsif ($line =~ / \(RBDIHS\) /) {push @rbdihs, $line;}
  }
  if (!defined($gmxdumpfile)) {
    die("gmxdumpfile not read\n");
  }


  
  my (%atoms, %bonds, %angles, %dihedrals, %bond_types, %angle_types, %dihedral_types, %nonbonded_par, %nonbonded_index, %ljsigeps);
  my %lj14_types;
  my (%posres_types, %posres);
  my %pairs;
  my $fudgeQQ; my $fudgeLJ;

  if ($gmxdumpfile =~ /moltype \(\d+\):/) {
    die("entry for 'moltype' found in gmxdump output; run gmxdump with option '-sys'\n");
  }
   
  if ($gmxdumpfile =~ /functype\[\s*\d+\s*\]=CONSTR/) {
    die("gmxdump file contains constraints - reprocess without constraints for correct bond and angle stuff\n");
  }

  $gmxdumpfile =~ /\n   atoms:(.*?)\n   \S/ms or die("No atoms section\n");
  my @atomsection = split("\n",$1);

  $gmxdumpfile =~ /\n   excls:(.*?)\n   \S/ms or die("No excls section\n");
  my @exclsection = split("}\n",$1);
  for (my $i = 0; $i < scalar(@exclsection); $i++) {
    $exclsection[$i] =~ s/\n//msg;
    $exclsection[$i] .= "}\n";
  }
  
  my $namecount = 0;
  my %resnames;
  my %uniquetypes;
  my %uniquetypeBs;

  foreach my $line (@atomsection) {
  
    if ($line =~ /         atom\[\s*(\d+)\]={type=\s*(\d+),\s*typeB=\s*(\d+),\s*ptype=\s*[^,]+,\s*m=\s*([^,]+),\s*q=\s*([^,]+),\s*mB=\s*([^,]+),\s*qB=\s*([^,]+),\s*res(nr|ind)=\s*(\d+),/ ) {
      my $atomnr = $1 + 1;
      my $type = $2;
      my $typeB = $3;
      my $mass = $4;
      my $charge = $5;
      my $massB = $6;
      my $chargeB = $7;
      my $resnr = $9;
      ${$atoms{$atomnr}}{atom_type_index} = $type+1;
      ${$atoms{$atomnr}}{mass} = $mass;
      ${$atoms{$atomnr}}{charge} = $charge;
      ${$atoms{$atomnr}}{resid} = $resnr+1;
      ${$atoms{$atomnr}}{atom_type_indexB} = $typeB+1;
      ${$atoms{$atomnr}}{massB} = $massB;
      ${$atoms{$atomnr}}{chargeB} = $chargeB;
    }
  
    if ($line =~ /         atom\[\s*(\d+)\]={name="([^"]+)"}/) {
      my $name = $2;
      my $nameB;
      my $nr = $1;
      if ($name =~ /(\S+)-x-(\S+)/) {
        $name = $1;
        $nameB = $2;
      }
      ${$atoms{$nr + 1}}{name} = $name;
      if (defined($nameB)) {
        ${$atoms{$nr + 1}}{nameB} = $nameB;
      }
      $namecount ++;
    }

    if ($line =~ /         type\[\s*(\d+)\]={name="([^"]+)",nameB="([^"]+)"}/) {
      my $type = $2;
      if (length($type) > 4) {
        $type =~ /^(.).*(...)$/;
        my $newtype = $1.$2;
        if (!defined($uniquetypes{$newtype})) {
          $uniquetypes{$newtype} = $type;
        }
        else {
          if (! ($uniquetypes{$newtype} eq $type)) {
            print "type $type $uniquetypes{$newtype}\n";
            die("clash with long atom types\n");
          }
        }
        #print "type $type newtype $newtype\n";
        $type = $newtype;
      }
  
      my $typeB = $3;
      if (length($typeB) > 4) {
        $typeB =~ /^(.).*(...)$/;
        my $newtypeB = $1.$2;
        if (!defined($uniquetypeBs{$newtypeB})) {
          $uniquetypeBs{$newtypeB} = $typeB;
        }
        else {
          if (! ($uniquetypeBs{$newtypeB} eq $typeB)) {
            print "typeB $typeB $uniquetypeBs{$newtypeB}\n";
            die("clash with long atom typeBs\n");
          }
        }
        #print "typeB $typeB newtypeB $newtypeB\n";
        $typeB = $newtypeB;
      }
  
      ${$atoms{$1 + 1}}{amber_atom_type} = $type;
      ${$atoms{$1 + 1}}{amber_atom_typeB} = $typeB;
    }
  

    if ($line =~ /         residue\[\s*(\d+)\]={name="([^"]+)"/) {
      $resnames{$1+1} = $2;
    }
  }

  print scalar(keys(%atoms)), " atoms\n";
  if ($namecount != scalar(keys(%atoms))) {die "More atoms than atom names\n";}
  
  foreach my $atom (keys(%atoms)) {
    ${$atoms{$atom}}{resname} = $resnames{${$atoms{$atom}}{resid}} or die ("No resname for atom $atom\n");
  }
  
  while ($gmxdumpfile =~ /      cgs\[\s*(\d+)\]={(\d+)\.\.(\d+)}/g) {
    for (my $i = $2 + 1; $i <= $3 + 1; $i++ ) {
      ${$atoms{$i}}{chargegroup} = $1 + 1;
    }
  }
  foreach my $atom (keys(%atoms)) {
    if (!defined(${$atoms{$atom}}{chargegroup})) {die("No charge group for atom $atom\n");}
    ${$atoms{$atom}}{number_excluded_atoms} = 0;
  }

  my %exclusions;
  foreach my $line (@exclsection) {
    $line =~ /      excls\[\s*(\d+)\]\[\s*\d+\.\.\d+\]={([^}]+)}/;
    $exclusions{$1 + 1} = $2;
  }
  foreach my $atom (keys(%exclusions)) {
    my $list = $exclusions{$atom};
    while ($list =~ /(\d+)/g) {
      ${${$atoms{$atom}}{excluded}}{$1 + 1} = 1;
      ${$atoms{$atom}}{number_excluded_atoms} ++;
    }
  }
  
  #foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
    #print "atom $atom name $atoms{$atom}{name} charge $atoms{$atom}{charge} mass $atoms{$atom}{mass} atom+type_index $atoms{$atom}{atom_type_index} amber_atom_tpe $atoms{$atom}{amber_atom_type} resid $atoms{$atom}{resid} resname $atoms{$atom}{resname} cg $atoms{$atom}{chargegroup}\n";
    #print "atom $atom exclusions ",join(" ",sort{$a <=> $b} keys(%{${$atoms{$atom}}{excluded}})),"\n";
  #}

  my $atomtypecount = 0;
  while ($gmxdumpfile =~ /      functype\[\s*(\d+)\]=LJ_SR,\s*c6=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*c12=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)/msg) {
    $atomtypecount ++;
  }
  if (sqrt($atomtypecount) != int(sqrt($atomtypecount))) {
    die("atom type count $atomtypecount is not a square\n");
  }
  $atomtypecount = sqrt($atomtypecount);

  print $atomtypecount, " LJ sigma / epsilon pairs\n";
  my $matchcount = 0;
  my $selfpaircount = 0;
  while ($gmxdumpfile =~ /      functype\[\s*(\d+)\]=LJ_SR,\s*c6=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*c12=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)/msg) {
    # if (1)
    if (($matchcount % $atomtypecount) == int($matchcount / $atomtypecount)) {
      my $atomtype = ($matchcount-int($matchcount / $atomtypecount)) / $atomtypecount + 1;
      my $A = $4;
      my $B = $2;
      #print "i:i at $1 A $A B $B\n";
      my $sigma; my $epsilon;
      if ($A == 0 || $B == 0) {
        $sigma = 0.1; $epsilon = 0;
      }
      else {
        $sigma = sqrt($A ** (1/3)) / sqrt($B ** (1/3));
        #$sigma = ($A/$B)**(1/6);
        $epsilon = $B*$B / (4*$A);
      }
      #printf( "atomtype %s sigma %11.5e epsilon %11.5e\n",$atomtype,$sigma,$epsilon);
      #if (($matchcount % $atomtypecount) == int($matchcount / $atomtypecount)) { print "*"; $selfpaircount ++;}

      ${$ljsigeps{$atomtype}}{sigma} = $sigma * 10;
      ${$ljsigeps{$atomtype}}{epsilon} = $epsilon / 4.184;
    }
    $matchcount ++;
  }
  #foreach my $atomtype (sort {$a <=> $b} keys(%ljsigeps)) {
  #  print "type $atomtype sig ${$ljsigeps{$atomtype}}{sigma} eps ${$ljsigeps{$atomtype}}{epsilon}\n";
  #}

  my $lastlj14_type = 0;
  my %lj14type_index;
  while ($gmxdumpfile =~ /      functype\[\s*(\d+)\]=LJ1[45][^,]*,\s*c6A=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*c12A=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*c6B=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*c12B=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)/msg) {
    $lj14type_index{$1} = $lastlj14_type;
    ${$lj14_types{$lastlj14_type}}{c6A} = $2;
    ${$lj14_types{$lastlj14_type}}{c12A} = $4;
    ${$lj14_types{$lastlj14_type}}{c6B} = $6;
    ${$lj14_types{$lastlj14_type}}{c12B} = $8;
    $lastlj14_type ++;
  }

  while ($gmxdumpfile =~ /\s+fudgeQQ\s+=\s*(\S*)/g) {
    if (defined($fudgeQQ)) {die("more than one fudgeQQ found\n");}
    $fudgeQQ = $1;
  }


  my $lastbondtype = 1;
  my $bondoutput = 0;
  my %bondindex;
  while ($gmxdumpfile =~ /      functype\[\s*(\d+)\]=BONDS,\s*b0A=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*cbA=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*b0B=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*cbB=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)/msg) {
    if (! $bondoutput) {print "Bond types .";} else {print".";} 
    $bondoutput ++; 
    #print "bondtype $1 $2 $4 $6 $8\n";
    $bondindex{$1} = $lastbondtype;
    ${$bond_types{$lastbondtype}}{equil_value} = $2 * 10;
    ${$bond_types{$lastbondtype}}{force_constant} = $4 / (4.184*200);
    ${$bond_types{$lastbondtype}}{equil_valueB} = $6 * 10;
    ${$bond_types{$lastbondtype}}{force_constantB} = $8 / (4.184*200);
    $lastbondtype ++;
  }
  if ($bondoutput) {print " ($bondoutput)\n";}

  my $lastsettletype = $lastbondtype;
  my $settleoutput = 0;
  while ($gmxdumpfile =~ /      functype\[\s*(\d+)\]=SETTLE,\s*doh=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*dhh=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)/msg) {
    if (! $settleoutput) {print "Settle types .";} else {print".";} 
    $settleoutput ++; 
    #print "settletype $1 $2 $4 $6 $8\n";
    $bondindex{$1} = $lastsettletype;
    ${$bond_types{$lastsettletype}}{equil_value} = $2 * 10;
    ${$bond_types{$lastsettletype}}{equil_valueB} = $2 * 10;
    ${$bond_types{$lastsettletype}}{force_constant} = 1;
    ${$bond_types{$lastsettletype}}{force_constantB} = 1;
    $lastsettletype ++;
  }
  if ($settleoutput) {print " ($settleoutput)\n";}

  my $lastangletype = 1;
  my $angleoutput = 0;
  my %angleindex;
  while ($gmxdumpfile =~ /      functype\[\s*(\d+)\]=ANGLES,\s*thA=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*ctA=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*thB=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*ctB=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)/msg) {
    if (! $angleoutput) {print "Angle types .";} else {print".";} 
    $angleoutput ++; 
    #print "angletype $1 $2 $4 $6 $8\n";
    $angleindex{$1} = $lastangletype;
    ${$angle_types{$lastangletype}}{equil_value} = $2 / $degtorad;
    ${$angle_types{$lastangletype}}{force_constant} = $4 / (4.184*2);
    ${$angle_types{$lastangletype}}{equil_valueB} = $6 / $degtorad;
    ${$angle_types{$lastangletype}}{force_constantB} = $8 / (4.184*2);
    $lastangletype ++;
  }
  if ($angleoutput) {print " ($angleoutput)\n";}
  
  my $lastdihedraltype = 1;
  my $pdihoutput = 0;
  my %dihedralindex;
  while ($gmxdumpfile =~ /      functype\[\s*(\d+)\]=PDIHS,\s*phiA=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*cpA=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*phiB=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*cpB=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*mult=\s*(\d+)/msg) {
    if (! $pdihoutput) {print "Periodic dihedral types .";} else {print".";} 
    $pdihoutput ++;
    #print "dihedraltype $1 $2 $4 $6 $8 $10\n";
    $dihedralindex{$1} = $lastdihedraltype;
    ${$dihedral_types{$lastdihedraltype}}{phase} = $2 / $degtorad;
    ${$dihedral_types{$lastdihedraltype}}{force_constant} = $4 / 4.184;
    ${$dihedral_types{$lastdihedraltype}}{periodicity} = $10;
    ${$dihedral_types{$lastdihedraltype}}{phaseB} = $6 / $degtorad;
    ${$dihedral_types{$lastdihedraltype}}{force_constantB} = $8 / 4.184;
    ${$dihedral_types{$lastdihedraltype}}{periodicityB} = $10;
    $lastdihedraltype ++;
  }
  if ($pdihoutput) {print " ($pdihoutput)\n";}
  
  my $rbdihoutput = 0;
  while ($gmxdumpfile =~ /      functype\[\s*(\d+)\]=RBDIHS,\s*rbcA\[0\]=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*rbcA\[1\]=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*rbcA\[2\]=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*rbcA\[3\]=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*rbcA\[4\]=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*rbcA\[5\]=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\n\s*rbcB\[0\]=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*rbcB\[1\]=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*rbcB\[2\]=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*rbcB\[3\]=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*rbcB\[4\]=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*rbcB\[5\]=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)/msg) {
    if (! $rbdihoutput) {print "Converting RB dihedrals to periodic .";} else {print".";}
    $rbdihoutput ++;
    #print "rbdihedraltype $1 $2 $4 $6 $8 $10 $12 $14 $16 $18 $20 $22 $24\n";
    my $rbnumber = $1;
    my @rbdihA = ($2,$4*(-1),$6,$8*(-1),$10,$12*(-1));
    my @pdihsequivA;
    pdihsfromrb(\@rbdihA, \@pdihsequivA);
    my @rbdihB = ($14,$16*(-1),$18,$20*(-1),$22,$24*(-1));
    my @pdihsequivB;
    pdihsfromrb(\@rbdihB, \@pdihsequivB);
    for (my $i = 1; $i <= 5; $i++) {
      if (($pdihsequivA[$i] != 0) || ($pdihsequivA[$i] != 0)) {
        if ($pdihsequivA[$i] <= 0) {
          ${$dihedral_types{$lastdihedraltype}}{phase} = 3.1415926535;
          ${$dihedral_types{$lastdihedraltype}}{force_constant} = sprintf("%10.4e",(0 - $pdihsequivA[$i]) / 4.184);
        }
        else {
          ${$dihedral_types{$lastdihedraltype}}{phase} = 0;          
          ${$dihedral_types{$lastdihedraltype}}{force_constant} = sprintf("%10.4e",$pdihsequivA[$i] / 4.184);
        }
        if ($pdihsequivB[$i] <= 0) {
          ${$dihedral_types{$lastdihedraltype}}{phaseB} = 3.1415926535;
          ${$dihedral_types{$lastdihedraltype}}{force_constantB} = sprintf("%10.4e",(0 - $pdihsequivB[$i]) / 4.184);
        }
        else {
          ${$dihedral_types{$lastdihedraltype}}{phaseB} = 0; 
          ${$dihedral_types{$lastdihedraltype}}{force_constantB} = sprintf("%10.4e",$pdihsequivB[$i] / 4.184);
        }
        ${$dihedral_types{$lastdihedraltype}}{periodicity} = $i;
        ${$dihedral_types{$lastdihedraltype}}{periodicityB} = $i;
        if (! defined(${$dihedral_types{$lastdihedraltype}}{phaseB})) {
          ${$dihedral_types{$lastdihedraltype}}{phaseB} = ${$dihedral_types{$lastdihedraltype}}{phase};
        }
        if (! defined(${$dihedral_types{$lastdihedraltype}}{force_constantB})) {
          ${$dihedral_types{$lastdihedraltype}}{force_constantB} = ${$dihedral_types{$lastdihedraltype}}{force_constant};
        }
        #print "  new dihdral type $lastdihedraltype fcA ${$dihedral_types{$lastdihedraltype}}{force_constant} fcB ${$dihedral_types{$lastdihedraltype}}{force_constantB} phaseA ${$dihedral_types{$lastdihedraltype}}{phase} phaseB ${$dihedral_types{$lastdihedraltype}}{phaseB}       periodicity ${$dihedral_types{$lastdihedraltype}}{periodicity}\n";
        push @{$dihedralindex{$rbnumber}}, $lastdihedraltype;
        $lastdihedraltype ++;
      }
    }
    #print "\n";
  }
  if ($rbdihoutput) {print " ($rbdihoutput)\n";}

  my $lastposrestype = 1;
  my $posresoutput = 0;
  my %posresindex;
  foreach my $line (@posrestypes) {
    if ($line =~ /      functype\[\s*(\d+)\]=POSRES, pos0A=\(\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\)\s*,\s*fcA=\(\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\)\s*,\s*pos0B=\(\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\)\s*,\s*fcB=\(\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?),\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\)/msg) {
      if (! $posresoutput) {print "Postion restraint types .";} else {print".";} 
      $posresoutput ++; 
      #print "posrestype $1 $2 $4 $6 $8 $10 $12 $14 $16 $18 $20 $22 $24\n";
      $posresindex{$1} = $lastposrestype;
      if (($8 != $10) || ($8 != $12) || ($20 != $22) || ($20 != $24)) {
      #if (($8 != $10) || ($8 != $12)) {
        if (! ($20 < -1000000)) {
          die("X, Y, Z in restraint force constant not equal - not supported here\n");
        }
      }
      ${$posres_types{$lastposrestype}}{force_constant} = $8;
      ${$posres_types{$lastposrestype}}{force_constantB} = $20;
      #print "lastposrestype $lastposrestype index $1\n";
      $lastposrestype ++;
    }
  }
  if ($posresoutput) {print " ($posresoutput)\n";}




  my $bondcount = 0;
  foreach my $line (@bonds) {
    $line =~ /\s*(\d+)\s+type\s*=\s*(\d+)\s+\(BONDS\)\s+(\d+)\s+(\d+)/;
    ${$bonds{$bondcount}}{a1} = $3 + 1;
    ${$bonds{$bondcount}}{a2} = $4 + 1;
    if (!defined($bondindex{$2})) {die ("No bondindex entry for bond type $2\n");}
    if ((${$bond_types{$bondindex{$2}}}{force_constant} != ${$bond_types{$bondindex{$2}}}{force_constantB}) || (${$bond_types{$bondindex{$2}}}{equil_value} != ${$bond_types{$bondindex{$2}}}{equil_valueB})) {
      ${$bonds{$bondcount}}{force_constantA} = ${$bond_types{$bondindex{$2}}}{force_constant};
      ${$bonds{$bondcount}}{equil_valueA} = ${$bond_types{$bondindex{$2}}}{equil_value};
      ${$bonds{$bondcount}}{force_constantB} = ${$bond_types{$bondindex{$2}}}{force_constantB};
      ${$bonds{$bondcount}}{equil_valueB} = ${$bond_types{$bondindex{$2}}}{equil_valueB};
    }
    else {
      ${$bonds{$bondcount}}{type} = $bondindex{$2}; 
    }
    $bondcount ++;
  }
  print "$bondcount non-water bonds\n";

  foreach my $line (@settles) {
    $line =~ /\s*(\d+)\s+type\s*=\s*(\d+)\s+\(SETTLE\)\s+(\d+)/;
    if (!defined($bondindex{$2})) {die ("No bondindex entry for settle type $2\n");}
    if ((${$bond_types{$bondindex{$2}}}{force_constant} != ${$bond_types{$bondindex{$2}}}{force_constantB}) || (${$bond_types{$bondindex{$2}}}{equil_value} != ${$bond_types{$bondindex{$2}}}{equil_valueB})) {
      die("settle entry differs between A and B states - not planned for\n");
    }
    ${$bonds{$bondcount}}{a1} = $3 + 1;
    ${$bonds{$bondcount}}{a2} = $3 + 2;
    ${$bonds{$bondcount}}{type} = $bondindex{$2}; 
    $bondcount ++;
    ${$bonds{$bondcount}}{a1} = $3 + 1;
    ${$bonds{$bondcount}}{a2} = $3 + 3;
    ${$bonds{$bondcount}}{type} = $bondindex{$2}; 
    $bondcount ++;
  }
  print "$bondcount bonds including settle (constraint bond force constant made up)\n";
  
  
  my %bondscountbyatom;
  foreach my $bond (keys(%bonds)) {
    $bondscountbyatom{${$bonds{$bond}}{a1}} ++;
    $bondscountbyatom{${$bonds{$bond}}{a2}} ++;
    #print "bond atoms ${$bonds{$bond}}{a1} ${$bonds{$bond}}{a2}\n";
  }
  foreach my $atom (keys(%atoms)) {
    ${$atoms{$atom}}{ishydrogen} = 0;
    if ((${$atoms{$atom}}{mass}) < 5 || (defined(${$atoms{$atom}}{massB}) && (${$atoms{$atom}}{massB} < 5) )) {
      if ($bondscountbyatom{$atom} == 1) {
        ${$atoms{$atom}}{ishydrogen} = 1;
      }
      else {
        print "suspicious light particle (atom $atom) not bound to exactly one other atom: mass ${$atoms{$atom}}{mass}, $bondscountbyatom{$atom} bonds\n";
        print "  ";
        foreach my $bond (keys(%bonds)) {
          if (${$bonds{$bond}}{a1} == $atom || ${$bonds{$bond}}{a2} == $atom) {
            print "${$bonds{$bond}}{a1}:${$bonds{$bond}}{a2}   ";
          }
        }
        print "\n";
      }
    }
  }
  
  my $anglecount = 0;
  foreach my $line (@angles) {
    $line =~ /\s*(\d+)\s+type\s*=\s*(\d+)\s+\(ANGLES\)\s+(\d+)\s+(\d+)\s+(\d+)/;
    ${$angles{$anglecount}}{a1} = $3 + 1;
    ${$angles{$anglecount}}{a2} = $4 + 1;
    ${$angles{$anglecount}}{a3} = $5 + 1;
    if (!defined($angleindex{$2})) {die ("No angleindex entry for angle type $2\n");}
    if ((${$angle_types{$angleindex{$2}}}{force_constant} != ${$angle_types{$angleindex{$2}}}{force_constantB}) || (${$angle_types{$angleindex{$2}}}{equil_value} != ${$angle_types{$angleindex{$2}}}{equil_valueB})) {
      ${$angles{$anglecount}}{force_constantA} = ${$angle_types{$angleindex{$2}}}{force_constant};
      ${$angles{$anglecount}}{equil_valueA} = ${$angle_types{$angleindex{$2}}}{equil_value};
      ${$angles{$anglecount}}{force_constantB} = ${$angle_types{$angleindex{$2}}}{force_constantB};
      ${$angles{$anglecount}}{equil_valueB} = ${$angle_types{$angleindex{$2}}}{equil_valueB};
    }
    else {
      ${$angles{$anglecount}}{type} = $angleindex{$2}; 
    }
    $anglecount ++;
  }
  print "$anglecount angles\n";
  
  my $dihedralcount = 0;
  foreach my $line (@pdihs) {
    $line =~ /\s*(\d+)\s+type\s*=\s*(\d+)\s+\(PDIHS\)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/;
    ${$dihedrals{$dihedralcount}}{a1} = $3 + 1;
    ${$dihedrals{$dihedralcount}}{a2} = $4 + 1;
    ${$dihedrals{$dihedralcount}}{a3} = $5 + 1;
    ${$dihedrals{$dihedralcount}}{a4} = $6 + 1;
    if (!defined($dihedralindex{$2})) {die ("No dihedralindex entry for dihedral type $2\n");}
    if (ref($dihedralindex{$2}))  {die "PDIH type is a reference (to array of multiple periodic dihedrals?); this isn't right\n";}
    if ((${$dihedral_types{$dihedralindex{$2}}}{force_constant} != ${$dihedral_types{$dihedralindex{$2}}}{force_constantB}) || (${$dihedral_types{$dihedralindex{$2}}}{periodicity} != ${$dihedral_types{$dihedralindex{$2}}}{periodicityB}) || (${$dihedral_types{$dihedralindex{$2}}}{phase} != ${$dihedral_types{$dihedralindex{$2}}}{phaseB})) {
      ${$dihedrals{$dihedralcount}}{force_constantA} = ${$dihedral_types{$dihedralindex{$2}}}{force_constant};
      ${$dihedrals{$dihedralcount}}{periodicityA} = ${$dihedral_types{$dihedralindex{$2}}}{periodicity};
      ${$dihedrals{$dihedralcount}}{phaseA} = ${$dihedral_types{$dihedralindex{$2}}}{phase};
      ${$dihedrals{$dihedralcount}}{force_constantB} = ${$dihedral_types{$dihedralindex{$2}}}{force_constantB};
      ${$dihedrals{$dihedralcount}}{periodicityB} = ${$dihedral_types{$dihedralindex{$2}}}{periodicityB};
      ${$dihedrals{$dihedralcount}}{phaseB} = ${$dihedral_types{$dihedralindex{$2}}}{phaseB};
    }
    else {
      ${$dihedrals{$dihedralcount}}{type} = $dihedralindex{$2}; 
    }
    $dihedralcount ++;
  }
  print "$dihedralcount dihedrals (PDIH type only)\n";

  foreach my $line (@rbdihs) {
    $line =~ /\s*(\d+)\s+type\s*=\s*(\d+)\s+\(RBDIHS\)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/;
    my $a1 = $3 + 1;
    my $a2 = $4 + 1;
    my $a3 = $5 + 1;
    my $a4 = $6 + 1;
    if (! (ref($dihedralindex{$2}) eq 'ARRAY')) {die "meltdown - RB dihedral type conversion expects an array of periodic dihedrals\n";}
    foreach my $pdih (@{$dihedralindex{$2}}) {
      ${$dihedrals{$dihedralcount}}{a1} = $a1;
      ${$dihedrals{$dihedralcount}}{a2} = $a2;
      ${$dihedrals{$dihedralcount}}{a3} = $a3;
      ${$dihedrals{$dihedralcount}}{a4} = $a4;
      if ((${$dihedral_types{$pdih}}{force_constant} != ${$dihedral_types{$pdih}}{force_constantB}) || (${$dihedral_types{$pdih}}{periodicity} != ${$dihedral_types{$pdih}}{periodicityB}) || (${$dihedral_types{$pdih}}{phase} != ${$dihedral_types{$pdih}}{phaseB})) {
        ${$dihedrals{$dihedralcount}}{force_constantA} = ${$dihedral_types{$pdih}}{force_constant};
        ${$dihedrals{$dihedralcount}}{periodicityA} = ${$dihedral_types{$pdih}}{periodicity};
        ${$dihedrals{$dihedralcount}}{phaseA} = ${$dihedral_types{$pdih}}{phase};
        ${$dihedrals{$dihedralcount}}{force_constantB} = ${$dihedral_types{$pdih}}{force_constantB};
        ${$dihedrals{$dihedralcount}}{periodicityB} = ${$dihedral_types{$pdih}}{periodicityB};
        ${$dihedrals{$dihedralcount}}{phaseB} = ${$dihedral_types{$pdih}}{phaseB};
      }
      else {
        ${$dihedrals{$dihedralcount}}{type} = $pdih;
      }
      $dihedralcount ++;
    }
  }
  print "$dihedralcount dihedrals (PDIH+RB types)\n";


  my $posrescount = 0;
  print "scaalar aray posres ",scalar(@posres),"\n";
  foreach my $line (@posres) {
    $line =~ /\s*(\d+)\s+type\s*=\s*(\d+)\s+\(POSRES\)\s+(\d+)/;
    ${$posres{$posrescount}}{a1} = $3 + 1;
    if (!defined($posresindex{$2})) {die ("No posresindex entry for posres type $2\n");}
    #if ((${$posres_types{$posresindex{$2}}}{force_constant} != ${$posres_types{$posresindex{$2}}}{force_constantB})) {
      ${$posres{$posrescount}}{force_constantA} = ${$posres_types{$posresindex{$2}}}{force_constant} / 836.8;
      if (${$posres_types{$posresindex{$2}}}{force_constantB} < -1000000) {
        ${$posres{$posrescount}}{fbhw} = ${$posres_types{$posresindex{$2}}}{force_constantB};
        ${$posres{$posrescount}}{force_constantB} = ${$posres_types{$posresindex{$2}}}{force_constant} / 836.8;
      }
      else {
        ${$posres{$posrescount}}{force_constantB} = ${$posres_types{$posresindex{$2}}}{force_constantB} / 836.8;
      }
    #}
    #else {
    #  ${$posres{$posrescount}}{type} = $posresindex{$2}; 
    #}
    $posrescount ++;
  }
  print "$posrescount position restraints\n";

  foreach my $r (keys(%posres)) {
    my $a = ${$posres{$r}}{a1};
    my $fcA = ${$posres{$r}}{force_constantA};
    my $fcB = ${$posres{$r}}{force_constantB};
    if (!defined($atoms{$a})) {
      die("tried to add a restraint to non-defined atom $a\n");
    }
    ${$atoms{$a}}{restraint_k} = $fcA;
    if (defined(${$posres{$r}}{fbhw})) {
      ${$atoms{$a}}{restraint_fbhw} = ${$posres{$r}}{fbhw};
    }
    else {
      ${$atoms{$a}}{restraint_kB} = $fcB;
    }
  }
  

  
  my $paircount = 0;
  my $specialtypes = 0;
  foreach my $line (@lj14s) {
    $line =~ /\s*(\d+)\s+type\s*=\s*(\d+)\s+\((LJ1[45][^\)]*)\)\s+(\d+)\s+(\d+)/;
    # enter explicit pairs, and also just check consistency between listed LJ14 ftypes and atom sigma/epsilon
    # I assume there's a single fudgeLJ parameter governing everything; if not, throw toys out of the pram
    my $f = $2;
    my $type = $3;
    my $a1 = $4 + 1;
    my $a2 = $5 + 1;
    my $pairkey = join(":",(sort {$a <=> $b} ($a1, $a2)));
    if ($type eq 'LJ14') {$type = 1;}
    elsif ($type eq 'LJ14_ASTATEONLY') {$type = 3;}
    elsif ($type eq 'LJ14_BSTATEONLY') {$type = 4;}
    elsif ($type eq 'LJ15_ASTATEONLY') {$type = 5;}
    elsif ($type eq 'LJ15_BSTATEONLY') {$type = 6;}
    else {die("unrecognised pair type $type\n");}
    if ($type > 3) {$specialtypes ++;}
    my $type1A = ${$atoms{$a1}}{atom_type_index};
    my $type2A = ${$atoms{$a2}}{atom_type_index};
    my $type1B = ${$atoms{$a1}}{atom_type_indexB};
    my $type2B = ${$atoms{$a2}}{atom_type_indexB};
    my $sigma1A = ${$ljsigeps{$type1A}}{sigma};
    my $epsilon1A = ${$ljsigeps{$type1A}}{epsilon};
    my $sigma2A = ${$ljsigeps{$type2A}}{sigma};
    my $epsilon2A = ${$ljsigeps{$type2A}}{epsilon};
    my $sigma1B = ${$ljsigeps{$type1B}}{sigma};
    my $epsilon1B = ${$ljsigeps{$type1B}}{epsilon};
    my $sigma2B = ${$ljsigeps{$type2B}}{sigma};
    my $epsilon2B = ${$ljsigeps{$type2B}}{epsilon};
    my $sijA = 0.5*($sigma1A+$sigma2A) / 10;
    my $eijA = sqrt($epsilon1A*$epsilon2A) * 4.184;
    my $c12A = 4*$eijA*($sijA**12);
    my $c6A = 4*$eijA*($sijA**6);
    my $sijB = 0.5*($sigma1B+$sigma2B) / 10;
    my $eijB = sqrt($epsilon1B*$epsilon2B) * 4.184;
    my $c12B = 4*$eijB*($sijB**12);
    my $c6B = 4*$eijB*($sijB**6);

    ${$pairs{$pairkey}}{$type} = 1;

    if ($c12A && $c6A) {
      my $div12A = sprintf("%5f",${$lj14_types{$lj14type_index{$2}}}{c12A} / $c12A);
      my $div6A = sprintf("%5f",${$lj14_types{$lj14type_index{$2}}}{c6A} / $c6A);
      #print "$div12A $div6A\n";
      if (!defined($fudgeLJ)) {$fudgeLJ = $div12A;}
      if (($type == 1) || ($type == 3)) {
        if ($div12A != $fudgeLJ) {die("Inconsistency in fudgeLJ A ($div12A $fudgeLJ $type) - might be ok (explicitly specified pairs or combination rule other than 2) but not coded for\n");}
        if ($div6A != $fudgeLJ) {die("Inconsistency in fudgeLJ A ($div6A $fudgeLJ $type) - might be ok (explicitly specified pairs or combination rule other than 2) but not coded for\n");}
      }
      if ($type == 5) { # should be full-strength A state
        if ((abs($div12A-1) > 0.00001) || (abs($div6A-1) > 0.00001)) {
          die("Inconsistency in fudgeLJ A ($div12A $fudgeLJ $type) for LJ15 type\n");
        }
      }
    }
    
    if ($c12B && $c6B) {
      my $div12B = sprintf("%5f",${$lj14_types{$lj14type_index{$2}}}{c12B} / $c12B);
      my $div6B = sprintf("%5f",${$lj14_types{$lj14type_index{$2}}}{c6B} / $c6B);
      #print "$div12B $div6B\n";
      if (!defined($fudgeLJ)) {$fudgeLJ = $div12B;}
      if (($type == 1) || ($type == 4)) {
        if ($div12B != $fudgeLJ) {die("Inconsistency in fudgeLJ B ($div12B $fudgeLJ $type) - might be ok (explicitly specified pairs or combination rule other than 2) but not coded for\n");}
        if ($div6B != $fudgeLJ) {die("Inconsistency in fudgeLJ B ($div6B $fudgeLJ $type) - might be ok (explicitly specified pairs or combination rule other than 2) but not coded for\n");}
      }
      if ($type == 6) { # should be full-strength B state
        if ((abs($div12B-1) > 0.00001) || (abs($div6B-1) > 0.00001)) {
          die("Inconsistency in fudgeLJ A ($div12B $fudgeLJ $type) for LJ15 type\n");
        }
      }
    }

    

    #print "from sigeps $c12A $c6A $c12B $c6B\n";
    #print "from ftype ${$lj14_types{$lj14type_index{$2}}}{c12A} ${$lj14_types{$lj14type_index{$2}}}{c6A} ${$lj14_types{$lj14type_index{$2}}}{c12B} ${$lj14_types{$lj14type_index{$2}}}{c6B}\n";
    #print "  div ",$c12A / ${$lj14_types{$lj14type_index{$2}}}{c12A};
    #print "  ", $c6A / ${$lj14_types{$lj14type_index{$2}}}{c6A} ;
    #print "  ", $c12B / ${$lj14_types{$lj14type_index{$2}}}{c12B};
    #print "  ", $c6B / ${$lj14_types{$lj14type_index{$2}}}{c6B}, "\n";
    #my $sigmafromftype = sqrt(${$lj14_types{$lj14type_index{$2}}}{c12A} ** (1/3)) / sqrt(${$lj14_types{$lj14type_index{$2}}}{c6A} ** (1/3));
    #my $epsilonfromftype = ${$lj14_types{$lj14type_index{$2}}}{c6A}*${$lj14_types{$lj14type_index{$2}}}{c6A} / (4*${$lj14_types{$lj14type_index{$2}}}{c12A});
    #print "from ftype sigma $sigmafromftype epsilon $epsilonfromftype\n";
    
  }
  print "fudgeLJ $fudgeLJ\n";
  #print "still need to add the pairs\n";


  my @return;
  push @return, \%atoms, \%bonds, \%angles, \%dihedrals, \%bond_types, \%angle_types, \%dihedral_types, \%nonbonded_par, \%nonbonded_index, \%ljsigeps, "gmxdump.mdp";
  #if ($specialtypes) {
  #  print "pair types LJ14_ASTATEONLY, LJ14_BSTATEONLY, LJ15_ASTATEONLY or LJ15_BSTATEONLY found; these will take precedence over the automatically generated ones in WriteGTop.pm, take care...\n";
    push @return, \%pairs;
  #}
  return @return;

#  exit(0);
}

sub pdihsfromrb {
  # i'm too stupid to do this analytically so here's a messy brute force iterative scheme
  # to figure out the periodic dihedral equivalents for a Ryckaert-Bellemans dihedral
  # Working backwards from an RB to pdihs routine I found from a charmm to gtop conversion script,
  # guess and refine periodic dihedrals until we match the input RB
  my $pi = 3.1415926535;
  my @rbin = @{$_[0]};
  my @pdout = @{$_[1]};
  for (my $i = 1; $i <= 5; $i++) {
    $pdout[$i] = 0;
  }
  my $ndihs;
  my $nlevels = 1;
  
  my $diff = 99999;
  my $stepsize = 0.1;
  my $failcount = 0;
  my $failcount2 = 0;
  while (abs($diff) > 1e-9) {
    #print "diff $diff\n";
    my @new;
    for (my $i = 1; $i <= 5; $i++) {
      $new[$i] = $pdout[$i];
    }
    my $changei = int(rand()*($nlevels))+1;
    my $sign = 1; if (rand()>0.5) {$sign = -1;}
    $new[$changei] += $sign*$stepsize;
    my $newdiff = rbfrompdihs(\@new, \@rbin);
    #print "diff $diff newdiff $newdiff\n";
    if ($newdiff < $diff) {
      $diff = $newdiff;
      for (my $i = 1; $i <= 5; $i++) {
        $pdout[$i] = $new[$i];
      }
      $failcount = 0;
    }
    else {$failcount ++;}
    if ($failcount > 20) {
      $stepsize /= 3;
      $failcount = 0;
    }
    if ($stepsize < 1e-10) {
      $diff = 99999;
      $stepsize = 0.1;
      if ($nlevels < 5) {$nlevels ++;}
      for (my $i = 1; $i <= 5; $i++) {
        $new[$i] = 0;
      }
      for (my $i = 1; $i <= $nlevels; $i++) {
        $new[$i] = rand()*10;
      }
    }
  }
  for (my $i = 1; $i <= 5; $i++) {
    if (abs($pdout[$i]) < 0.0001) {$pdout[$i] = 0;}
    ${$_[1]}[$i] = $pdout[$i];
  }
}

sub rbfrompdihs {
  my @newp = @{$_[0]};
  my @newrb = (0,0,0,0,0,0);
  my @in = @{$_[1]};
  for (my $i = 1; $i <= 5; $i++) {
    if ( $i == 1 ) {
      $newrb[1] += $newp[$i];
    } elsif ( $i == 2 ) {
      $newrb[2] += 2*$newp[$i];
    } elsif ( $i == 3 ) {
      $newrb[1] += -3*$newp[$i];
      $newrb[3] += 4*$newp[$i];
    } elsif ( $i == 4 ) {
      $newrb[2] += -8*$newp[$i];
      $newrb[4] += 8*$newp[$i];
    } elsif ( $i == 5 ) {
      $newrb[1] += 5*$newp[$i];
      $newrb[3] += -20*$newp[$i];
      $newrb[5] += 16*$newp[$i];
    }
  }

  my $diff = 0;
  for (my $i = 1; $i <= 5; $i++) {
    $diff += ($newrb[$i] - $in[$i]) ** 2;
  }
  return abs($diff);  
}
1;
