package ParseParm7;
use strict;

use lib "/home/fbuelen/fec";
use ToputilsShared;

use Exporter;  
our @ISA = qw(Exporter);
our @EXPORT = qw(ParseParm7);

sub ParseParm7 {
  my $parmfile = $_[0];
  #print $parmfile;
  $parmfile =~ /FLAG ATOM_NAME\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG ATOM_NAME not found\n";
  my $ATOM_NAME = $1;
  $parmfile =~ /FLAG CHARGE\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG CHARGE not found\n";
  my $CHARGE = $1;
  $parmfile =~ /FLAG MASS\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG MASS not found\n";
  my $MASS = $1;
  $parmfile =~ /FLAG ATOM_TYPE_INDEX\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG ATOM_TYPE_INDEX not found\n";
  my $ATOM_TYPE_INDEX = $1;
  $parmfile =~ /FLAG NUMBER_EXCLUDED_ATOMS\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG NUMBER_EXCLUDED_ATOMS not found\n";
  my $NUMBER_EXCLUDED_ATOMS = $1;
  $parmfile =~ /FLAG NONBONDED_PARM_INDEX\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG NONBONDED_PARM_INDEX not found\n";
  my $NONBONDED_PARM_INDEX = $1;
  $parmfile =~ /FLAG LENNARD_JONES_ACOEF\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG LENNARD_JONES_ACOEF not found\n";
  my $LENNARD_JONES_ACOEF = $1;
  $parmfile =~ /FLAG LENNARD_JONES_BCOEF\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG LENNARD_JONES_BCOEF not found\n";
  my $LENNARD_JONES_BCOEF = $1;
  $parmfile =~ /FLAG AMBER_ATOM_TYPE\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG AMBER_ATOM_TYPE not found\n";
  my $AMBER_ATOM_TYPE = $1;
  $parmfile =~ /FLAG BONDS_INC_HYDROGEN\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG BONDS_INC_HYDROGEN not found\n";
  my $BONDS_INC_HYDROGEN = $1;
  $parmfile =~ /FLAG BONDS_WITHOUT_HYDROGEN\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG BONDS_WITHOUT_HYDROGEN not found\n";
  my $BONDS_WITHOUT_HYDROGEN = $1;
  $parmfile =~ /FLAG BOND_FORCE_CONSTANT\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG BOND_FORCE_CONSTANT not found\n";
  my $BOND_FORCE_CONSTANT = $1;
  $parmfile =~ /FLAG BOND_EQUIL_VALUE\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG BOND_EQUIL_VALUE not found\n";
  my $BOND_EQUIL_VALUE = $1;
  $parmfile =~ /FLAG ANGLES_INC_HYDROGEN\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG ANGLES_INC_HYDROGEN not found\n";
  my $ANGLES_INC_HYDROGEN = $1;
  $parmfile =~ /FLAG ANGLES_WITHOUT_HYDROGEN\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG ANGLES_WITHOUT_HYDROGEN not found\n";
  my $ANGLES_WITHOUT_HYDROGEN = $1;
  $parmfile =~ /FLAG ANGLE_FORCE_CONSTANT\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG ANGLE_FORCE_CONSTANT not found\n";
  my $ANGLE_FORCE_CONSTANT = $1;
  $parmfile =~ /FLAG ANGLE_EQUIL_VALUE\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG ANGLE_EQUIL_VALUE not found\n";
  my $ANGLE_EQUIL_VALUE = $1;
  $parmfile =~ /FLAG DIHEDRALS_INC_HYDROGEN\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG DIHEDRALS_INC_HYDROGEN not found\n";
  my $DIHEDRALS_INC_HYDROGEN = $1;
  $parmfile =~ /FLAG DIHEDRALS_WITHOUT_HYDROGEN\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG DIHEDRALS_WITHOUT_HYDROGEN not found\n";
  my $DIHEDRALS_WITHOUT_HYDROGEN = $1;
  $parmfile =~ /FLAG DIHEDRAL_FORCE_CONSTANT\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG DIHEDRAL_FORCE_CONSTANT not found\n";
  my $DIHEDRAL_FORCE_CONSTANT = $1;
  $parmfile =~ /FLAG DIHEDRAL_PERIODICITY\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG DIHEDRAL_PERIODICITY not found\n";
  my $DIHEDRAL_PERIODICITY = $1;
  $parmfile =~ /FLAG DIHEDRAL_PHASE\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG DIHEDRAL_PHASE not found\n";
  my $DIHEDRAL_PHASE = $1;
  $parmfile =~ /FLAG RESIDUE_LABEL\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG RESIDUE_LABEL not found\n";
  my $RESIDUE_LABEL = $1;
  $parmfile =~ /FLAG RESIDUE_POINTER\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG RESIDUE_POINTER not found\n";
  my $RESIDUE_POINTER = $1;
  $parmfile =~ /FLAG EXCLUDED_ATOMS_LIST\s*[\n]*\n.*\n([^%]+)\%/m or die "FLAG EXCLUDED_ATOMS_LIST not found\n";
  my $EXCLUDED_ATOMS_LIST = $1;

  my (%atoms, %bonds, %angles, %dihedrals, %bond_types, %angle_types, %dihedral_types, %nonbonded_par, %nonbonded_index, %ljsigeps);
  my $count = 1;
  while ($ATOM_NAME =~ /(....)\n?/mg) {
    my $var = trim($1);
    if ($var =~ /^\S+$/) {
      ${$atoms{$count}}{name} = $var;
      $count ++;
    }
  }
  $count = 1;
  while ($CHARGE =~ /(................)\n?/mg) {
    my $var = trim($1);
    if ($var =~ /^\S+$/) {
      ${$atoms{$count}}{charge} = $var / 18.2223;
      $count ++;
    }
  }
  $count = 1;
  while ($MASS =~ /(................)\n?/mg) {
    my $var = trim($1);
    if ($var =~ /^\S+$/) {
      ${$atoms{$count}}{mass} = $var;
      $count ++;
    }
  }
  $count = 1;
  while ($ATOM_TYPE_INDEX =~ /(........)\n?/mg) {
    my $var = trim($1);
    if ($var =~ /^\S+$/) {
      ${$atoms{$count}}{atom_type_index} = $var;
      #if (!defined( $ljsigeps{$var})) {
      #  my %blankhash = ();
      #  $ljsigeps{$var} = \%blankhash;
      #}
      $count ++;
    }
  }
  $count = 1;
  while ($NUMBER_EXCLUDED_ATOMS =~ /(........)\n?/mg) {
    my $var = trim($1);
    if ($var =~ /^\S+$/) {
      ${$atoms{$count}}{number_excluded_atoms} = $var;
      $count ++;
    }
  }
  $count = 1;
  while ($AMBER_ATOM_TYPE =~ /(....)\n?/mg) {
    my $var = trim($1);
    if ($var =~ /^\S+$/) {
      ${$atoms{$count}}{amber_atom_type} = $var;
      $count ++;
    }
  }
  foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
    ${$atoms{$atom}}{ishydrogen} = 1;
    # assuming that all atoms are hydrogen until we find an entry for them in BONDS_WITHOUT_HYDROGEN
    # this wouldn't pick up non-hydrogen atoms not bound to any other non-hydrogen atoms...
    # e.g. ions... this shouldn't matter as i think we only have to know what's hydrogen and what isn't for bonded interactions
  }
  $count = 1;
  my %bondsbyatom_h;
  while ($BONDS_INC_HYDROGEN =~ /(........)\n?(........)\n?(........)\n?/mg) {
    my $var1 = trim($1); my $var2 = trim($2); my $var3 = trim($3);
    if (($var1 =~ /^\S+$/) && ($var2 =~ /^\S+$/)&& ($var3 =~ /^\S+$/)) {
      ${bonds{$count}}{a1} = $var1 / 3 + 1;
      ${bonds{$count}}{a2} = $var2 / 3 + 1;
      ${bonds{$count}}{type} = $var3;
      $count ++;
      push @{$bondsbyatom_h{$var1 / 3 + 1}}, ($var2 / 3 + 1);
      push @{$bondsbyatom_h{$var2 / 3 + 1}}, ($var1 / 3 + 1);
    }
  }

  # delete triangular H-H bonds in water
  my $delcount;
  foreach my $bond (keys(%bonds)) {
    if (scalar(@{$bondsbyatom_h{$bonds{$bond}{a1}}}) >= 2) {
      if (${$atoms{${$bonds{$bond}}{a1}}}{mass} == ${$atoms{${$bonds{$bond}}{a2}}}{mass}) {
        if (${$atoms{${$bonds{$bond}}{a1}}}{mass} > 5) {
          print "suspicious: deleting what looks a bit like a H-H bond but the atoms involved both have mass of ${$atoms{${$bonds{$bond}}{a1}}}{mass}\n";
        }
        delete($bonds{$bond});
        $delcount ++;
      }
    }
  }
  my %newbonds;
  $count = 1;
  foreach my $bond (sort {$a <=> $b} keys(%bonds)) {
    $newbonds{$count} = $bonds{$bond};
    $count ++;
  }
  %bonds = %newbonds;
  undef(%newbonds);
  if ($delcount) {print "Deleted $delcount H-H bonds (from Amber waters, for cosmetic reasons)\n";}
  
  
  while ($BONDS_WITHOUT_HYDROGEN =~ /(........)\n?(........)\n?(........)\n?/mg) {
    my $var1 = trim($1); my $var2 = trim($2); my $var3 = trim($3);
    if (($var1 =~ /^\S+$/) && ($var2 =~ /^\S+$/)&& ($var3 =~ /^\S+$/)) {
      my $a1 = $var1 / 3  + 1;
      my $a2 = $var2 / 3  + 1;
      ${bonds{$count}}{a1} = $a1;
      ${bonds{$count}}{a2} = $a2;
      ${bonds{$count}}{type} = $var3;
      ${$atoms{$a1}}{ishydrogen} = 0;
      ${$atoms{$a2}}{ishydrogen} = 0;
      $count ++;
    }
  }
  $count = 1;
  while ($BOND_FORCE_CONSTANT =~ /(................)\n?/mg) {
    my $var1 = trim($1);
    if ($var1 =~ /^\S+$/) {
      ${$bond_types{$count}}{force_constant} = $var1;
      $count ++;
    }
  }
  $count = 1;
  while ($BOND_EQUIL_VALUE =~ /\s*(\S+)/mg) {
    ${$bond_types{$count}}{equil_value} = $1;
    $count ++;
  }
  $count = 1;
  while ($ANGLES_INC_HYDROGEN =~ /\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)/mg) {
    ${angles{$count}}{a1} = $1 / 3 + 1;
    ${angles{$count}}{a2} = $2 / 3 + 1;
    ${angles{$count}}{a3} = $3 / 3 + 1;
    ${angles{$count}}{type} = $4;
    $count ++;
  }
  while ($ANGLES_WITHOUT_HYDROGEN =~ /\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)/mg) {
    ${angles{$count}}{a1} = $1 / 3 + 1;
    ${angles{$count}}{a2} = $2 / 3 + 1;
    ${angles{$count}}{a3} = $3 / 3 + 1;
    ${angles{$count}}{type} = $4;
    $count ++;
  }
  $count = 1;
  while ($ANGLE_FORCE_CONSTANT =~ /\s*(\S+)/mg) {
    ${$angle_types{$count}}{force_constant} = $1;
    $count ++;
  }
  $count = 1;
  while ($ANGLE_EQUIL_VALUE =~ /\s*(\S+)/mg) {
    ${$angle_types{$count}}{equil_value} = $1;
    $count ++;
  }
  $count = 1;
  while ($DIHEDRALS_INC_HYDROGEN =~ /\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)/mg) {
    my $a1 = $1; my $a2 = $2; my $a3 = $3; my $a4 = $4;
    if ($a1 >= 0) {${dihedrals{$count}}{a1} = $a1 / 3 + 1} else {${dihedrals{$count}}{a1} = $a1 / 3 - 1};
    if ($a2 >= 0) {${dihedrals{$count}}{a2} = $a2 / 3 + 1} else {${dihedrals{$count}}{a2} = $a2 / 3 - 1};
    if ($a3 >= 0) {${dihedrals{$count}}{a3} = $a3 / 3 + 1} else {${dihedrals{$count}}{a3} = $a3 / 3 - 1};
    if ($a4 >= 0) {${dihedrals{$count}}{a4} = $a4 / 3 + 1} else {${dihedrals{$count}}{a4} = $a4 / 3 - 1};
    ${dihedrals{$count}}{type} = $5;
    $count ++;
  }
  while ($DIHEDRALS_WITHOUT_HYDROGEN =~ /\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)/mg) {
    my $a1 = $1; my $a2 = $2; my $a3 = $3; my $a4 = $4;
    if ($a1 >= 0) {${dihedrals{$count}}{a1} = $a1 / 3 + 1} else {${dihedrals{$count}}{a1} = $a1 / 3 - 1};
    if ($a2 >= 0) {${dihedrals{$count}}{a2} = $a2 / 3 + 1} else {${dihedrals{$count}}{a2} = $a2 / 3 - 1};
    if ($a3 >= 0) {${dihedrals{$count}}{a3} = $a3 / 3 + 1} else {${dihedrals{$count}}{a3} = $a3 / 3 - 1};
    if ($a4 >= 0) {${dihedrals{$count}}{a4} = $a4 / 3 + 1} else {${dihedrals{$count}}{a4} = $a4 / 3 - 1};
    ${dihedrals{$count}}{type} = $5;
    $count ++;
  }
  $count = 1;
  while ($DIHEDRAL_FORCE_CONSTANT =~ /\s*(\S+)/mg) {
    ${$dihedral_types{$count}}{force_constant} = $1;
    $count ++;
  }
  $count = 1;
  while ($DIHEDRAL_PERIODICITY =~ /\s*(\S+)/mg) {
    ${$dihedral_types{$count}}{periodicity} = $1;
    $count ++;
  }
  $count = 1;
  while ($DIHEDRAL_PHASE =~ /\s*(\S+)/mg) {
    ${$dihedral_types{$count}}{phase} = $1;
    $count ++;
  }
  $count = 1;
  while ($LENNARD_JONES_ACOEF =~ /\s*(\S+)/mg) {
    ${$nonbonded_par{$count}}{acoef} = $1;
    $count ++;
  }
  $count = 1;
  while ($LENNARD_JONES_BCOEF =~ /\s*(\S+)/mg) {
    ${$nonbonded_par{$count}}{bcoef} = $1;
    $count ++;
  }
  $count = 1;
  while ($NONBONDED_PARM_INDEX =~ /\s*(\S+)/mg) {
    $nonbonded_index{$count} = $1;
    $count ++;
  }

  my $nljtypes = sqrt(scalar(keys(%nonbonded_index)));
  for (my $i = 1; $i <= $nljtypes; $i++ ) {
    my %blankhash = ();
    $ljsigeps{$i} = \%blankhash;
  }
  foreach my $ljtype (keys(%ljsigeps)) {
    # LJ interaction of each atom with itself for sigma and epsilon
    my $index_into_nonbonded_index = scalar(keys(%ljsigeps)) * ($ljtype-1) + $ljtype;
    my $index_into_nonbonded_par = $nonbonded_index{$index_into_nonbonded_index};
    my $A = ${$nonbonded_par{$index_into_nonbonded_par}}{acoef};
    my $B = ${$nonbonded_par{$index_into_nonbonded_par}}{bcoef};
    my $sigma; my $epsilon;
    if ($A == 0 || $B == 0) {
      $sigma = 1; $epsilon = 0;
    }
    else {
      $sigma = sqrt($A ** (1/3)) / sqrt($B ** (1/3));
      $epsilon = $B*$B / (4*$A);
    }
    #print "ljtype $ljtype A $A B $B sigma $sigma epsilon $epsilon\n";
    ${$ljsigeps{$ljtype}}{sigma} = $sigma;
    ${$ljsigeps{$ljtype}}{epsilon} = $epsilon;
  }


  #add exclusions to the %atoms hash
  foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
    my $excl_count = 0;
    if (!defined(${$atoms{$atom}}{number_excluded_atoms})) {
      die( "number_excluded_atoms for atom $atom not defined\n");
    }
    while ($excl_count < ${$atoms{$atom}}{number_excluded_atoms}) {
      $EXCLUDED_ATOMS_LIST =~ /\s*(\S+)/mg;
      ${${$atoms{$atom}}{excluded}}{$1} = 1;
      $excl_count ++;
    }
  }
  #add residue information
  my @reslabels; $reslabels[0]=0;
  while ($RESIDUE_LABEL =~ /(....)/mg) {
    my $res = $1;
    $res =~ s/\s//g;
    push @reslabels, $res;
  }
  #if (!scalar(@reslabels)) {$reslabels[0]=0;}
  my $rescounter = 0;
  my $previous;
  while ($RESIDUE_POINTER =~ /\s*(\S+)/mg) {
    $rescounter ++;
    my $current = $1;
    if(defined($previous)) {
      for (my $i = $previous; $i < $current; $i ++) {
        ${$atoms{$i}}{resid} = $rescounter - 1;
        ${$atoms{$i}}{resname} = $reslabels[$rescounter-1];
      }
    }
    $previous = $current;
  }
  $rescounter ++;
  for (my $i = $previous; $i <= scalar(keys(%atoms)); $i ++){ 
    ${$atoms{$i}}{resid} = $rescounter - 1;
    ${$atoms{$i}}{resname} = $reslabels[$rescounter - 1];
  }


#  foreach my $atom (sort {$a <=> $b} (keys(%atoms))) {
#    print "$atom: ",join(" ",keys(%{$atoms{$atom}})), "\n";
#  }
#print "\n***\n***\n***\n***\n***\n***\n";
#  foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
#    print $atom, ": ";
#    print (join " ",sort ({$a <=> $b} keys(%{${$atoms{$atom}}{excluded}})));
#    print "\n";
#  }
#  foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
#    if (${$atoms{$atom}}{ishydrogen} == 1) {print "Atom $atom is hydrogen\n";};
#  }
my @return;
push @return, \%atoms, \%bonds, \%angles, \%dihedrals, \%bond_types, \%angle_types, \%dihedral_types, \%nonbonded_par, \%nonbonded_index, \%ljsigeps;
return @return;
#  foreach my $bond (sort {$a <=> $b} keys(%bonds)) {
#    print ${$bonds{$bond}}{a1}; print " ";
#    print ${$bonds{$bond}}{a2}; print " ";
#    print ${$bonds{$bond}}{type}; print "\n";
#  }
  
}

1;
