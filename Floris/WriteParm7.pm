package WriteParm7;
use strict;

use lib "/home/fbuelen/fec";
use ToputilsShared;
use sigtrap;

use Exporter;  
our @ISA = qw(Exporter);
our @EXPORT = qw(WriteParm7);

sub WriteParm7 {
  my $doamberexclusions = 0; # pain in the ass and as far as i know not informative for NAMD or GROMACS
  Cleanup($_[1]);
  print "Writing parm7 format...\n";
  my @collected = @{$_[1]};
  my %atoms = %{$collected[0]}; my %bonds = %{$collected[1]};my %angles = %{$collected[2]};my %dihedrals = %{$collected[3]};my %bond_types = %{$collected[4]};my %angle_types = %{$collected[5]};my %dihedral_types = %{$collected[6]};my %nonbonded_par = %{$collected[7]};my %nonbonded_index = %{$collected[8]};
  my %ljsigeps; if (defined($collected[9] )) {%ljsigeps = %{$collected[9]};}
  my $AB = $_[2];
  my $forgiving = 0;
  # !! just to get WriteParm7 to spit out the new topology without stressing about parameter consistency. Resulting parm7 files should only be for viewing
  if ($forgiving) {
    print "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n! 'forgiving' output mode will write out topologies, but parameters will be wrong!\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  }


  if (%ljsigeps && scalar(keys(%ljsigeps)) ) {
    # rebuild nonbonded A and B based on sigma and epsilon
    # loses a bit on floating point precision with all the big exponents...
  
    # these two are now ignored
    undef(%nonbonded_index);
    undef(%nonbonded_par);
    
    my @ljtypessorted = sort( {$a <=> $b} keys(%ljsigeps));
    my $nonbonded_par_count = 1;
    my $nonbonded_index_count = 1;
    my %reverselookup;
    for (my $i = 0; $i < scalar(@ljtypessorted); $i++) {
      my $type_i = $ljtypessorted[$i];
      my $sig_i = ${$ljsigeps{$type_i}}{sigma};
      my $eps_i = ${$ljsigeps{$type_i}}{epsilon};
      for (my $j = 0; $j < scalar(@ljtypessorted); $j++) {
        my $type_j = $ljtypessorted[$j];
        my $sig_j = ${$ljsigeps{$type_j}}{sigma};
        my $eps_j = ${$ljsigeps{$type_j}}{epsilon};
        my $sig_ij = ($sig_i + $sig_j) / 2;
        my $eps_ij = sqrt($eps_i * $eps_j);
        my $index_into_nonbonded_index = scalar(@ljtypessorted) * ($type_i-1) + $type_j;
        my $index_into_nonbonded_par;
        if ($i <= $j) {
          $index_into_nonbonded_par = $nonbonded_par_count;
          $nonbonded_par_count ++;
          $nonbonded_index{$nonbonded_index_count} = $index_into_nonbonded_par;
          ${$reverselookup{$j}}{$i} = $index_into_nonbonded_par;
        }
        else {
          $index_into_nonbonded_par = ${$reverselookup{$i}}{$j};
        }
        $nonbonded_index{$nonbonded_index_count} = $index_into_nonbonded_par;

        ${$nonbonded_par{$index_into_nonbonded_par}}{acoef} = 4 * ($sig_ij ** 12) * $eps_ij;
        ${$nonbonded_par{$index_into_nonbonded_par}}{bcoef} = 4 * ($sig_ij ** 6) * $eps_ij;
        
        $nonbonded_index_count ++;
      }
    }
  }


  if (!$doamberexclusions) {
    foreach my $atom (keys(%atoms)) {
      ${$atoms{$atom}}{number_excluded_atoms} = 0;
      undef %{${$atoms{$atom}}{excluded}};
    }
  }


  # patch up bonded interactions which don't have a type defined

  foreach my $bond (sort {$a <=> $b} (keys(%bonds))) {
    my $a1 = abs(${$bonds{$bond}}{a1});
    my $a2 = abs(${$bonds{$bond}}{a2});
    if (! $forgiving) {
      if (defined(${$bonds{$bond}}{type})) {
      }
      else {
        if (! (defined(${$bonds{$bond}}{force_constantA}) && defined(${$bonds{$bond}}{equil_valueA}))) {
          print "bond $bond ",join(" ",%{$bonds{$bond}}),"\n";
          die("bond type designation or parameters missing for bond $bond\n"); 
        }
        my $fc; my $eq;
        if ($AB eq "B") {
          $fc = ${$bonds{$bond}}{force_constantB};
          $eq = ${$bonds{$bond}}{equil_valueB};
        }
        else {
          $fc = ${$bonds{$bond}}{force_constantA};
          $eq = ${$bonds{$bond}}{equil_valueA};
        }
        my $maxbt = -1;
        foreach my $bondtype (keys(%bond_types)) {
          if ((${$bond_types{$bondtype}}{force_constant} == $fc) && (${$bond_types{$bondtype}}{equil_value} == $eq)) {
            ${$bonds{$bond}}{type} = $bondtype;
            last;
          }
          if ($bondtype > $maxbt) {$maxbt = $bondtype;}
        }
        if (!defined(${$bonds{$bond}}{type})) {
          ${$bond_types{$maxbt+1}}{force_constant} = $fc;
          ${$bond_types{$maxbt+1}}{equil_value} = $eq;
          ${$bonds{$bond}}{type} = $maxbt+1;
        }
      }
    }
  }     

  foreach my $angle (sort {$a <=> $b} (keys(%angles))) {
    my $a1 = abs(${$angles{$angle}}{a1});
    my $a2 = abs(${$angles{$angle}}{a2});
    my $a3 = abs(${$angles{$angle}}{a3});
    if (! $forgiving) {
      if (defined(${$angles{$angle}}{type})) {
      }
      else {
        if (! (defined(${$angles{$angle}}{force_constantA}) && defined(${$angles{$angle}}{equil_valueA}))) {
          print "angle $angle ",join(" ",%{$angles{$angle}}),"\n";
          die("angle type designation or parameters missing for angle $angle\n"); 
        }
        my $fc; my $eq;
        if ($AB eq "B") {
          $fc = ${$angles{$angle}}{force_constantB};
          $eq = ${$angles{$angle}}{equil_valueB};
        }
        else {
          $fc = ${$angles{$angle}}{force_constantA};
          $eq = ${$angles{$angle}}{equil_valueA};
        }
        my $maxat = -1;
        foreach my $angletype (keys(%angle_types)) {
          if ((${$angle_types{$angletype}}{force_constant} == $fc) && (${$angle_types{$angletype}}{equil_value} == $eq)) {
            ${$angles{$angle}}{type} = $angletype;
            last;
          }
          if ($angletype > $maxat) {$maxat = $angletype;}
        }
        if (!defined(${$angles{$angle}}{type})) {
          ${$angle_types{$maxat+1}}{force_constant} = $fc;
          ${$angle_types{$maxat+1}}{equil_value} = $eq;
          ${$angles{$angle}}{type} = $maxat+1;
        }
      }
    }
  }     

  foreach my $dihedral (sort {$a <=> $b} (keys(%dihedrals))) {   
    my ($a1, $a2, $a3, $a4);
    my $id1 = ${$dihedrals{$dihedral}}{a1}; my $a1flag; if ($id1 < 0) {$a1flag = 1; $a1 = ($id1+1)*3; $id1 *= -1} else {$a1 = ($id1-1)*3};
    my $id2 = ${$dihedrals{$dihedral}}{a2}; my $a2flag; if ($id2 < 0) {$a2flag = 1; $a2 = ($id2+1)*3; $id2 *= -1} else {$a2 = ($id2-1)*3};
    my $id3 = ${$dihedrals{$dihedral}}{a3}; my $a3flag; if ($id3 < 0) {$a3flag = 1; $a3 = ($id3+1)*3; $id3 *= -1} else {$a3 = ($id3-1)*3};
    my $id4 = ${$dihedrals{$dihedral}}{a4}; my $a4flag; if ($id4 < 0) {$a4flag = 1; $a4 = ($id4+1)*3; $id4 *= -1} else {$a4 = ($id4-1)*3};
    if (! $forgiving) {
      if (defined(${$dihedrals{$dihedral}}{type})) {
      }
      else {
        my $fc; my $ph; my $periodicity;
        if (! (defined(${$dihedrals{$dihedral}}{force_constantA}) && defined(${$dihedrals{$dihedral}}{periodicityA}) && defined(${$dihedrals{$dihedral}}{phaseA}))) {
          print "dihedral $dihedral ",join(" ",%{$dihedrals{$dihedral}}),"\n";
          die("dihedral type designation or parameters missing for dihedral $dihedral\n"); 
        }
        if ($AB eq "B") {
          $fc = ${$dihedrals{$dihedral}}{force_constantB};
          $ph = ${$dihedrals{$dihedral}}{phaseB};
          $periodicity = ${$dihedrals{$dihedral}}{periodicityB};
        }
        else {
          $fc = ${$dihedrals{$dihedral}}{force_constantA};
          $ph = ${$dihedrals{$dihedral}}{phaseA};
          $periodicity = ${$dihedrals{$dihedral}}{periodicityA};
        }
        my $maxdt = -1;
        foreach my $dihedraltype (keys(%dihedral_types)) {
          if ((${$dihedral_types{$dihedraltype}}{force_constant} == $fc) && (${$dihedral_types{$dihedraltype}}{phase} == $ph) && (${$dihedral_types{$dihedraltype}}{periodicity} == $periodicity)) {
            ${$dihedrals{$dihedral}}{type} = $dihedraltype;
            last;
          }
          if ($dihedraltype > $maxdt) {$maxdt = $dihedraltype;}
        }
        if (!defined(${$dihedrals{$dihedral}}{type})) {
          ${$dihedral_types{$maxdt+1}}{force_constant} = $fc;
          ${$dihedral_types{$maxdt+1}}{phase} = $ph;
          ${$dihedral_types{$maxdt+1}}{periodicity} = $periodicity;
          ${$dihedrals{$dihedral}}{type} = $maxdt+1;
        }
      }
    }
  }     
  

  open OUTFILE, (">$_[0]") or die("couldn't open outfile $_[0]\n");
  printf OUTFILE ("%-80s\n", "\%VERSION  VERSION_STAMP = V0001.000  DATE = ");
  printf OUTFILE ("%-80s\n", "\%FLAG TITLE");
  printf OUTFILE ("%-80s\n", "\%FORMAT(20a4)");
  printf OUTFILE ("%-80s\n", "WRITEPARM7.PM");


  
  my @pointers;
  for (my $i = 0; $i <= 30; $i++) {$pointers[$i] =0;}
#  POINTERS SECTION
#  0  NATOM  : total number of atoms 
  $pointers[0] = scalar(keys(%atoms));
#  1  NTYPES : total number of distinct atom types
  my %temp;
  if ($AB eq "B") {
    foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
      if (defined($temp{${atoms{$atom}}{atom_type_indexB}})) {
        $temp{${atoms{$atom}}{atom_type_indexB}} = 1;
      }
      else {
        $temp{${atoms{$atom}}{atom_type_index}} = 1;
      }
    }
  }
  else {
    foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
      $temp{${atoms{$atom}}{atom_type_index}} = 1;
    }
  }
  $pointers[1] = scalar(keys(%temp));
#  2  NBONH  : number of bonds containing hydrogen
#  3  MBONA  : number of bonds not containing hydrogen
  foreach my $bond (keys(%bonds)) {
    my $a1 = ${$bonds{$bond}}{a1};
    my $a2 = ${$bonds{$bond}}{a2};
    if (${atoms{$a1}}{ishydrogen} || ${$atoms{$a2}}{ishydrogen}) {$pointers[2] ++;}
    else {$pointers[3] ++;}
  }
#  4  NTHETH : number of angles containing hydrogen
#  5  MTHETA : number of angles not containing hydrogen
  foreach my $angle (keys(%angles)) {
    my $a1 = abs(${$angles{$angle}}{a1});
    my $a2 = abs(${$angles{$angle}}{a2});
    my $a3 = abs(${$angles{$angle}}{a3});
    if (${$atoms{$a1}}{ishydrogen} || ${$atoms{$a2}}{ishydrogen} || ${$atoms{$a3}}{ishydrogen}) {$pointers[4] ++;}
    else {$pointers[5] ++;}
  }
#  6  NPHIH  : number of dihedrals containing hydrogen
#  7  MPHIA  : number of dihedrals not containing hydrogen
  foreach my $dihedral (keys(%dihedrals)) {
    my $a1 = abs(${$dihedrals{$dihedral}}{a1});
    my $a2 = abs(${$dihedrals{$dihedral}}{a2});
    my $a3 = abs(${$dihedrals{$dihedral}}{a3});
    my $a4 = abs(${$dihedrals{$dihedral}}{a4});
    if (${atoms{$a1}}{ishydrogen} || ${$atoms{$a2}}{ishydrogen} || ${$atoms{$a3}}{ishydrogen} || ${$atoms{$a4}}{ishydrogen}) {$pointers[6] ++;}
    else {$pointers[7] ++;}
  }
#  8  NHPARM : currently not used
  $pointers[8] = 0;
#  9  NPARM  : currently not used
  $pointers[9] = 0;
#  10 NEXT   : number of excluded atoms
#  11 NRES   : number of residues
  my $nres = 0; my $prev = "";
  foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
    $pointers[10] += scalar(keys(  %{${$atoms{$atom}}{excluded}} ));
    if (! (${$atoms{$atom}}{resid} eq $prev)) {
      $nres ++;
      $prev = ${$atoms{$atom}}{resid};
    }
  }
  #$pointers[11] = ${$atoms{scalar(keys(%atoms))}}{resid};
  $pointers[11] = $nres;
#  12 NBONA  : MBONA + number of constraint bonds
  $pointers[12] = $pointers[3];
#  13 NTHETA : MTHETA + number of constraint angles
  $pointers[13] = $pointers[5];
#  14 NPHIA  : MPHIA + number of constraint dihedrals
  $pointers[14] = $pointers[7];
#  15 NUMBND : number of unique bond types
  $pointers[15] = scalar(keys(%bond_types));
#  16 NUMANG : number of unique angle types
  $pointers[16] = scalar(keys(%angle_types));
#  17 NPTRA  : number of unique dihedral types
  $pointers[17] = scalar(keys(%dihedral_types));
#  18 NATYP  : number of atom types in parameter file, see SOLTY below
  $pointers[18] = $pointers[1];
#  19 NPHB   : number of distinct 10-12 hydrogen bond pair types
  $pointers[19] = 1;
#  20 IFPERT : set to 1 if perturbation info is to be read in
  $pointers[20] = 0;
#  21 NBPER  : number of bonds to be perturbed
  $pointers[21] = 0;
#  22 NGPER  : number of angles to be perturbed
  $pointers[22] = 0;
#  23 NDPER  : number of dihedrals to be perturbed
  $pointers[23] = 0;
#  24 MBPER  : number of bonds with atoms completely in perturbed group
  $pointers[24] = 0;
#  25 MGPER  : number of angles with atoms completely in perturbed group
  $pointers[25] = 0;
#  26 MDPER  : number of dihedrals with atoms completely in perturbed groups
  $pointers[26] = 0;
#  27 IFBOX  : set to 1 if standard periodic box, 2 when truncated octahedral
  $pointers[27] = 1;
#  28 NMXRS  : number of atoms in the largest residue
  my $biggest = 0;
  my $previous = 0;
  my $current_count = 0;
  foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
    my $resid = ${$atoms{$atom}}{resid};
    if (!($resid eq $previous)) {
      if ($current_count > $biggest) {$biggest = $current_count;}
      $current_count = 0;
      $previous = $resid;
    }
    $current_count ++;
  }
    
    
  $pointers[28] = $biggest;
#  29 IFCAP  : set to 1 if the CAP option from edit was specified
  $pointers[29] = 0;
  $pointers[30] = 0;

  if (0) {
    print "Pointers:\n";
    for (my $i = 0; $i < scalar (@pointers); $i++) {
      print "[$i] $pointers[$i]\n";
    }
  }

  foreach my $val (@pointers) {if (!defined($val)) {die("missing value in array  POINTERS\n");}}
  print OUTFILE output("POINTERS", "10I8",\@pointers);
  my @temparray;
  if ($AB eq "B") {
    foreach my $atom (sort {$a <=> $b} (keys(%atoms))) {
      if (defined(${$atoms{$atom}}{nameB})) {
        push @temparray, ${$atoms{$atom}}{nameB};
      }
      else {
        push @temparray, ${$atoms{$atom}}{name};
      }
    }
  }
  else {
    foreach my $atom (sort {$a <=> $b} (keys(%atoms))) {
      push @temparray, ${$atoms{$atom}}{name};
    }
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  ATOM_NAME\n");}}
  print OUTFILE output("ATOM_NAME", "20a4",\@temparray);
  undef(@temparray);
  if ($AB eq "B") {
    foreach my $atom (sort {$a <=> $b} (keys(%atoms))) {
      if (defined(${$atoms{$atom}}{chargeB})) {
        push @temparray, ${$atoms{$atom}}{chargeB} * 18.2223;
      }
      else {
        push @temparray, ${$atoms{$atom}}{charge} * 18.2223;
      }
    }
  }
  else {
    foreach my $atom (sort {$a <=> $b} (keys(%atoms))) {
      push @temparray, ${$atoms{$atom}}{charge} * 18.2223;
    }
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  CHARGE\n");}}
  print OUTFILE output("CHARGE", "5E16.8",\@temparray);
  undef(@temparray);
  if ($AB eq "B") {
    foreach my $atom (sort {$a <=> $b} (keys(%atoms))) {
      if (defined(${$atoms{$atom}}{massB})) {
        push @temparray, ${$atoms{$atom}}{massB};
      }
      else {
        push @temparray, ${$atoms{$atom}}{mass};
      }
    }
  }
  else {
    foreach my $atom (sort {$a <=> $b} (keys(%atoms))) {
      push @temparray, ${$atoms{$atom}}{mass};
    }
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  MASS\n");}}
  print OUTFILE output("MASS", "5E16.8",\@temparray);
  undef(@temparray);
  if ($AB eq "B") {
    foreach my $atom (sort {$a <=> $b} (keys(%atoms))) {
      if (defined(${$atoms{$atom}}{atom_type_indexB})) {
        push @temparray, ${$atoms{$atom}}{atom_type_indexB};
      }
      else {
        push @temparray, ${$atoms{$atom}}{atom_type_index};
      }
    }
  }
  else {
    foreach my $atom (sort {$a <=> $b} (keys(%atoms))) {
      push @temparray, ${$atoms{$atom}}{atom_type_index};
    }
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  ATOM_TYPE_INDEX\n");}}
  print OUTFILE output("ATOM_TYPE_INDEX", "10I8",\@temparray);
  undef(@temparray);
  foreach my $atom (sort {$a <=> $b} (keys(%atoms))) {
    if (${$atoms{$atom}}{number_excluded_atoms} != scalar(keys(%{${$atoms{$atom}}{excluded}})) ){
      if (! $forgiving) {
        die ("number_excluded_atoms value in %atoms hash does not agree with the number of entries in sub-hash 'excluded' for atom $atom\n");
        #print "number_excluded_atoms value in %atoms hash does not agree with the number of entries in sub-hash 'excluded' for atom $atom\n";
        #${$atoms{$atom}}{number_excluded_atoms} = scalar(keys(%{${$atoms{$atom}}{excluded}}));
      }
    }
    push @temparray, ${$atoms{$atom}}{number_excluded_atoms};
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  NUMBER_EXCLUDED_ATOMS\n");}}
  print OUTFILE output("NUMBER_EXCLUDED_ATOMS", "10I8",\@temparray);
  undef(@temparray);
  foreach my $val (sort {$a <=> $b} keys(%nonbonded_index)) {
    push @temparray, $nonbonded_index{$val};
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  NONBONDED_PARM_INDEX\n");}}
  print OUTFILE output("NONBONDED_PARM_INDEX", "10I8",\@temparray);
  undef(@temparray);

  my @temparray2;
  $previous = 0;
  foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
    my $resid = ${$atoms{$atom}}{resid};
    if (!($resid eq $previous)) {
      $previous = $resid;
      if ($AB eq "B") {
        if (defined(${$atoms{$atom}}{resnameB})) {
          push @temparray,  ${$atoms{$atom}}{resnameB};
        }
        else {
          push @temparray,  ${$atoms{$atom}}{resname};
        }
      }
      else {
        push @temparray,  ${$atoms{$atom}}{resname};
      }
      push @temparray2,  $atom;    
    }
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  RESIDUE_LABEL\n");}}
  print OUTFILE output("RESIDUE_LABEL", "20a4",\@temparray);
  foreach my $val (@temparray2) {if (!defined($val)) {die("missing value in array  RESIDUE_POINTER\n");}}
  print OUTFILE output("RESIDUE_POINTER", "10I8",\@temparray2);
  undef(@temparray); undef(@temparray2);

#  print OUTFILE "%FLAG RESIDUE_LABEL                                                             \n%FORMAT(20a4)                                                                   \naaaaaaaa\n";
#  print OUTFILE "%FLAG RESIDUE_POINTER                                                           \n%FORMAT(10I8)                                                                   \naaaaaaaa\n";

  foreach my $bondtype (sort {$a <=> $b} keys(%bond_types)) {
    push @temparray, ${$bond_types{$bondtype}}{force_constant};
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  BOND_FORCE_CONSTANT\n");}}
  print OUTFILE output("BOND_FORCE_CONSTANT", "5E16.8",\@temparray);
  undef(@temparray);
  foreach my $bondtype (sort {$a <=> $b} keys(%bond_types)) {
    push @temparray, ${$bond_types{$bondtype}}{equil_value};
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  BOND_EQUIL_VALUE\n");}}
  print OUTFILE output("BOND_EQUIL_VALUE", "5E16.8",\@temparray);
  undef(@temparray);

  foreach my $angletype (sort {$a <=> $b} keys(%angle_types)) {
    push @temparray, ${$angle_types{$angletype}}{force_constant};
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  ANGLE_FORCE_CONSTANT\n");}}
  print OUTFILE output("ANGLE_FORCE_CONSTANT", "5E16.8",\@temparray);
  undef(@temparray);
  foreach my $angletype (sort {$a <=> $b} keys(%angle_types)) {
    push @temparray, ${$angle_types{$angletype}}{equil_value};
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  ANGLE_EQUIL_VALUE\n");}}
  print OUTFILE output("ANGLE_EQUIL_VALUE", "5E16.8",\@temparray);
  undef(@temparray);

  foreach my $dihedraltype (sort {$a <=> $b} keys(%dihedral_types)) {
    push @temparray, ${$dihedral_types{$dihedraltype}}{force_constant};
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  DIHEDRAL_FORCE_CONSTANT\n");}}
  print OUTFILE output("DIHEDRAL_FORCE_CONSTANT", "5E16.8",\@temparray);
  undef(@temparray);
  foreach my $dihedraltype (sort {$a <=> $b} keys(%dihedral_types)) {
    push @temparray, ${$dihedral_types{$dihedraltype}}{periodicity};
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  DIHEDRAL_PERIODICITY\n");}}
  print OUTFILE output("DIHEDRAL_PERIODICITY", "5E16.8",\@temparray);
  undef(@temparray);
  foreach my $dihedraltype (sort {$a <=> $b} keys(%dihedral_types)) {
    push @temparray, ${$dihedral_types{$dihedraltype}}{phase};
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  DIHEDRAL_PHASE\n");}}
  print OUTFILE output("DIHEDRAL_PHASE", "5E16.8",\@temparray);
  undef(@temparray);

  for (my $i = 0; $i < $pointers[18]; $i ++) {
    push @temparray, "  0.00000000E+00";
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  SOLTY\n");}}
  print OUTFILE output("SOLTY", "5E16.8",\@temparray);
  undef(@temparray);

  foreach my $val (sort {$a <=> $b} keys(%nonbonded_par)) {
    push @temparray, ${$nonbonded_par{$val}}{acoef};
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  LENNARD_JONES_ACOEF\n");}}
  print OUTFILE output("LENNARD_JONES_ACOEF", "5E16.8",\@temparray);
  undef(@temparray);
  foreach my $val (sort {$a <=> $b} keys(%nonbonded_par)) {
    push @temparray, ${$nonbonded_par{$val}}{bcoef};
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  LENNARD_JONES_BCOEF\n");}}
  print OUTFILE output("LENNARD_JONES_BCOEF", "5E16.8",\@temparray);
  undef(@temparray);

  foreach my $bond (sort {$a <=> $b} (keys(%bonds))) {
    my $a1 = abs(${$bonds{$bond}}{a1});
    my $a2 = abs(${$bonds{$bond}}{a2});
    my $type;
    if (! $forgiving) {
      if (defined(${$bonds{$bond}}{type})) {
        $type = abs(${$bonds{$bond}}{type});
      }
      else {
        die("no bond type entry\n");
      }
    }
    else {
      $type = -1;
    }
    if (!defined($bond_types{$type})) {
      die("bond $bond type $type not present in \%bond_types\n");
    }
    if (${$atoms{$a1}}{ishydrogen} || ${$atoms{$a2}}{ishydrogen} ) {
      push @temparray, ($a1-1)*3, ($a2-1)*3, $type;
    }
    else {
      push @temparray2, ($a1-1)*3, ($a2-1)*3, $type;
    }
  }     

  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  BONDS_INC_HYDROGEN\n");}}
  print OUTFILE output("BONDS_INC_HYDROGEN", "10I8",\@temparray);
  foreach my $val (@temparray2) {if (!defined($val)) {die("missing value in array  BONDS_WITHOUT_HYDROGEN\n");}}
  print OUTFILE output("BONDS_WITHOUT_HYDROGEN", "10I8",\@temparray2);
  undef(@temparray);undef(@temparray2);

  foreach my $angle (sort {$a <=> $b} (keys(%angles))) {
    my $a1 = abs(${$angles{$angle}}{a1});
    my $a2 = abs(${$angles{$angle}}{a2});
    my $a3 = abs(${$angles{$angle}}{a3});
    my $type;
    if (! $forgiving) {
      if (defined(${$angles{$angle}}{type})) {
        $type = abs(${$angles{$angle}}{type});
      }
      else {
        die("no angle type entry\n");
      }
    }
    else {
      $type = -1;
    }
    if (!defined($angle_types{$type})) {
      die("angle $angle type $type not present in \%angle_types\n");
    }
    if (${$atoms{$a1}}{ishydrogen} || ${$atoms{$a2}}{ishydrogen} || ${$atoms{$a3}}{ishydrogen}) {
      push @temparray, ($a1-1)*3, ($a2-1)*3, ($a3-1)*3, $type;
    }
    else {
      push @temparray2, ($a1-1)*3, ($a2-1)*3, ($a3-1)*3, $type;
    }
  }     



  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  ANGLES_INC_HYDROGEN\n");}}
  print OUTFILE output("ANGLES_INC_HYDROGEN", "10I8",\@temparray);
  foreach my $val (@temparray2) {if (!defined($val)) {die("missing value in array  ANGLES_WITHOUT_HYDROGEN\n");}}
  print OUTFILE output("ANGLES_WITHOUT_HYDROGEN", "10I8",\@temparray2);
  undef(@temparray);undef(@temparray2);

  foreach my $dihedral (sort {$a <=> $b} (keys(%dihedrals))) {   
    my ($a1, $a2, $a3, $a4);
    my $id1 = ${$dihedrals{$dihedral}}{a1}; my $a1flag; if ($id1 < 0) {$a1flag = 1; $a1 = ($id1+1)*3; $id1 *= -1} else {$a1 = ($id1-1)*3};
    my $id2 = ${$dihedrals{$dihedral}}{a2}; my $a2flag; if ($id2 < 0) {$a2flag = 1; $a2 = ($id2+1)*3; $id2 *= -1} else {$a2 = ($id2-1)*3};
    my $id3 = ${$dihedrals{$dihedral}}{a3}; my $a3flag; if ($id3 < 0) {$a3flag = 1; $a3 = ($id3+1)*3; $id3 *= -1} else {$a3 = ($id3-1)*3};
    my $id4 = ${$dihedrals{$dihedral}}{a4}; my $a4flag; if ($id4 < 0) {$a4flag = 1; $a4 = ($id4+1)*3; $id4 *= -1} else {$a4 = ($id4-1)*3};
    my $type;
    if (! $forgiving) {
      if (defined(${$dihedrals{$dihedral}}{type})) {
        $type = abs(${$dihedrals{$dihedral}}{type});
      }
      else {
        die("no dihedral type entry\n");
      }
    }
    else {
      $type = -1;
    }
    if (!defined($dihedral_types{$type})) {
      die("dihedral $dihedral type $type not present in \%dihedral_types\n");
    }
    if (${$atoms{$id1}}{ishydrogen} || ${$atoms{$id2}}{ishydrogen} || ${$atoms{$id3}}{ishydrogen} || ${$atoms{$id4}}{ishydrogen}) {
      push @temparray, $a1, $a2, $a3, $a4, $type;
    }
    else {
      push @temparray2,$a1, $a2, $a3, $a4, $type;
    }
  }     


  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  DIHEDRALS_INC_HYDROGEN\n");}}
  print OUTFILE output("DIHEDRALS_INC_HYDROGEN", "10I8",\@temparray);
  foreach my $val (@temparray2) {if (!defined($val)) {die("missing value in array  DIHEDRALS_WITHOUT_HYDROGEN\n");}}
  print OUTFILE output("DIHEDRALS_WITHOUT_HYDROGEN", "10I8",\@temparray2);
  undef(@temparray);undef(@temparray2);

  if ($doamberexclusions) {
    foreach my $atom (sort {$a <=> $b} (keys(%atoms))) {
      foreach my $excl (sort {$a <=> $b} keys(%{${$atoms{$atom}}{excluded}})) {
        push @temparray, $excl;
      }
    }
    foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  EXCLUDED_ATOMS_LIST\n");}}
    print OUTFILE output("EXCLUDED_ATOMS_LIST", "10I8",\@temparray);
    undef(@temparray);
  }
  else {
    print OUTFILE "%FLAG EXCLUDED_ATOMS_LIST                                                       \n%FORMAT(10I8)                                                                   \n       0       0       0       0       0       0       0       0       0       0\n";
  }


  print OUTFILE "%FLAG HBOND_ACOEF                                                               \n%FORMAT(5E16.8)                                                                 \n  0.00000000E+00\n%FLAG HBOND_BCOEF                                                               \n%FORMAT(5E16.8)                                                                 \n  0.00000000E+00\n%FLAG HBCUT                                                                     \n%FORMAT(5E16.8)                                                                 \n  0.00000000E+00\n";
  foreach my $atom (sort {$a <=> $b} (keys(%atoms))) {
    if ($AB eq "B") {
      if (defined(${$atoms{$atom}}{amber_atom_typeB})) {
        push @temparray, ${$atoms{$atom}}{amber_atom_typeB};
      }
      else {
        push @temparray, ${$atoms{$atom}}{amber_atom_type};
      }
    }
    else {
      push @temparray, ${$atoms{$atom}}{amber_atom_type};
    }
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  AMBER_ATOM_TYPE\n");}}
  print OUTFILE output("AMBER_ATOM_TYPE", "20a4",\@temparray);
  undef(@temparray);
  foreach my $atom (sort {$a <=> $b} (keys(%atoms))) {
    push @temparray, "BLA";
  }
  foreach my $val (@temparray) {if (!defined($val)) {die("missing value in array  TREE_CHAIN_CLASSIFICATION\n");}}
  print OUTFILE output("TREE_CHAIN_CLASSIFICATION", "20a4",\@temparray);

  print OUTFILE "%FLAG JOIN_ARRAY                                                                \n%FORMAT(10I8)                                                                   \n       0       0       0       0       0       0       0       0       0       0\n%FLAG IROTAT                                                                    \n%FORMAT(10I8)                                                                   \n       0       0       0       0       0       0       0       0       0       0\n%FLAG SOLVENT_POINTERS                                                          \n%FORMAT(3I8)                                                                    \n     111    1111     111\n%FLAG ATOMS_PER_MOLECULE                                                        \n%FORMAT(10I8)                                                                   \n       1       1       1       1       1       1       1       1       1       1\n%FLAG BOX_DIMENSIONS                                                            \n%FORMAT(5E16.8)                                                                 \n  9.00000000E+01  7.64092210E+01  8.26312710E+01  8.25145990E+01\n%FLAG RADIUS_SET                                                                \n%FORMAT(1a80)                                                                   \nmodified Bondi radii (mbondi)                                                   \n%FLAG RADII                                                                     \n%FORMAT(5E16.8)                                                                 \n  1.70000000E+00  1.70000000E+00  1.70000000E+00  1.70000000E+00  1.70000000E+00\n%FLAG SCREEN                                                                    \n%FORMAT(5E16.8)                                                                 \n  7.20000000E-01  7.20000000E-01  7.20000000E-01  7.20000000E-01  7.20000000E-01\n";
  close(OUTFILE);
}
sub output {   #1) FLAG, 2) output type, 3) array reference with output data
  my $out = sprintf("%-80s\n","\%FLAG $_[0]");
  $out .= sprintf("%-80s\n","\%FORMAT($_[1])");
  my @data = @{$_[2]};
  if (!scalar(@data)) {
    print "$_[0] (!!!!---empty---!!!!)     ";
    $out .= "\n";
    return $out;
  }
  #else {
    #print "$_[0] (",scalar(@data),")     ";
  #}
  if ($_[1] eq "10I8") {
    my $linecount = 0;
    foreach my $value (@data) {
      if ($linecount == 10) {$out .= "\n"; $linecount = 0;}
      $out .= sprintf("%8i",$value);
      $linecount ++;
    }
  }
  if ($_[1] eq "20a4") {
    my $linecount = 0;
    foreach my $value (@data) {
      if ($linecount == 20) {$out .= "\n"; $linecount = 0;}
      $out .= sprintf("%-4s",$value);
      $linecount ++;
    }
  }
  if ($_[1] eq "5E16.8") {
    my $linecount = 0;
    foreach my $value (@data) {
      if ($linecount == 5) {$out .= "\n"; $linecount = 0;}
      $out .= sprintf("%16.8E",$value);
      $linecount ++;
    }
  }
  $out .= "\n";
  return $out;
}
1;
