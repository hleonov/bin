package WriteGTop;
use strict;

use lib "/home/fbuelen/fec";
use ToputilsShared;

use Exporter;  
our @ISA = qw(Exporter);
our @EXPORT = qw(WriteGTop);

no warnings 'recursion';

my $caltoj = 4.184;
my $degtorad = 57.2957795;
my $newformatPosres = 0;  # if true, write wellWidth and alchGroup columns in position restraints
my $greedymol = 1; # lump everything that's not solvent and ion together into one molecule

sub WriteGTop {
  Cleanup($_[1]);
  print "Writing gromacs topology format...\n";
  my @collected = @{$_[1]};
  my %atoms = %{$collected[0]};
  my %bonds = %{$collected[1]};
  my %angles = %{$collected[2]};
  my %dihedrals = %{$collected[3]};
  my %bond_types = %{$collected[4]};
  my %angle_types = %{$collected[5]};
  my %dihedral_types = %{$collected[6]};
  my %nonbonded_par = %{$collected[7]};
  my %nonbonded_index = %{$collected[8]};
  my %ljsigeps; if (defined( $collected[9] )) {%ljsigeps = %{$collected[9]};} else {die("No lj sigma and epsilon, this shouldn't be...\n");}
  my $parmfilewas = $collected[10];
  my %pairs; if (defined( $collected[11] )) {%pairs = %{$collected[11]};} 


  # check if we have any A / B state information for free energy
  my $perturbedNonbonded = 0;
  my $perturbedBonds = 0;
  my $perturbedAngles = 0;
  my $perturbedDihedrals = 0;
  foreach my $atom (keys(%atoms)) {
    if (defined(${$atoms{$atom}}{chargeB})) {
      $perturbedNonbonded = 1;
      last;
    }
  }
  foreach my $bond (keys(%bonds)) {
    if (defined(${$bonds{$bond}}{force_constantB})) {
      $perturbedBonds = 1;
      last;
    }
  }
  foreach my $angle (keys(%angles)) {
    if (defined(${$angles{$angle}}{force_constantB})) {
      $perturbedAngles = 1;
      last;
    }
  }
  foreach my $dihedral (keys(%dihedrals)) {
    if (defined(${$dihedrals{$dihedral}}{force_constantB})) {
      $perturbedDihedrals = 1;
      last;
    }
  }
  if ($perturbedNonbonded) {print "Alchemical A and B state information found for nonbonded interactions\n";}
  if ($perturbedBonds) {print "Alchemical A and B state information found for bonds\n";}
  if ($perturbedAngles) {print "Alchemical A and B state information found for angles\n";}
  if ($perturbedDihedrals) {print "Alchemical A and B state information found for dihedrals\n";}


  # try to identify molecules, including repeating ones, as generally as possible...?

  my @molecules;
  my %bondsbyatom; # this gets dismantled later...
  my %bondsbyatom2; # I want a second copy for group assignment later
  foreach my $bond (keys(%bonds)) {
    ${$bondsbyatom{${$bonds{$bond}}{a1}}}{${$bonds{$bond}}{a2}} = 1;
    ${$bondsbyatom{${$bonds{$bond}}{a2}}}{${$bonds{$bond}}{a1}} = 1;
    ${$bondsbyatom2{${$bonds{$bond}}{a1}}}{${$bonds{$bond}}{a2}} = 1;
    ${$bondsbyatom2{${$bonds{$bond}}{a2}}}{${$bonds{$bond}}{a1}} = 1;
  }

  #find single atoms (no bonds)
  foreach my $atom (keys(%atoms)) {
    if (! defined($bondsbyatom{$atom})) {
      ${$bondsbyatom{$atom}}{-1} = 1;
    }
  }

  my $startatom = 1; 
  my $lastperc = 0;
  while (scalar(keys(%bondsbyatom))) {
    my %newmol;
    # my $perc = 100 - int(scalar(keys(%bondsbyatom)) / scalar(keys(%atoms)) * 100);
    # if ($perc > $lastperc) {printf ("\rIdentifying molecules...%3i%%",$perc); $lastperc = $perc;}
    AddAtomsToMol(\%bondsbyatom, \@molecules, $startatom, \%newmol);
    $startatom ++;
    if (scalar(keys(%newmol))) {push @molecules, \%newmol;}
  }

  print "\nFound ",scalar(@molecules), " molecules\n";
  if (!scalar(@molecules)) {die("No molecules found\n");}
  my %molsizes;
  foreach my $molecule (@molecules) {$molsizes{scalar(keys(%{$molecule}))} ++;}
  my $checksum;
  print "Fragments: ";
  my @fragments;
  foreach my $size (sort ({$a <=> $b} keys(%molsizes))) {
    push @fragments, "$molsizes{$size} of $size";
    $checksum += $size * $molsizes{$size};
  }
  print join(", ",@fragments),"\n";
  if (scalar(keys(%atoms)) != $checksum) {
    print "checksum ",scalar(keys(%atoms)), " $checksum\n";
    die("");
  }

  # find repeating units among these molecules...
  my %moleculefingerprints;
  my %moleculenumber_atomcount;
  my $fpcount = 0; my %fpenum;
  my $molnumbercount = 0;
  foreach my $molecule (@molecules) {
    my @molatoms = sort ({$a <=> $b} keys(%{$molecule}));
    my $key; 
    $molnumbercount++;
    # concatenate a load of information to use as a hash key
    # includes: atom types, charges, masses, exclusions (numbered relative to 1st atom) 
    # all in a defined order. Should give a unique fingerprint for all practical purposes as far as I'm concerned...
    # ...but bonded parameters aren't checked, so in theory two otherwise identical molecules with same connectivity
    # but different bonded parameters could end up incorrectly lumped together...
    # correct this one day if full generality is desired...

    # now including FEP and restraint flags from NAMD fep and constr.pdb
    # which means that identical molecules with differing lambda coupling / restraints will be handled separately
    
    my $lowestid = $molatoms[0];
    foreach my $molatom (sort ({$a <=> $b} @molatoms)) {
      ${$atoms{$molatom}}{moleculenumber} = $molnumbercount;
      $moleculenumber_atomcount{$molnumbercount} ++;
      #if (defined(${$atoms{$molatom}}{fepflag}) && defined(${$atoms{$molatom}}{restraint_k}) && defined(${$atoms{$molatom}}{restraint_fbhw})) {
      #  $key .= ${$atoms{$molatom}}{atom_type_index}.${$atoms{$molatom}}{charge}.${$atoms{$molatom}}{ishydrogen}.${$atoms{$molatom}}{mass}.${$atoms{$molatom}}{fepflag}.${$atoms{$molatom}}{restraint_k}.${$atoms{$molatom}}{restraint_fbhw};
      #}
      #else {
      #  $key .= ${$atoms{$molatom}}{atom_type_index}.${$atoms{$molatom}}{charge}.${$atoms{$molatom}}{ishydrogen}.${$atoms{$molatom}}{mass};
      #}
      
      $key .= ${$atoms{$molatom}}{amber_atom_type}.${$atoms{$molatom}}{atom_type_index}.":".${$atoms{$molatom}}{charge}.":".${$atoms{$molatom}}{ishydrogen}.":".${$atoms{$molatom}}{mass};
      if (defined(${$atoms{$molatom}}{fepflag})) {$key .= ":".${$atoms{$molatom}}{fepflag};}
      if (defined(${$atoms{$molatom}}{restraint_k})) {$key .= ":".${$atoms{$molatom}}{restraint_k};}
      if (defined(${$atoms{$molatom}}{restraint_fbhw})) {$key .= ":".${$atoms{$molatom}}{restraint_fbhw};}
      if (defined(${$atoms{$molatom}}{atom_type_indexB})) {$key .= ":".${$atoms{$molatom}}{atom_type_indexB};}
      if (defined(${$atoms{$molatom}}{chargeB})) {$key .= ":".${$atoms{$molatom}}{chargeB};}
      if (defined(${$atoms{$molatom}}{massB})) {$key .= ":".${$atoms{$molatom}}{massB};}
      foreach my $exclatom (sort({$a <=> $b} keys(%{${$atoms{$molatom}}{excluded}}))) {
        if ($exclatom) { $key .= sprintf(":*%i*",$exclatom-$lowestid); }
      }
    }
    if (!defined ($moleculefingerprints{$key})) {$fpcount ++; $fpenum{$key} = $fpcount;}
    ${$moleculefingerprints{$key}}{count} ++;
    if (! defined(${$moleculefingerprints{$key}}{atoms}) ) {
      foreach my $molatom (@molatoms) {
        push @{${$moleculefingerprints{$key}}{atoms}}, $molatom;
      }
    }
    foreach my $molatom (@molatoms) {
      ${$atoms{$molatom}}{moleculeid} = $fpenum{$key};
    }
  }

  my %molecules_representedonce;
  foreach my $fingerprint (keys(%moleculefingerprints)) {
    if (${$moleculefingerprints{$fingerprint}}{count} == 1) {
      $molecules_representedonce{$fpenum{$fingerprint}} = 1;
    }
  }

  # a 'group' in this context is a contiguous collection of unique molecules
  # ...but with different lambda-coupled groups separate
  # e.g. LG1 and LG2 following each other, lambda-coupled in different directions, 
  # are classified 'unique' because of FEP flags, but need to be two separate groups
  
  my %moleculenumbertogroup;
  my %moleculeidtogroup;
  my $groupid = 0;
  my $lastgroupid = 0;
  my $merge = 0; # consecutive unique molecules packaged into one group
  my $saveAtomsToRepresentGroup = 0;
  my %OneMoleculeRepresentsEachGroup;
  foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
    if ($atom > 1) {
      if ($merge == 0 && (${$atoms{$atom-1}}{moleculenumber} != ${$atoms{$atom}}{moleculenumber}) && defined( $OneMoleculeRepresentsEachGroup{$groupid} )) {
        $saveAtomsToRepresentGroup = 0;
      }
      if (defined(${$atoms{$atom-1}}{fepflag}) && defined (${$atoms{$atom}}{fepflag})) {
        if (${$atoms{$atom-1}}{fepflag} != ${$atoms{$atom}}{fepflag}) {$merge = 0;}
      }
    }

    # moleculenumber: every contiguous segment of atoms, bonded together, has its moleculenumber
    # moleculeid: same, except if there were identical segments earlier they have the same moleculeid
    
    if (! defined( $moleculenumbertogroup{${$atoms{$atom}}{moleculenumber}} )) {
      if (defined ($molecules_representedonce{${$atoms{$atom}}{moleculeid}})) {
        # any unique molecule (not bonded to rest of structure) gets merged in with previous
        $merge++;
      }
      elsif (${$atoms{$atom}}{resname} eq "ZN") {  # hack to force ZNs into previous group, for distance restraints
        print "!!!HACK\n!!!HACK force ZN into previous group to allow for distance restraints\n!!!HACK\n";
        $merge++;
      }
      elsif ($greedymol && ($moleculenumber_atomcount{${$atoms{$atom}}{moleculenumber}} > 3)) {
        $merge ++;
      }
      else {
        $merge = 0;
      }
      if ($merge < 2) {
        if (!defined($moleculeidtogroup{${$atoms{$atom}}{moleculeid}})) { # we haven't seen this moleculeid before
          $groupid = $lastgroupid + 1;
          $lastgroupid = $groupid;
          $saveAtomsToRepresentGroup = 1;
        }
        else { # already had a molecule like this one before - reuse the old groupid, don't add the atoms to OneMoleculeRepresentsEachGroup again
          $groupid = $moleculeidtogroup{${$atoms{$atom}}{moleculeid}};
          $saveAtomsToRepresentGroup = 0;
        }
        #$groupid ++;
      }
      $moleculeidtogroup{${$atoms{$atom}}{moleculeid}} = $groupid;
      $moleculenumbertogroup{${$atoms{$atom}}{moleculenumber}} = $groupid;
    }

    # line below commented out and seems to fix a problem but might have been there for a reason that i've forgotten...
    #my $groupid = $moleculenumbertogroup{${$atoms{$atom}}{moleculeid}};

    ${$atoms{$atom}}{groupid} = $groupid;
    if ($saveAtomsToRepresentGroup) {
      ${$OneMoleculeRepresentsEachGroup{$groupid}}{$atom} = 1;
    }
  }

  #foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
  #  print "$atom ",${$atoms{$atom}}{moleculeid},"\n";
  #}

  my %groupcounts;
  foreach my $group (sort {$a <=> $b} keys (%OneMoleculeRepresentsEachGroup)) {
    $groupcounts{$group} = scalar(keys(%{$OneMoleculeRepresentsEachGroup{$group}}));
    # check there are no discontinuities within molecules
    my $lastatom;
    foreach my $groupatom (sort {$a <=> $b} keys(%{$OneMoleculeRepresentsEachGroup{$group}})) {
      if (defined($lastatom)) {
        if ($groupatom != ($lastatom+1)) {die("Discontinuity within molecule (non-contiguous atoms) - bug\n");}
      }
      $lastatom = $groupatom;
    }
  }

  # try to fish out some group names...?
  my %group_resnames;
  foreach my $atom (sort ({$a <=> $b} keys(%atoms))) {
    ${$group_resnames{${$atoms{$atom}}{groupid}}}{${$atoms{$atom}}{resname}} = 1;
  }

  my %groupnames;
  my %name_dupcheck;
  my $macromolcount = 1;
  print "Naming of groups: anything containing 5 or more different residue names\n";
  print "is 'MACROMOL', otherwise individual residue names are concatenated\n";
  foreach my $group (sort {$a <=> $b} keys(%group_resnames)) {
    if (scalar(keys(%{$group_resnames{$group}})) > 4) {
      my $groupname = sprintf("%s","MACROMOL$macromolcount");
      if (defined($name_dupcheck{$groupname})) {
        my $letter = "A";
        my $groupnametest = $groupname."_$letter";
        while (defined($name_dupcheck{$groupnametest})) {
          $letter ++;
        }
        $groupname .= "_$letter";
      }
      $groupnames{$group} = $groupname;
      $name_dupcheck{$groupname} ++;
      print "Group $group, ",$groupcounts{$group}," atoms: ", $groupnames{$group}, "\n";
      $macromolcount ++;
    }
    else {
      my $groupname = sprintf("%s",join("_",sort({$a cmp $b} keys(%{$group_resnames{$group}}))));
      if (defined($name_dupcheck{$groupname})) {
        my $letter = "A";
        my $groupnametest = $groupname."_$letter";
        while (defined($name_dupcheck{$groupnametest})) {
          $letter ++;
        }
        $groupname .= "_$letter";
      }
      $groupnames{$group} = $groupname;
      $name_dupcheck{$groupname} ++;
      print "Group $group, ",$groupcounts{$group}," atoms: ", $groupnames{$group}, "\n";
    }
  }
  foreach my $name (keys(%name_dupcheck)) {
    if ($name_dupcheck{$name} != 1) {
      die ("Found $name_dupcheck{$name} entries for name $name, which won't work\n");
    }
  }

  my @structureByGroups;
  my @structureByGroupsCount;
  for (my $i = 1; $i <= scalar(keys(%atoms)); ) {
    my $groupname = $groupnames{${$atoms{$i}}{groupid}};
    if ($i == 1) {push (@structureByGroups, $groupname);}
    elsif (! ($groupname eq $structureByGroups[scalar(@structureByGroups) - 1] )) {
      push (@structureByGroups, $groupname); 
    }
    $structureByGroupsCount[scalar(@structureByGroups) - 1] ++;
    #print $groupnames{${$atoms{$i}}{groupid}},"\n";
    $i += $groupcounts{${$atoms{$i}}{groupid}};
  }

  # need to identifty water to write a 'settles' section

  my %groupLooksLikeWater;  
  foreach my $group (keys(%groupcounts)) {
    if ($groupcounts{$group} == 3) {
      my $masssum; my $chargesum;
      my %masses;
      foreach my $atom (keys(%{$OneMoleculeRepresentsEachGroup{$group}})) {
        $masssum += ${$atoms{$atom}}{mass};
        $masses{${$atoms{$atom}}{mass}} ++;
        $chargesum += ${$atoms{$atom}}{charge};
      }
      if ($masssum < 18.1 && $masssum > 17.9 && (abs($chargesum) < 0.001) ) {
        my @massessorted = sort {$a <=> $b} keys(%masses); 
        if ($masses{$massessorted[0]} == 2 && $masses{$massessorted[1]} == 1) {
          if ($massessorted[0] < 1.05 && $massessorted[0] > 0.95 && $massessorted[1] < 16.05 && $massessorted[1] > 15.95) {
            $groupLooksLikeWater{$group} = 1;
          }
        }
      }
    }
  }

  if (scalar(keys(%groupLooksLikeWater)) == 0) {print "******************************\nNo groups identified as water\n******************************\n";}
  if (scalar(keys(%groupLooksLikeWater)) == 1) {print "*************************************************************************************************\n*** Group ",keys(%groupLooksLikeWater)," looks like water (mass is ~18, charge 0, two atoms mass ~1 and one atom mass ~16) ***\n*************************************************************************************************\n";}
  if (scalar(keys(%groupLooksLikeWater)) > 1) {print "The following groups were identified as water: ",join(" ",keys(%groupLooksLikeWater)),"\n"; die ("This might not be a problem but dying anyway as this is unanticipated in this script\n");}

  # print summary
  print "Structure layout by group: \n";
  for (my $i = 0; $i < scalar(@structureByGroups); $i++ ) {
    printf("%-12s x%i\n",$structureByGroups[$i],$structureByGroupsCount[$i]);
  }

  # Assign 'charge' groups - actually just something similar to NAMD's hydrogen groups
  # insofar as they won't necessarily be neutral, just a unit to speed up neighbour searching

  foreach my $group (sort {$a <=> $b} keys(%groupnames)) {
    my $groupcount = 0;
    # this scheme doesn't work for 'triangular' water with 3 bonds...
    # but a special arrangement below identifies water and assigns it all to one charge group
    foreach my $atom (sort {$a <=> $b} keys(%{$OneMoleculeRepresentsEachGroup{$group}})) {
      if (! defined($bondsbyatom2{$atom})) {
        $groupcount ++;
        ${$atoms{$atom}}{cgrp} = $groupcount;
      }
      elsif (scalar(keys(%{$bondsbyatom2{$atom}})) == 1) {
        # terminal atom... only work with this if it's only bound to another terminal atom (2-atom molecule)
        my $boundatom = keys(%{$bondsbyatom2{$atom}});
        if (scalar(keys(%{$bondsbyatom2{$boundatom}})) == 1) {
          $groupcount ++;
          ${$atoms{$atom}}{cgrp} = $groupcount;
          ${$atoms{$boundatom}}{cgrp} = $groupcount;
        }
      }
      else {
        $groupcount ++;
        ${$atoms{$atom}}{cgrp} = $groupcount;
        # start out from every non-hydrogen atom; form group including any bonded atoms which have no further bonds of their own
        foreach my $boundatom (keys(%{$bondsbyatom2{$atom}})) {
          if (scalar(keys(%{$bondsbyatom2{$boundatom}})) == 1) {
            ${$atoms{$boundatom}}{cgrp} = $groupcount;
            # not sure about non-contiguous charge groups... probably safer to weed these out later...?
          }
        }
      }
    }

    # disallow non-contiguous charge groups, just in case
    my $forceContiguousChargeGroups = 1;
    if ($forceContiguousChargeGroups) {
      my $lastchargegroup = 0;
      my $newcgrp;
      foreach my $atom (sort {$a <=> $b} keys(%{$OneMoleculeRepresentsEachGroup{$group}})) {
        if (${$atoms{$atom}}{cgrp} != $lastchargegroup) {
          $newcgrp ++;
        }
        ${$atoms{$atom}}{newcgrp} = $newcgrp;
        $lastchargegroup = ${$atoms{$atom}}{cgrp};
      }
      foreach my $atom (sort {$a <=> $b} keys(%{$OneMoleculeRepresentsEachGroup{$group}})) {
        ${$atoms{$atom}}{cgrp} = ${$atoms{$atom}}{newcgrp};
        delete(${$atoms{$atom}}{newcgrp});
      }
    }
    #print "************Group $groupnames{$group}\n";
    #foreach my $atom (sort {$a <=> $b} keys(%{$OneMoleculeRepresentsEachGroup{$group}})) {
    #  print "$atom ${$atoms{$atom}}{name} ${$atoms{$atom}}{cgrp}\n";
    #}
    
  }



  
  open OUTFILE, (">$_[0]");
  
  print OUTFILE "; Written by WriteTop.pm on ",`date`,"\n";
  print OUTFILE "[ defaults ]\n";
  print OUTFILE "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n";
  print OUTFILE "1               2               yes             0.5     0.8333\n";

  ### print atom types
  my %atomtypes;
  foreach my $atom (sort ({$a <=> $b} keys(%atoms))) {
    my $sig = ${$ljsigeps{${$atoms{$atom}}{atom_type_index}}}{sigma} / 10;
    my $eps = ${$ljsigeps{${$atoms{$atom}}{atom_type_index}}}{epsilon} * $caltoj;
    my $name = ${$atoms{$atom}}{amber_atom_type};
    if (defined($atomtypes{$name})) {
      if ((${$atomtypes{$name}}{sig} != $sig) || (${$atomtypes{$name}}{eps} != $eps)) {
        print "name $name sig $sig (expect ${$atomtypes{$name}}{sig}) eps $eps (expect ${$atomtypes{$name}}{eps}) type ${$atoms{$atom}}{atom_type_index}\n";
        die ("Contradictory LJ information found for amber_atom_tpe $name\n");
      }
    }
    ${$atomtypes{$name}}{sig} = $sig;
    ${$atomtypes{$name}}{eps} = $eps;

    if (defined(${$atoms{$atom}}{atom_type_indexB})) {
      my $sigB = ${$ljsigeps{${$atoms{$atom}}{atom_type_indexB}}}{sigma} / 10;
      my $epsB = ${$ljsigeps{${$atoms{$atom}}{atom_type_indexB}}}{epsilon} * $caltoj;
      my $nameB = ${$atoms{$atom}}{amber_atom_typeB};
      if (defined($atomtypes{$nameB})) {
        if ((${$atomtypes{$nameB}}{sig} != $sigB) || (${$atomtypes{$nameB}}{eps} != $epsB)) {
          print "B state name $nameB sig $sigB (expect ${$atomtypes{$nameB}}{sig}) eps $epsB (expect ${$atomtypes{$nameB}}{eps}) type ${$atoms{$atom}}{atom_type_indexB}\n";
          die ("Contradictory LJ information found for amber_atom_typeB $nameB\n");
        }
      }
      ${$atomtypes{$nameB}}{sig} = $sigB;
      ${$atomtypes{$nameB}}{eps} = $epsB;
    }
  }


  # Re. "Overriding atomtype" warnings... seems gromacs doesn't care about atom type case
  # while GAFF has all lower-case atom types. Mailing list seems to think this isn't an issue
  # cos all but a couple of random types are the same...
  # however it seems like good hygiene to code around this potential pitfall... and get rid of
  # grompp warnings.
  # quick test shows grompp appears alright with 3-letter atom types... let's go down that road.
  
  print OUTFILE "\n[ atomtypes ]\n";
  print OUTFILE ";name  bond_type    mass    charge   ptype          sigma      epsilon\n";
  my %atomtype_avoidduplicate;
  my %dontwrite;
  foreach my $type (sort ({$a cmp $b} keys(%atomtypes))) {
    $dontwrite{$type} = 0;
    if ( !defined ($atomtype_avoidduplicate{lc($type)})) {
      $atomtype_avoidduplicate{lc($type)} = $type;
      $atomtype_avoidduplicate{$type} = $type;
    }
    else {
      my $conflictingtype = $atomtype_avoidduplicate{lc($type)};
      #print "atom type $type already exists as $conflictingtype";
      if ( (${$atomtypes{$type}}{sig} == ${$atomtypes{$conflictingtype}}{sig}) && (${$atomtypes{$type}}{eps} == ${$atomtypes{$conflictingtype}}{eps}) ) {
        print "Found duplicate atom types with same sigma and epsilon, translating $type->$conflictingtype\n";
        $atomtype_avoidduplicate{lc($type)} = $conflictingtype;
        $atomtype_avoidduplicate{$type} = $conflictingtype;
        $dontwrite{$type} = 1;
      }
      else {
        my $newtype = $type."'";
        if (defined($atomtypes{$newtype})) {
          die("trying to define new type $newtype on top of $type, which conflicted with $conflictingtype, but $newtype already exists\n");
        }
        ${$atomtypes{$newtype}}{sig} = ${$atomtypes{$type}}{sig};
        ${$atomtypes{$newtype}}{eps} = ${$atomtypes{$type}}{eps};
        $dontwrite{$type} = 1;
        $dontwrite{$newtype} = 0;
        $atomtype_avoidduplicate{$type} = $newtype;
        print "Created new type $newtype for type $type, which conflicted with previous type $conflictingtype\n";
        #die ("Duplicate atom types ($type $conflictingtype) found with non-identical sigma and epsilon (probably case insensitivity issue). Easily solved (3-letter atom types) but not yet implemented - give up\n");
      }
    }
  }
  foreach my $type (sort ({$a cmp $b} keys(%atomtypes))) {
    if (! $dontwrite{$type}) {
      printf OUTFILE ("%-10s%6s      0.0000  0.0000  A %13.5e%13.5e\n", $type, $type, ${$atomtypes{$type}}{sig}, ${$atomtypes{$type}}{eps});
    }
  }
  
  #foreach my $atomtype (sort {$a <=> $b} keys(%ljsigeps)) {
  #  print " !zz2! type $atomtype sig ${$ljsigeps{$atomtype}}{sigma} eps ${$ljsigeps{$atomtype}}{epsilon}\n";
  #}


  #if (scalar(keys(%groupLooksLikeWater)) == 1) {
  #  print OUTFILE "\n  [ bondtypes ]\n  ; i    j      func       b0          kb\n  OW    HW         1    0.09572   462750.4 ; TIP3P water\n  HW    HW         1    0.15136   462750.4 ; TIP3P water\n  \n  [ angletypes ]\n  ;  i    j    k  func       th0       cth\n  HW  OW  HW           1   104.520    836.800 ; TIP3P water\n  HW  HW  OW           1   127.740      0.000 ; (found in crystallographic water with 3 bonds)\n";
  #}

  ### on to molecules

  my %atomToGroup;
  foreach my $group (keys(%groupnames)) {
    foreach my $atom (keys(%{$OneMoleculeRepresentsEachGroup{$group}})) {
      $atomToGroup{$atom} = $group;
    }
  }

  my %relevantbonds;
  foreach my $bond (keys(%bonds)) {
    my $a1 = ${$bonds{$bond}}{a1};
    if (defined($atomToGroup{$a1})) {
      push (@{$relevantbonds{$atomToGroup{$a1}}}, $bond);
    }
  }
  my %relevantangles;
  foreach my $angle (keys(%angles)) {
    my $a1 = ${$angles{$angle}}{a1};
    if (defined($atomToGroup{$a1})) {
      push (@{$relevantangles{$atomToGroup{$a1}}}, $angle);
    }
  }
  my %relevantdihedrals;
  foreach my $dihedral (keys(%dihedrals)) {
    my $a1 = ${$dihedrals{$dihedral}}{a1};
    if (defined($atomToGroup{$a1})) {
      push (@{$relevantdihedrals{$atomToGroup{$a1}}}, $dihedral);
    }
  }
    
  
  foreach my $group (sort {$a <=> $b} keys(%groupnames)) {
    my $iswater = 0;
    if ( defined($groupLooksLikeWater{$group})) {$iswater = 1;}

    print OUTFILE "\n\n[ moleculetype ]\n; Name            nrexcl\n";
    if ($iswater) {printf OUTFILE ("%-16s  1\n",$groupnames{$group});}
    elsif (scalar(keys(%{$OneMoleculeRepresentsEachGroup{$group}})) == 1) {printf OUTFILE ("%-16s  0\n",$groupnames{$group});}
    else {printf OUTFILE ("%-16s  3\n",$groupnames{$group});}

    # atoms
    print OUTFILE "\n[ atoms ]\n;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB\n";
    my $firstatom = -1; my $firstres = -1;
    my $posres;
    foreach my $atom (sort {$a <=> $b} keys(%{$OneMoleculeRepresentsEachGroup{$group}})) {
      if ($firstatom == -1) {$firstatom = $atom;}
      if ($firstres == -1) {$firstres = ${$atoms{$atom}}{resid};}
      my $atomrenum = $atom - $firstatom + 1;
      my $resnumrenum = ${$atoms{$atom}}{resid} - $firstres + 1;
      my $safetype = $atomtype_avoidduplicate{${$atoms{$atom}}{amber_atom_type}};
      my $name = ${$atoms{$atom}}{name};
      my $resname = ${$atoms{$atom}}{resname};
      #if (defined(${$atoms{$atom}}{nameB}) && ($name eq "DMY")) {
      #  $name = ${$atoms{$atom}}{nameB};
      #}
      my $safetypeB;
      if (defined(${$atoms{$atom}}{amber_atom_typeB})) {
        $safetypeB = $atomtype_avoidduplicate{${$atoms{$atom}}{amber_atom_typeB}};
        #print "defined safetypeB: type ${$atoms{$atom}}{amber_atom_typeB} safetype $safetypeB\n";
      }
      if ($iswater) {
        printf OUTFILE "%6d %10s %6d %6s %6s %6d %10.5f %10f\n", $atomrenum, $safetype, $resnumrenum, ${$atoms{$atom}}{resname}, $name, 1, ${$atoms{$atom}}{charge}, ${$atoms{$atom}}{mass};
      }
      else {
        my $pertcount = 0;
        if (defined(${$atoms{$atom}}{chargeB})) {$pertcount ++;}
        if (defined(${$atoms{$atom}}{amber_atom_typeB})) {$pertcount ++;}
        if (defined(${$atoms{$atom}}{massB})) {$pertcount ++;}
        if ( ($pertcount != 0) && ($pertcount != 3) ) {
          die ("Found one or more B state parameters (chargeB, amber_atom_typeB, massB) for atom $atom, but not all 3\n");
        }
        if ($pertcount == 3) {
          my $parametersDiffer = 0;
          if (defined(${$atoms{$atom}}{nameB}) && (! (${$atoms{$atom}}{name} eq ${$atoms{$atom}}{nameB}))) {
            $name .= "-x-".${$atoms{$atom}}{nameB};
          }
          if (defined(${$atoms{$atom}}{resnameB}) && (! (${$atoms{$atom}}{resname} eq ${$atoms{$atom}}{resnameB}))) {
            $resname .= "-x-".${$atoms{$atom}}{resnameB};
          }
          printf OUTFILE "%6d %10s %6d %6s %6s %6d %10.5f %10f", $atomrenum, $safetype, $resnumrenum, $resname, $name, ${$atoms{$atom}}{cgrp}, ${$atoms{$atom}}{charge}, ${$atoms{$atom}}{mass};
          if (${$atoms{$atom}}{chargeB} != ${$atoms{$atom}}{charge}) {$parametersDiffer++;}
          if (${$atoms{$atom}}{massB} != ${$atoms{$atom}}{mass}) {$parametersDiffer++;}
          if (! (${$atoms{$atom}}{amber_atom_typeB} eq ${$atoms{$atom}}{amber_atom_type})) {$parametersDiffer++;}
          if ($parametersDiffer) {
            #print "atom $atom type ${$atoms{$atom}}{amber_atom_type} safetype $safetype typeB ${$atoms{$atom}}{amber_atom_typeB} safetypeB $safetypeB\n";
            printf OUTFILE "%10s %10.5f %10f", $safetypeB, ${$atoms{$atom}}{chargeB}, ${$atoms{$atom}}{massB};
            #if (defined(${$atoms{$atom}}{nameB})) {
            #  printf OUTFILE " ;   nameB %s",${$atoms{$atom}}{nameB};
            #}
          }
        }
        else {
          printf OUTFILE "%6d %10s %6d %6s %6s %6d %10.5f %10f", $atomrenum, $safetype, $resnumrenum, ${$atoms{$atom}}{resname}, $name, ${$atoms{$atom}}{cgrp}, ${$atoms{$atom}}{charge}, ${$atoms{$atom}}{mass};
        }
        print OUTFILE "\n";
      }
      if (defined(${$atoms{$atom}}{restraint_kB})) {
        if ((!defined(${$atoms{$atom}}{restraint_k})) || (abs(${$atoms{$atom}}{restraint_k} - ${$atoms{$atom}}{restraint_kB})>0.001)) {
          print "atom $atom restraint_k ${$atoms{$atom}}{restraint_k} restraint_kB ${$atoms{$atom}}{restraint_kB}\n";
          die("Different restraint_k entries for A and B states not supported yet\n");
        }
      }
      if (defined(${$atoms{$atom}}{restraint_k})) {
        if ($newformatPosres) {
          if (${$atoms{$atom}}{restraint_fbhw}) {
            if (${$atoms{$atom}}{restraint_fbhw} < 0) {
              $posres .= sprintf("%6d %6d %6d %8.3f %8.3f %8.3f %8.3f\n",$atomrenum, 1, 1, ${$atoms{$atom}}{restraint_fbhw} * (-0.1), ${$atoms{$atom}}{restraint_k} * 200 * $caltoj,${$atoms{$atom}}{restraint_k} * 200 * $caltoj, ${$atoms{$atom}}{restraint_k} * 200 * $caltoj);
            }
            else {
              $posres .= sprintf("%6d %6d %6d %8.3f %8.3f %8.3f %8.3f\n",$atomrenum, 1, 0, ${$atoms{$atom}}{restraint_fbhw} * (0.1), ${$atoms{$atom}}{restraint_k} * 200 * $caltoj,${$atoms{$atom}}{restraint_k} * 200 * $caltoj, ${$atoms{$atom}}{restraint_k} * 200 * $caltoj);
            }
          }
          else {
            $posres .= sprintf("%6d %6d %6d %8.3f %8.3f %8.3f %8.3f\n",$atomrenum, 1, 0, 0, ${$atoms{$atom}}{restraint_k} * 200 * $caltoj,${$atoms{$atom}}{restraint_k} * 200 * $caltoj, ${$atoms{$atom}}{restraint_k} * 200 * $caltoj);
          }
        }
        else {
          if (defined(${$atoms{$atom}}{restraint_fbhw})) {
            $posres .= sprintf("%6d %6d %8.3f %8.3f %8.3f %i 0 0\n",$atomrenum, 1, ${$atoms{$atom}}{restraint_k} * 200 * $caltoj,${$atoms{$atom}}{restraint_k} * 200 * $caltoj, ${$atoms{$atom}}{restraint_k} * 200 * $caltoj,${$atoms{$atom}}{restraint_fbhw});
          }
          else {
            $posres .= sprintf("%6d %6d %8.3f %8.3f %8.3f\n",$atomrenum, 1, ${$atoms{$atom}}{restraint_k} * 200 * $caltoj,${$atoms{$atom}}{restraint_k} * 200 * $caltoj, ${$atoms{$atom}}{restraint_k} * 200 * $caltoj);
          }
        }
      }
    }

    # bonds
    if (! $iswater) {
      print OUTFILE "\n[ bonds ]\n;  ai    aj funct  r  k\n";
      if (defined($relevantbonds{$group})) {
        foreach my $bond (sort {${$bonds{$a}}{a1} <=> ${$bonds{$b}}{a1} || ${$bonds{$a}}{a2} <=> ${$bonds{$b}}{a2}} @{$relevantbonds{$group}}) {
          # amb2gmx multiplies bond force constant by 200, no idea why
          if (defined(${$bonds{$bond}}{type})) {
            printf OUTFILE "%5d%6d%6d%12.4e%12.4e\n",${$bonds{$bond}}{a1} - $firstatom + 1, ${$bonds{$bond}}{a2} - $firstatom + 1, 1, ${$bond_types{${$bonds{$bond}}{type}}}{equil_value} / 10, ${$bond_types{${$bonds{$bond}}{type}}}{force_constant} * $caltoj * 200;
          }
          else {
            if (! $perturbedBonds) {
              die("Bond $bond has no type entry - this is ok for alchemical topologies but no force_constantB entries were previously found\n");
            }
            if ( (!defined(${$bonds{$bond}}{force_constantA})) || (!defined(${$bonds{$bond}}{force_constantB})) || (!defined(${$bonds{$bond}}{equil_valueA})) || (!defined(${$bonds{$bond}}{equil_valueB})) ) {
              die("No bond type and missing or incomplete A and B force_constant / equil_value entries\n");
            }
            printf OUTFILE "%5d%6d%6d%12.4e%12.4e",${$bonds{$bond}}{a1} - $firstatom + 1, ${$bonds{$bond}}{a2} - $firstatom + 1, 1, ${$bonds{$bond}}{equil_valueA} / 10, ${$bonds{$bond}}{force_constantA} * $caltoj * 200, ${$bonds{$bond}}{equil_valueB} / 10, ${$bonds{$bond}}{force_constantB} * $caltoj * 200;
            if ( (${$bonds{$bond}}{equil_valueB} != ${$bonds{$bond}}{equil_valueA}) || (${$bonds{$bond}}{force_constantB} != ${$bonds{$bond}}{force_constantA}) ) {
              printf OUTFILE "%12.4e%12.4e",${$bonds{$bond}}{equil_valueB} / 10, ${$bonds{$bond}}{force_constantB} * $caltoj * 200;
            }
            print OUTFILE "\n";
          }
        }
      }
    }

    print OUTFILE "\n[ pairs ]\n;  ai    aj funct\n";
    if (scalar(keys(%pairs))) {
      # 1-4 pairs provided by calling script
      my $countall = 0;
      my $countnonstandard = 0;
      #print "[ pairs ] section: writing pairs as found in \%pairs hash passed to WriteGTop\n";
      my %sortkeys;
      my %sortkeys2;
      foreach my $pair (keys(%pairs)) {
        my @temp = split(":",$pair);
        $sortkeys{$pair} = $temp[0];
        $sortkeys2{$pair} = $temp[1];
      }
      foreach my $pair (sort {$sortkeys{$a} <=> $sortkeys{$b} || $sortkeys2{$a} <=> $sortkeys2{$b}} keys(%pairs)) {
        my @temp = split(":",$pair);
        if (defined($atomToGroup{$temp[0]}) || defined($atomToGroup{$temp[1]})) {
          if (! (defined($atomToGroup{$temp[0]}) && defined($atomToGroup{$temp[1]}))) {
            die("[ pairs ] section contains an entry ($pair) that crosses molecule types - this is a bug of some sort\n");
          }
          if ($atomToGroup{$temp[0]} == $group) {
            if ($atomToGroup{$temp[1]} != $group) {
              die("[ pairs ] section contains an entry ($pair) that crosses molecule types - this is a bug of some sort\n");
            }
            foreach my $ftype (keys(%{$pairs{$pair}})) {
              $countall ++;
              if ($ftype > 2) {$countnonstandard ++;}
              printf OUTFILE ("%6d %6d %6d\n",$temp[0] - $firstatom + 1, $temp[1] - $firstatom + 1, $ftype);
            }
          }
        }
      }
      if ($countall) {
        print "Group $groupnames{$group}, [ pairs ] section: writing pairs as found in \%pairs hash passed to WriteGTop\n";
        if ($countnonstandard) {
          print "$countnonstandard [ pairs ] of types 3, 4, 5 or 6 written. These are nonstandard types requiring FB's hacked gromacs\n";
        }
      }
    }
    else {
      # 1-4 pairs... getting them from connectivity might seem more general but might screw up dual-topology FEC later on...
      # so copy amb2gmx.pl again and just take 1s and 4s from dihedrals.
#      print "Constructing [ pairs ] section with 1-4s based on dihedrals\n";
#      my %pairs_entry;
#      if (defined($relevantdihedrals{$group})) {
#        foreach my $dihedral (sort {$a <=> $b} @{$relevantdihedrals{$group}}) {
#          my $a1 = ${$dihedrals{$dihedral}}{a1};
#          my $a2 = abs(${$dihedrals{$dihedral}}{a2});
#          my $a3 = abs(${$dihedrals{$dihedral}}{a3});
#          my $a4 = abs(${$dihedrals{$dihedral}}{a4});
#          # make sure it's a proper dihedral
#          if (defined(${$bondsbyatom2{$a1}}{$a2}) && defined(${$bondsbyatom2{$a2}}{$a3}) && defined(${$bondsbyatom2{$a3}}{$a4})) {
#            # doh another pitfall... check a1 and a4 have no bonded atoms in common, otherwise you can get spurious 1-3 pairs eg. in 5-membered ring
#            my %bondedtoa1anda4;
#            my $problem = 0; 
#            foreach my $bondedatom (keys(%{$bondsbyatom2{$a1}}) ) {$bondedtoa1anda4{$bondedatom} ++;}
#            foreach my $bondedatom (keys(%{$bondsbyatom2{$a4}}) ) {$bondedtoa1anda4{$bondedatom} ++;}
#            foreach my $val (values(%bondedtoa1anda4)) {
#              if ($val > 1) {$problem = 1;}
#            }
#            if (! $problem) {
#              my $a_lo; my $a_hi;
#              if ($a1 < $a4) {$a_lo = $a1; $a_hi = $a4;}
#              elsif ($a1 > $a4) {$a_lo = $a4; $a_hi = $a1;}
#              else {die ("huh? 1-4 of atom with itself\n");}
#              my $entry = sprintf("%6d %6d %6d\n",$a_lo - $firstatom + 1, $a_hi - $firstatom + 1, 1);
#              $pairs_entry{$entry} = $a_lo;
#            }
#          }
#        }
#      }
      # scratch that: forget dihedrals, traverse bonds to determine 1-4s
      my %pairs_entry;
      my %excl;
      my %sortkeys;
      my %sortkeys2;
      foreach my $a1 (sort {$a <=> $b} keys(%{$OneMoleculeRepresentsEachGroup{$group}})) {
        $excl{$a1} = 1; # avoid 1-4 pairs of atom with itself in triangular groups
        foreach my $a2 (keys(%{$bondsbyatom2{$a1}})) {
          my $excl12;
          if ($a1 < $a2) {$excl12 = sprintf("%6d %6d %6d\n",$a1 - $firstatom + 1, $a2 - $firstatom + 1, 1);}
          else {$excl12 = sprintf("%6d %6d %6d\n",$a2 - $firstatom + 1, $a1 - $firstatom + 1, 1);}
          $excl{$excl12} = 1;
          foreach my $a3 (keys(%{$bondsbyatom2{$a2}})) {
            my $excl13;
            if ($a1 < $a3) {$excl13 = sprintf("%6d %6d %6d\n",$a1 - $firstatom + 1, $a3 - $firstatom + 1, 1);}
            else {$excl13 = sprintf("%6d %6d %6d\n",$a3 - $firstatom + 1, $a1 - $firstatom + 1, 1);}
            $excl{$excl13} = 1;
          }
        }
      }
      foreach my $a1 (sort {$a <=> $b} keys(%{$OneMoleculeRepresentsEachGroup{$group}})) {
        foreach my $a2 (keys(%{$bondsbyatom2{$a1}})) {
          foreach my $a3 (keys(%{$bondsbyatom2{$a2}})) {
            foreach my $a4 (keys(%{$bondsbyatom2{$a3}})) {
              my $a_lo; my $a_hi;
              if ($a1 < $a4) {$a_lo = $a1; $a_hi = $a4;}
              elsif ($a1 > $a4) {$a_lo = $a4; $a_hi = $a1;}
              else {
                my $pair14 = sprintf("%6d %6d %6d\n",$a1 - $firstatom + 1, $a4 - $firstatom + 1, 1);
                if (!defined($excl{$pair14})) {
                  die "huh? a1 $a1 == a4 $a4, but no exclusion defined\n";
                }
                next;
              }
              my $pair14 = sprintf("%6d %6d %6d\n",$a_lo - $firstatom + 1, $a_hi - $firstatom + 1, 1);
              if (!defined($excl{$pair14})) {$pairs_entry{$pair14} = $a_lo;}
            }
          }
        }
      }
      foreach my $pair (keys(%pairs_entry)) {
        my @temp = split(" ",$pair);
        $sortkeys{$pair} = $temp[0];
        $sortkeys2{$pair} = $temp[1];
      }
      
      foreach my $entry (sort {$sortkeys{$a} <=> $sortkeys{$b} || $sortkeys2{$a} <=> $sortkeys2{$b}} keys(%pairs_entry)) {
        print OUTFILE $entry;
      }
    }

    # angles
    print OUTFILE "\n[ angles ]\n;  ai    aj    ak funct  theta   cth\n";
    if (defined($relevantangles{$group})) {
      #foreach my $angle (sort {$a <=> $b} @{$relevantangles{$group}}) {
      foreach my $angle (sort {${$angles{$a}}{a1} <=> ${$angles{$b}}{a1} || ${$angles{$a}}{a2} <=> ${$angles{$b}}{a2} || ${$angles{$a}}{a3} <=> ${$angles{$b}}{a3}} @{$relevantangles{$group}}) {
        my $a1 = ${$angles{$angle}}{a1};
        my $a2 = ${$angles{$angle}}{a2};
        my $a3 = ${$angles{$angle}}{a3};
        if (defined(${$angles{$angle}}{type})) {
          printf OUTFILE ("%5d%6d%6d%6d%12.4e%12.4e\n",$a1 - $firstatom + 1, $a2 - $firstatom + 1, $a3 - $firstatom + 1, 1, ${$angle_types{${$angles{$angle}}{type}}}{equil_value} * $degtorad, ${$angle_types{${$angles{$angle}}{type}}}{force_constant} * $caltoj * 2);
        }
        else {
          if (! $perturbedAngles) {
            die("Angle $angle has no type entry - this is ok for alchemical topologies but no force_constantB entries were previously found\n");
          }
          if ( (!defined(${$angles{$angle}}{force_constantA})) || (!defined(${$angles{$angle}}{force_constantB})) || (!defined(${$angles{$angle}}{equil_valueA})) || (!defined(${$angles{$angle}}{equil_valueB})) ) {
            die("No angle type and missing or incomplete A and B force_constant / equil_value entries\n");
          }
          printf OUTFILE ("%5d%6d%6d%6d%12.4e%12.4e",$a1 - $firstatom + 1, $a2 - $firstatom + 1, $a3 - $firstatom + 1, 1, ${$angles{$angle}}{equil_valueA} * $degtorad, ${$angles{$angle}}{force_constantA} * $caltoj * 2);
          if ( (${$angles{$angle}}{equil_valueB} != ${$angles{$angle}}{equil_valueA}) || (${$angles{$angle}}{force_constantB} != ${$angles{$angle}}{force_constantA}) ) {
            printf OUTFILE ("%12.4e%12.4e",${$angles{$angle}}{equil_valueB} * $degtorad, ${$angles{$angle}}{force_constantB} * $caltoj * 2);
          }
          print OUTFILE "\n";
        }
      }
    }
    
    # dihedrals
    print OUTFILE "\n[ dihedrals ]\n;  ai    aj    ak    al   funct  theta    k     multiplicity\n";
    if (defined($relevantdihedrals{$group})) {
      #foreach my $dihedral (sort {$a <=> $b} @{$relevantdihedrals{$group}}) {
      foreach my $dihedral (sort {abs(${$dihedrals{$a}}{a1}) <=> abs(${$dihedrals{$b}}{a1}) || abs(${$dihedrals{$a}}{a2}) <=> abs(${$dihedrals{$b}}{a2}) || abs(${$dihedrals{$a}}{a3}) <=> abs(${$dihedrals{$b}}{a3}) || abs(${$dihedrals{$a}}{a4}) <=> abs(${$dihedrals{$b}}{a4})} @{$relevantdihedrals{$group}}) {
        my $a1 = abs(${$dihedrals{$dihedral}}{a1});
        my $a2 = abs(${$dihedrals{$dihedral}}{a2});
        my $a3 = abs(${$dihedrals{$dihedral}}{a3});
        my $a4 = abs(${$dihedrals{$dihedral}}{a4});
        if (($a1 - $firstatom + 1) < 0) {print "below zero ",$a1 - $firstatom + 1," a $a1 $a2 $a3 $a4 firstatom $firstatom\n"; die();} 
        if (($a2 - $firstatom + 1) < 0) {print "below zero ",$a2 - $firstatom + 1," a $a1 $a2 $a3 $a4 firstatom $firstatom\n"; die();} 
        if (($a3 - $firstatom + 1) < 0) {print "below zero ",$a3 - $firstatom + 1," a $a1 $a2 $a3 $a4 firstatom $firstatom\n"; die();} 
        if (($a4 - $firstatom + 1) < 0) {print "below zero ",$a4 - $firstatom + 1," a $a1 $a2 $a3 $a4 firstatom $firstatom\n"; die();} 
        if (defined(${$dihedrals{$dihedral}}{type})) {
          printf OUTFILE ("%5d%6d%6d%6d%6d%12.4e%12.4e%6d\n",$a1 - $firstatom + 1, $a2 - $firstatom + 1, $a3 - $firstatom + 1, $a4 - $firstatom + 1, 1, ${$dihedral_types{${$dihedrals{$dihedral}}{type}}}{phase} * $degtorad, ${$dihedral_types{${$dihedrals{$dihedral}}{type}}}{force_constant} * $caltoj, ${$dihedral_types{${$dihedrals{$dihedral}}{type}}}{periodicity});
        }
        else {
          if (! $perturbedDihedrals) {
            die("Dihedral $dihedral has no type entry - this is ok for alchemical topologies but no force_constantB entries were previously found\n");
          }
          if ( (!defined(${$dihedrals{$dihedral}}{force_constantA})) || (!defined(${$dihedrals{$dihedral}}{force_constantB})) 
          || (!defined(${$dihedrals{$dihedral}}{periodicityA})) || (!defined(${$dihedrals{$dihedral}}{periodicityA})) 
          || (!defined(${$dihedrals{$dihedral}}{phaseA})) || (!defined(${$dihedrals{$dihedral}}{phaseB})) ) {
            die("No dihedral type and missing or incomplete A and B force_constant / periodicity / phase entries\n");
          }
          printf OUTFILE ("%5d%6d%6d%6d%6d%12.4e%12.4e%6d",$a1 - $firstatom + 1, $a2 - $firstatom + 1, $a3 - $firstatom + 1, $a4 - $firstatom + 1, 1, ${$dihedrals{$dihedral}}{phaseA} * $degtorad, ${$dihedrals{$dihedral}}{force_constantA} * $caltoj, ${$dihedrals{$dihedral}}{periodicityA});
          if ( (${$dihedrals{$dihedral}}{periodicityB} != ${$dihedrals{$dihedral}}{periodicityA}) || (${$dihedrals{$dihedral}}{force_constantB} != ${$dihedrals{$dihedral}}{force_constantA}) || (${$dihedrals{$dihedral}}{phaseB} != ${$dihedrals{$dihedral}}{phaseA}) ) {
            printf OUTFILE ("%12.4e%12.4e%6d",${$dihedrals{$dihedral}}{phaseB} * $degtorad, ${$dihedrals{$dihedral}}{force_constantB} * $caltoj, ${$dihedrals{$dihedral}}{periodicityB});
          }
          print OUTFILE "\n";
        }
      }
    }

    if ($posres) {
      print OUTFILE "\n[ position_restraints ]\n";
      print OUTFILE $posres;
    }

    if ($iswater) {
      print OUTFILE "\n[ settles ]\n; OW    funct   doh     dhh\n1       1       0.09572 0.15139\n";
      # turns out mdrun needs these explicit exclusions for fast water loops
      print OUTFILE "\n[ exclusions ]\n1   2   3\n2   1   3\n3   1   2\n";
    }
    
  }

  print OUTFILE "\n\n[ system ]\nWriteGtop FB\n";
  print OUTFILE "\n\n[ molecules ]\n; Compound        nmols\n";
  for (my $i = 0; $i < scalar(@structureByGroups); $i++ ) {
    printf OUTFILE ("%-18s %i\n",$structureByGroups[$i],$structureByGroupsCount[$i]);
  }
  
  
  close(OUTFILE);
}

sub AddAtomsToMol { 
  #AddAtomsToMol(\%bondsbyatom, \@molecules, $startatom, \%newmol);
  if (defined( ${$_[0]}{$_[2]})) {
    ${$_[3]}{$_[2]} = 1;
    if(defined(${${$_[0]}{$_[2]}}{-1})) {delete(${${$_[0]}{$_[2]}}{-1});}
    my @bonds; 
    foreach my $bond (keys(%{${$_[0]}{$_[2]}})) {
      if (defined(${$_[0]}{$bond})) {push (@bonds, $bond);}
    }
    delete ${$_[0]}{$_[2]};
    foreach my $bond (@bonds) {
      AddAtomsToMol($_[0], $_[1], $bond, $_[3]);
    }
  }
}

1;
