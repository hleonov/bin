package ToputilsShared;
use strict;
use Exporter;  
our @ISA = qw(Exporter);
our @EXPORT = qw(Cleanup addGap trim);

sub Cleanup {
  my @collected = @{$_[0]};
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


  # atom renumbering - make sure everything is contiguous

  if (1) {
    my %oldtonew;
    my $offset = 0;
    my $last = 0;
    my $gapcount = 0;
    foreach my $atom (sort { $a <=> $b } keys(%atoms)) {
      if ($atom - $last != 1) {
        $offset += ($atom - $last - 1);
        $gapcount ++;
      }
      $oldtonew{$atom} = $atom - $offset;
      $last = $atom;
    }
    print "$gapcount gaps, renumbering... ";
    if ($offset != 0) {
      print "atoms... ";
      foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
        my @newexclusions;
        foreach my $excludedatom (keys(%{${$atoms{$atom}}{excluded}})) {
          push @newexclusions, $oldtonew{$excludedatom};
          delete(${${$atoms{$atom}}{excluded}}{$excludedatom});
        }
        foreach my $newexcludedatom (@newexclusions) {
          ${${$atoms{$atom}}{excluded}}{$newexcludedatom} = 1;
        }
        if ($oldtonew{$atom} != $atom) {
            $atoms{$oldtonew{$atom}} = $atoms{$atom};
            delete($atoms{$atom});
        }
      }
  
      print "bonds... ";
      #per-bond stuff
      foreach my $bond (keys(%bonds)) {
        ${$bonds{$bond}}{a1} = $oldtonew{${$bonds{$bond}}{a1}};
        ${$bonds{$bond}}{a2} = $oldtonew{${$bonds{$bond}}{a2}};
      }

      if (scalar(keys(%pairs))) {
        my %newpairs;
        print "pairs... ";
        foreach my $pair (keys(%pairs)) {
          my @pairarray  = split(":",$pair);
          $pairarray[0] = $oldtonew{$pairarray[0]};
          $pairarray[1] = $oldtonew{$pairarray[1]};
          my $newpair = join(":",sort {$a <=> $b} @pairarray);
          $newpairs{$newpair} = $pairs{$pair};
        }
        undef %pairs;
        %pairs = %newpairs;
      }
  
      print "angles... ";
      #per-angle stuff
      foreach my $angle (keys(%angles)) {
        ${$angles{$angle}}{a1} = $oldtonew{${$angles{$angle}}{a1}};
        ${$angles{$angle}}{a2} = $oldtonew{${$angles{$angle}}{a2}};
        ${$angles{$angle}}{a3} = $oldtonew{${$angles{$angle}}{a3}};
      }
  
      print "dihedrals... ";
      #per-dihedral stuff
      foreach my $dihedral (keys(%dihedrals)) {
        my $a1neg = 0; my $a2neg = 0; my $a3neg = 0; my $a4neg = 0; 
        if (${$dihedrals{$dihedral}}{a1} < 0) {$a1neg = 1; ${$dihedrals{$dihedral}}{a1} *= -1;}
        if (${$dihedrals{$dihedral}}{a2} < 0) {$a2neg = 1; ${$dihedrals{$dihedral}}{a2} *= -1;}
        if (${$dihedrals{$dihedral}}{a3} < 0) {$a3neg = 1; ${$dihedrals{$dihedral}}{a3} *= -1;}
        if (${$dihedrals{$dihedral}}{a4} < 0) {$a4neg = 1; ${$dihedrals{$dihedral}}{a4} *= -1;}

        ${$dihedrals{$dihedral}}{a1} = $oldtonew{${$dihedrals{$dihedral}}{a1}};
        ${$dihedrals{$dihedral}}{a2} = $oldtonew{${$dihedrals{$dihedral}}{a2}};
        ${$dihedrals{$dihedral}}{a3} = $oldtonew{${$dihedrals{$dihedral}}{a3}};
        ${$dihedrals{$dihedral}}{a4} = $oldtonew{${$dihedrals{$dihedral}}{a4}};
  
        if ($a1neg) {${$dihedrals{$dihedral}}{a1} *= -1;}
        if ($a2neg) {${$dihedrals{$dihedral}}{a2} *= -1;}
        if ($a3neg) {${$dihedrals{$dihedral}}{a3} *= -1;}
        if ($a4neg) {${$dihedrals{$dihedral}}{a4} *= -1;}
      }
      print "\n";
    }
  }
  else {
    my $reachedtheend = 0;
    while (!$reachedtheend) {
      $reachedtheend = 1;
      my $last = 0;
      foreach my $atom (sort { $a <=> $b } keys(%atoms)) {
        if ($atom - $last != 1) { # identified a gap
          $reachedtheend = 0;
          my $lastbeforegap = $last;
          my $firstaftergap = $atom;
          my $gapsize = $atom - $last - 1;
          print "Toputils Cleanup: found gap, between atoms ",$last," and ",$atom," size ",$atom - $last - 1," - renumbering\n";
          # cycle through bonds, angles, dihedrals, reducing atom numbers as necessary
          # also go through each individual atom and renumber excluded atoms as necessary
          # example: atom 5 is followed by atom 10
          # gap size = 4
          # need to subtract 4 from every atom ID > 5
  
          # per-atom stuff (ID and exclusions)
          foreach my $atom (sort { $a <=> $b } keys(%atoms)) {
            my @newexclusions;
            foreach my $excludedatom (keys(%{${$atoms{$atom}}{excluded}})) {
              if ($excludedatom > $lastbeforegap) {
                push @newexclusions, $excludedatom - $gapsize;
                delete(${${$atoms{$atom}}{excluded}}{$excludedatom});
              }
            }
            foreach my $newexcludedatom (@newexclusions) {
              ${${$atoms{$atom}}{excluded}}{$newexcludedatom} = 1;
            }
            if ($atom > $lastbeforegap) {
                $atoms{$atom - $gapsize} = $atoms{$atom};
                #print "move to $atom to ",$atom-$gapsize," ";
                delete($atoms{$atom});
                #print "and delete atom $atom. Atoms left: ",scalar(keys(%atoms)),"\n"; 
            }
          }
  
          #per-bond stuff
          foreach my $bond (keys(%bonds)) {
            if (${$bonds{$bond}}{a1} > $lastbeforegap) {${$bonds{$bond}}{a1} = ${$bonds{$bond}}{a1} - $gapsize;}
            if (${$bonds{$bond}}{a2} > $lastbeforegap) {${$bonds{$bond}}{a2} = ${$bonds{$bond}}{a2} - $gapsize;}
          }
  
          #per-angle stuff
          foreach my $angle (keys(%angles)) {
            if (${$angles{$angle}}{a1} > $lastbeforegap) {${$angles{$angle}}{a1} = ${$angles{$angle}}{a1} - $gapsize;}
            if (${$angles{$angle}}{a2} > $lastbeforegap) {${$angles{$angle}}{a2} = ${$angles{$angle}}{a2} - $gapsize;}
            if (${$angles{$angle}}{a3} > $lastbeforegap) {${$angles{$angle}}{a3} = ${$angles{$angle}}{a3} - $gapsize;}
          }
  
          #per-dihedral stuff
          foreach my $dihedral (keys(%dihedrals)) {
            my $a1neg = 0; my $a2neg = 0; my $a3neg = 0; my $a4neg = 0; 
            if (${$dihedrals{$dihedral}}{a1} < 0) {$a1neg = 1; ${$dihedrals{$dihedral}}{a1} *= -1;}
            if (${$dihedrals{$dihedral}}{a2} < 0) {$a2neg = 1; ${$dihedrals{$dihedral}}{a2} *= -1;}
            if (${$dihedrals{$dihedral}}{a3} < 0) {$a3neg = 1; ${$dihedrals{$dihedral}}{a3} *= -1;}
            if (${$dihedrals{$dihedral}}{a4} < 0) {$a4neg = 1; ${$dihedrals{$dihedral}}{a4} *= -1;}
  
            if (${$dihedrals{$dihedral}}{a1} > $lastbeforegap) {${$dihedrals{$dihedral}}{a1} = ${$dihedrals{$dihedral}}{a1} - $gapsize;}
            if (${$dihedrals{$dihedral}}{a2} > $lastbeforegap) {${$dihedrals{$dihedral}}{a2} = ${$dihedrals{$dihedral}}{a2} - $gapsize;}
            if (${$dihedrals{$dihedral}}{a3} > $lastbeforegap) {${$dihedrals{$dihedral}}{a3} = ${$dihedrals{$dihedral}}{a3} - $gapsize;}
            if (${$dihedrals{$dihedral}}{a4} > $lastbeforegap) {${$dihedrals{$dihedral}}{a4} = ${$dihedrals{$dihedral}}{a4} - $gapsize;}
  
            if ($a1neg) {${$dihedrals{$dihedral}}{a1} *= -1;}
            if ($a2neg) {${$dihedrals{$dihedral}}{a2} *= -1;}
            if ($a3neg) {${$dihedrals{$dihedral}}{a3} *= -1;}
            if ($a4neg) {${$dihedrals{$dihedral}}{a4} *= -1;}
          }
          last;
        }
        $last = $atom;
      }
    }
  }

  ## sort out reversed exclusions (excluded atom id should always be higher than origin atom)
  #foreach my $atom (sort ({$a <=> $b} keys(%atoms))) {
  #  #foreach my $excludedatom (keys(%{${$atoms{$atom}}{excluded}})) {
  #  #  #print "$atom $excludedatom\n";
  #  #  if (($excludedatom) && ($excludedatom < $atom)) {
  #  #    #print "$atom > $excludedatom\n";
  #  #    ${${$atoms{$excludedatom}}{excluded}}{$atom} = 1;
  #  #    delete ${${$atoms{$atom}}{excluded}}{excludedatom};
  #  #  }
  #  #}
  #  if ($atom ) {
  #    foreach my $excludedatom (keys(%{${$atoms{$atom}}{excluded}})) {
  #      print "*2* $atom $excludedatom\n";
  #      if (($excludedatom) && ($excludedatom < $atom)) {
  #        print "*2* $atom > $excludedatom\n";
  #      }
  #    }
  #  }
  #}


  if (1) {

  # merge duplicate and delete unused bond, angle and dihedral types
  my %fingerprints_bond;
  my %bondtypes_oldtonew;
  foreach my $bondtype (sort {$a <=> $b} keys(%bond_types)) {
    my $fingerprint = ${$bond_types{$bondtype}}{force_constant}.":".${$bond_types{$bondtype}}{equil_value};
    if (defined($fingerprints_bond{$fingerprint})) {
      $bondtypes_oldtonew{$bondtype} = $fingerprints_bond{$fingerprint};
      delete($bond_types{$bondtype});
    }
    else {
      $fingerprints_bond{$fingerprint} = $bondtype;
      $bondtypes_oldtonew{$bondtype} = $bondtype;
    }
  }
  my %bondtypes_used;
  foreach my $bond (keys(%bonds)) {
    if (defined(${$bonds{$bond}}{type})) {
      if (defined($bondtypes_oldtonew{${$bonds{$bond}}{type}})) {
        ${$bonds{$bond}}{type} = $bondtypes_oldtonew{${$bonds{$bond}}{type}};
      }
      $bondtypes_used{${$bonds{$bond}}{type}} ++;
    }
  }
  foreach my $bondtype (sort {$a <=> $b} keys(%bond_types)) {
    if (!defined($bondtypes_used{$bondtype})) {
      delete($bond_types{$bondtype});
    }
  }
  undef(%bondtypes_oldtonew);
  my @bondtypesleft = sort {$a <=> $b} keys(%bond_types);
  my %newbondtypes;
  for (my $i = 0; $i < scalar(@bondtypesleft); $i++) {
    $bondtypes_oldtonew{$bondtypesleft[$i]} = $i+1;
    $newbondtypes{$i+1} = $bond_types{$bondtypesleft[$i]};
    delete($bond_types{$bondtypesleft[$i]});
  }
  %bond_types = %newbondtypes;
  foreach my $bond (keys(%bonds)) {
    if (defined(${$bonds{$bond}}{type})) {
      ${$bonds{$bond}}{type} = $bondtypes_oldtonew{${$bonds{$bond}}{type}};
    }
  }
  
  my %fingerprints_angle;
  my %angletypes_oldtonew;
  #xx print "604type ",${$angles{604}}{type},"\n";
  foreach my $angletype (sort {$a <=> $b} keys(%angle_types)) {
    my $fingerprint = ${$angle_types{$angletype}}{force_constant}.":".${$angle_types{$angletype}}{equil_value};
    if (defined($fingerprints_angle{$fingerprint})) {
      $angletypes_oldtonew{$angletype} = $fingerprints_angle{$fingerprint};
      delete($angle_types{$angletype});
    }
    else {
      $fingerprints_angle{$fingerprint} = $angletype;
      $angletypes_oldtonew{$angletype} = $angletype;
    }
  }
  #xx print join(":",values(%fingerprints_angle)),"\n";
  my %angletypes_used;
  foreach my $angle (keys(%angles)) {
    if (defined(${$angles{$angle}}{type})) {
      if (defined($angletypes_oldtonew{${$angles{$angle}}{type}})) {
        ${$angles{$angle}}{type} = $angletypes_oldtonew{${$angles{$angle}}{type}};
      }
      $angletypes_used{${$angles{$angle}}{type}} ++;
    }
  }
  foreach my $angletype (sort {$a <=> $b} keys(%angle_types)) {
    if (!defined($angletypes_used{$angletype})) {
      delete($angle_types{$angletype});
    }
  }
  undef(%angletypes_oldtonew);
  my @angletypesleft = sort {$a <=> $b} keys(%angle_types);
  my %newangletypes;
  for (my $i = 0; $i < scalar(@angletypesleft); $i++) {
    $angletypes_oldtonew{$angletypesleft[$i]} = $i+1;
    $newangletypes{$i+1} = $angle_types{$angletypesleft[$i]};
    delete($angle_types{$angletypesleft[$i]});
  }
  #print join(":",%angletypes_oldtonew),"\n";
  #xx foreach (sort {$a <=> $b } keys(%angletypes_oldtonew)) {print "$_:$angletypes_oldtonew{$_}\n";}
  %angle_types = %newangletypes;
  foreach my $angle (keys(%angles)) {
    if (defined(${$angles{$angle}}{type})) {
      ${$angles{$angle}}{type} = $angletypes_oldtonew{${$angles{$angle}}{type}};
    }
  }
  
  my %fingerprints_dihedral;
  my %dihedraltypes_oldtonew;
  foreach my $dihedraltype (sort {$a <=> $b} keys(%dihedral_types)) {
    my $fingerprint = ${$dihedral_types{$dihedraltype}}{force_constant}.":".${$dihedral_types{$dihedraltype}}{periodicity}.":".${$dihedral_types{$dihedraltype}}{phase};
    if (defined($fingerprints_dihedral{$fingerprint})) {
      $dihedraltypes_oldtonew{$dihedraltype} = $fingerprints_dihedral{$fingerprint};
      delete($dihedral_types{$dihedraltype});
    }
    else {
      $fingerprints_dihedral{$fingerprint} = $dihedraltype;
      $dihedraltypes_oldtonew{$dihedraltype} = $dihedraltype;
    }
  }
  my %dihedraltypes_used;
  foreach my $dihedral (keys(%dihedrals)) {
    if (defined(${$dihedrals{$dihedral}}{type})) {
      if (defined($dihedraltypes_oldtonew{${$dihedrals{$dihedral}}{type}})) {
        ${$dihedrals{$dihedral}}{type} = $dihedraltypes_oldtonew{${$dihedrals{$dihedral}}{type}};
      }
      $dihedraltypes_used{${$dihedrals{$dihedral}}{type}} ++;
    }
  }
  foreach my $dihedraltype (sort {$a <=> $b} keys(%dihedral_types)) {
    if (!defined($dihedraltypes_used{$dihedraltype})) {
      delete($dihedral_types{$dihedraltype});
    }
  }
  undef(%dihedraltypes_oldtonew);
  my @dihedraltypesleft = sort {$a <=> $b} keys(%dihedral_types);
  my %newdihedraltypes;
  for (my $i = 0; $i < scalar(@dihedraltypesleft); $i++) {
    $dihedraltypes_oldtonew{$dihedraltypesleft[$i]} = $i+1;
    $newdihedraltypes{$i+1} = $dihedral_types{$dihedraltypesleft[$i]};
    delete($dihedral_types{$dihedraltypesleft[$i]});
  }
  %dihedral_types = %newdihedraltypes;
  foreach my $dihedral (keys(%dihedrals)) {
    if (defined(${$dihedrals{$dihedral}}{type})) {
      ${$dihedrals{$dihedral}}{type} = $dihedraltypes_oldtonew{${$dihedrals{$dihedral}}{type}};
    }
  }
  }
  


  ${$_[0]}[0] = \%atoms;
  ${$_[0]}[1] = \%bonds;
  ${$_[0]}[2] = \%angles;
  ${$_[0]}[3] = \%dihedrals;
  ${$_[0]}[4] = \%bond_types;
  ${$_[0]}[5] = \%angle_types;
  ${$_[0]}[6] = \%dihedral_types;
  ${$_[0]}[7] = \%nonbonded_par;
  ${$_[0]}[8] = \%nonbonded_index;
  if (scalar(keys(%ljsigeps))) {${$_[0]}[9] = \%ljsigeps;}
  ${$_[0]}[10] = \$parmfilewas;
  if (scalar(keys(%pairs))) {${$_[0]}[11] = \%pairs;}
}

sub addGap {
  my @collected = @{$_[0]};
  my $gapstart = $_[1];
  my $gapsize = $_[2];
  my $atoms = $collected[0]; my $bonds = $collected[1];my $angles = $collected[2];my $dihedrals = $collected[3];my $bond_types = $collected[4];my $angle_types = $collected[5];my $dihedral_types = $collected[6];my $nonbonded_par = $collected[7];my $nonbonded_index = $collected[8];
  my %oldtonew;
  my %newtoold;
  #print scalar(keys(%{$atoms})), " atoms in state A\n";
  #print "gapstart $gapstart gapsize $gapsize\n";
  my $highest = (sort {$b <=> $a} keys(%{$atoms}))[0];
  for (my $i = 0; $i <= $highest; $i++) {
    if (defined(${$atoms}{$i})) {
      if ($i >= $gapstart) {
        $oldtonew{$i} = $i + $gapsize;
        $newtoold{$i+$gapsize} = $i;
      }
      else {
        $oldtonew{$i} = $i;
        $newtoold{$i} = $i;
      }
    }
  }
  foreach my $atom (sort {$b <=> $a} keys(%{$atoms})) {
    foreach my $excl (sort {$b <=> $a} keys(%{${${$atoms}{$atom}}{excluded}})) {
      #print "excl $excl\n";
      if (defined (${$atoms}{$excl})) {
        ${${${$atoms}{$atom}}{excluded}}{$oldtonew{$excl}} = delete ${${${$atoms}{$atom}}{excluded}}{$excl};
      }
    }
    ${$atoms}{$oldtonew{$atom}} = delete ${$atoms}{$atom};
  }

  foreach my $bond (keys(%{$bonds})) {
    my $newa1 = $oldtonew{${${$bonds}{$bond}}{a1}};
    ${${$bonds}{$bond}}{a1} = $newa1;
    my $newa2 = $oldtonew{${${$bonds}{$bond}}{a2}};
    ${${$bonds}{$bond}}{a2} = $newa2;
  }
  foreach my $angle (keys(%{$angles})) {
    my $newa1 = $oldtonew{${${$angles}{$angle}}{a1}};
    ${${$angles}{$angle}}{a1} = $newa1;
    my $newa2 = $oldtonew{${${$angles}{$angle}}{a2}};
    ${${$angles}{$angle}}{a2} = $newa2;
    my $newa3 = $oldtonew{${${$angles}{$angle}}{a3}};
    ${${$angles}{$angle}}{a3} = $newa3;
  }
  foreach my $dihedral (keys(%{$dihedrals})) {
    my $a1 = ${${$dihedrals}{$dihedral}}{a1};
    my $sign = $a1 / abs($a1);
    $a1 *= $sign;
    my $newa1 = $oldtonew{$a1};
    ${${$dihedrals}{$dihedral}}{a1} = $newa1 * $sign;
    my $a2 = ${${$dihedrals}{$dihedral}}{a2};
    $sign = $a2 / abs($a2);
    $a2 *= $sign;
    my $newa2 = $oldtonew{$a2};
    ${${$dihedrals}{$dihedral}}{a2} = $newa2 * $sign;
    my $a3 = ${${$dihedrals}{$dihedral}}{a3};
    $sign = $a3 / abs($a3);
    $a3 *= $sign;
    my $newa3 = $oldtonew{$a3};
    ${${$dihedrals}{$dihedral}}{a3} = $newa3 * $sign;
    my $a4 = ${${$dihedrals}{$dihedral}}{a4};
    $sign = $a4 / abs($a4);
    $a4 *= $sign;
    my $newa4 = $oldtonew{$a4};
    ${${$dihedrals}{$dihedral}}{a4} = $newa4 * $sign;
  }
}

sub trim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}



