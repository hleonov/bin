#!/usr/bin/perl -w
my $scripdir = ".";
BEGIN {
  $scriptdir = $0;
  $scriptdir =~ s/\/[^\/]*$//;
}
use lib "$scriptdir";
use strict;
use lib ".";
use Toputils;
use Getopt::Long;
$|=1;

my (@in, @out);
my (@interestingsel, @bufferdist);
if (! GetOptions('i=s' => \@in, 'o=s' => \@out, 's=s', \@interestingsel, 'b=f', \@bufferdist)) {
  die("Usage:\n  -i [input files] topology.[parm7|top] and coordinates.[gro|pdb|rst7]\n  -o [output files] topology.[parm7|top] coordinates.[pdb] posres.[itp]\n  -s vmd selection for interesting zone\n  -b buffer distance\n");
}

if (!(scalar(@in) && scalar(@out) && scalar(@interestingsel) && scalar(@bufferdist))) {
  die("Usage:\n  -i [input files] topology.[parm7|top] and coordinates.[gro|pdb|rst7]\n  -o [output files] topology.[parm7|top] coordinates.[pdb] posres.[itp]\n  -s vmd selection for interesting zone\n  -b buffer distance\n");
}

if (scalar(@interestingsel) != scalar(@bufferdist)) {
  die("Multiple interesting regions are supported but need one selection text and one buffer distance for each\n");
}

#if (scalar(@ARGV) == 0) {
#  die("Arugments: 1) topology file in (top or parm7), 2) topology file out (top or parm7), 3) coordinates in, 4) coordinates out (trimmed), 5) position restraints itp, 6) vmd selection for \"interesting\" zone, 7) distance in A for buffer zone\n");
#}

my ($topin, $coorin, $coorout, $itpout);
my @topout; # allow multiple output topologies to be written
foreach my $infile (@in) {
  if ($infile =~ /\.parm7$/) {
    if (defined($topin)) {die("double defition of input topology\n");}
    $topin = $infile;
  }
  elsif ($infile =~ /\.top$/) {
    if (defined($topin)) {die("double defition of input topology\n");}
    $topin = $infile;
  }
  if ($infile =~ /\.gro$/) {
    if (defined($topin)) {die("double defition of input coordinates\n");}
    $coorin = $infile;
  }
  elsif ($infile =~ /\.pdb$/) {
    if (defined($coorin)) {die("double defition of input coordinates\n");}
    $coorin = $infile;
  }
  elsif ($infile =~ /\.rst7$/) {
    if (defined($coorin)) {die("double defition of input coordinates\n");}
    $coorin = $infile;
  }
}
foreach my $outfile (@out) {
  if ($outfile =~ /\.parm7$/) {
    push @topout, $outfile;
  }
  elsif ($outfile =~ /\.top$/) {
    push @topout, $outfile;
  }
  if ($outfile =~ /\.pdb$/) {
    if (defined($coorout)) {die("double defition of output coordinates\n");}
    $coorout = $outfile;
  }
  if ($outfile =~ /\.itp$/) {
    if (defined($itpout)) {die("double defition of output position restraints\n");}
    $itpout = $outfile;
  }
}

if (!defined($itpout)) {die("no output file specified for position restraints\n");}
if (!scalar(@topout)) {die("no output file specified for topology\n");}
if (!defined($coorout)) {die("no output file specified for coordinates\n");}
if (!defined($topin)) {die("no input file specified for topology\n");}
if (!defined($coorin)) {die("no input file specified for coordinates\n");}
if (!scalar(@bufferdist)) {die("no buffer distance specified\n");}
if (!scalar(@interestingsel)) {die("no atom selection for interesting zone specified\n");}


my %interesting;
my %buffer;
my @keepsel;
my @interestingatoms;
my @bufferatoms;

my $countsels;

for (my $i = 0; $i < scalar(@interestingsel); $i++) {
  my $sel = $interestingsel[$i];
  my $dist = $bufferdist[$i];
  @interestingatoms = (@interestingatoms, VmdSel({"seltext", $sel, "topol", $topin, "coor", $coorin}));
  my $buffersel = "(same resid as (within $dist of ($sel))) and not ($sel)\n";
  @bufferatoms = (@bufferatoms, VmdSel({"seltext", $buffersel, "topol", $topin, "coor", $coorin}));
  @keepsel = (@keepsel, @interestingatoms, @bufferatoms);
  $countsels += scalar(@interestingatoms) + scalar(@bufferatoms);
  %interesting = (%interesting, map { $_ => 1} @interestingatoms);
  %buffer = (%buffer, map { $_ => 1} @bufferatoms);
}

if ($countsels != scalar(@keepsel)) {
  die("multiple overlapping interesting / buffer regions\n");
}

my @selbeta33 = VmdSel({"seltext", "(serial ".join(" ",@interestingatoms).") or ((serial ".join(" ",@bufferatoms).") and (water or numbonds 0))", "topol", $topin, "coor", $coorin});
my @selbeta66 = VmdSel({"seltext", "(serial ".join(" ",@bufferatoms).") and not (water or numbonds 0)", "topol", $topin, "coor", $coorin});

#my %setbeta = (33, \@interestingatoms, 66, \@bufferatoms);

my %setbeta = (33, \@selbeta33, 66, \@selbeta66);

VmdWriteSubset({"atoms", \@keepsel, "topol", $topin, "coorin", $coorin, "coorout", $coorout, "setbeta", \%setbeta});

print scalar(keys(%interesting)), " atoms in interesting region, ",scalar(keys(%buffer))," in buffer region\n";
if (!scalar(keys(%interesting))) {die("no atoms selected for interesting region\n");}
my %merged = (%interesting, %buffer);

#print "\n",join(" ",sort {$a <=> $b} keys(%merged)),"\n\n";

my @return = ParseTopology($topin);

my %atoms = %{$return[0]};
my %bonds = %{$return[1]};
my %angles = %{$return[2]};
my %dihedrals = %{$return[3]};
my %bond_types = %{$return[4]};
my %angle_types = %{$return[5]};
my %dihedral_types = %{$return[6]};
my %nonbonded_par = %{$return[7]};
my %nonbonded_index = %{$return[8]};
my %ljsigeps = %{$return[9]};

print "Atoms: ",scalar(keys(%atoms)),"\n";
print "Bonds: ",scalar(keys(%bonds)),"\n";
print "Angles: ",scalar(keys(%angles)),"\n";
print "Dihedrals: ",scalar(keys(%dihedrals)),"\n";
print "Bond types: ",scalar(keys(%bond_types)),"\n";
print "Angle types: ",scalar(keys(%angle_types)),"\n";
print "Dihedral types: ",scalar(keys(%dihedral_types)),"\n";
print "LJ types: ",scalar(keys(%ljsigeps)),"\n";

foreach my $atom (keys(%atoms)) {
  if (! defined($merged{$atom})) {
    delete($atoms{$atom});
  }
  else {
    foreach my $excluded (keys(%{${$atoms{$atom}}{excluded}})) {
      if (!defined($merged{$excluded})) {
        delete ${${$atoms{$atom}}{excluded}}{$excluded};
        ${$atoms{$atom}}{number_excluded_atoms} --;
      }
    }
  }
}
foreach my $bondnr (keys(%bonds)) {
  my %bond = %{$bonds{$bondnr}};
  if (! defined($merged{$bond{a1}})) {delete($bonds{$bondnr});}
  elsif (! defined($merged{$bond{a2}})) {delete($bonds{$bondnr});}
}
foreach my $anglenr (keys(%angles)) {
  my %angle = %{$angles{$anglenr}};
  if (!defined($merged{$angle{a1}})) {delete($angles{$anglenr});}
  elsif (!defined($merged{$angle{a2}})) {delete($angles{$anglenr});}
  elsif (!defined($merged{$angle{a3}})) {delete($angles{$anglenr});}
}
foreach my $dihedralnr (sort {$a <=> $b} keys(%dihedrals)) {
  my %dihedral = %{$dihedrals{$dihedralnr}};
  if (!defined($merged{abs($dihedral{a1})})) {delete($dihedrals{$dihedralnr});}
  elsif (!defined($merged{abs($dihedral{a2})})) {delete($dihedrals{$dihedralnr});}
  elsif (!defined($merged{abs($dihedral{a3})})) {delete($dihedrals{$dihedralnr});}
  elsif (!defined($merged{abs($dihedral{a4})})) {delete($dihedrals{$dihedralnr});}
}

# correct residue numbers
my $previous = 0;
my $offset=0;
foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
  my $resid = ${$atoms{$atom}}{resid};
  if ($resid - $previous > 1) {
    $offset += $resid - $previous - 1;
  }
  ${$atoms{$atom}}{resid} -= $offset;
  $previous = $resid;
}

my $restrsel = "beta 66 and (not water) and (not numbonds 0)";
my @restratoms;
my $parm7ifpossible = $topout[0];
foreach my $tout (@topout) {
  if ($tout =~ /\.parm7$/) {
    $parm7ifpossible = $tout;
  }
}
@restratoms = VmdSel({"seltext", $restrsel, "topol", $parm7ifpossible, "coor", $coorout});

foreach my $atom (@selbeta66) {
  if (!defined($atoms{$atom})) {
    die("trying to restrain a non-existent atom $atom\n");
  }
  ${$atoms{$atom}}{restraint_k} = 10;
}



print "Writing parameter file...\n";
print "Atoms: ",scalar(keys(%atoms)),"\n";
print "Bonds: ",scalar(keys(%bonds)),"\n";
print "Angles: ",scalar(keys(%angles)),"\n";
print "Dihedrals: ",scalar(keys(%dihedrals)),"\n";
print "Bond types: ",scalar(keys(%bond_types)),"\n";
print "Angle types: ",scalar(keys(%angle_types)),"\n";
print "Dihedral types: ",scalar(keys(%dihedral_types)),"\n";
print "LJ types: ",scalar(keys(%ljsigeps)),"\n";
my @collected;
push @collected, \%atoms, \%bonds, \%angles, \%dihedrals, \%bond_types, \%angle_types, \%dihedral_types, \%nonbonded_par, \%nonbonded_index, \%ljsigeps, "trimsystem.pl";

foreach my $tout (@topout) {
  WriteTopology($tout, \@collected);
}



#open (ITPFILE, ">$itpout");
#print ITPFILE "; Position restraints for buffer zone from trimsystem.pl\n";
#print ITPFILE "; 8368 corresponds to 10 kcal/mol/A^2\n\n";
#print ITPFILE "[ position_restraints ]\n; atom  type      fx      fy      fz\n";
#foreach my $atom (sort {$a <=> $b} @restratoms) {
#  printf(ITPFILE "%6i     1  8368  8368  8368\n",$atom);
#}

#close(ITPFILE);

