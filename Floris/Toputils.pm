package Toputils;
use lib "/home/fbuelen/fec";

my $grompp = "grompp";
my $gmxdump = "gmxdump";

use strict;
use Exporter;  
our @ISA = qw(Exporter);
#our @EXPORT = qw(ParseParm7 ParseGtop WriteGTop WriteParm7 Cleanup addGap);
our @EXPORT = qw(ParseTopology WriteTopology VmdSel VmdWriteSubset addGap CoorFromFile WritePdb);
$| = 1;

use ParseParm7;
use ParseGmxdump;
use WriteGTop;
use WriteParm7;
use ToputilsShared;


# core data array: atoms bonds angles dihedrals bond_types angle_types dihedral_types nonbonded_par nonbonded_index 
# atoms: hash, key=id, values: 
#   - name
#   - charge
#   - mass 
#   - atom_type_index
#   - number_excluded_atoms
#   - amber_atom_type
#   - ishydrogen
#   - excluded: hash with keys=excluded atoms, values=1
#   - resid
#   - resname
# bonds: hash, key=index, values = hash with keys a1, a2, type
# angles: hash, key=index, values = hash with keys a1, a2, a3, type
# dihedrals: hash, key=index, values = hash with keys a1, a2, a3, a4, type
# bond_types: hash with key=index, values = hash with keys force_constant, equil_value
# angle_types: hash with key=index, values = hash with keys force_constant, equil_value
# dihedral_types: hash with key=index, values = hash with keys force_constant, periodicity, phase
# nonbonded_par: hash with key=index, values = acoef, bcoef
# nonbonded_index: hash with key=index, value = pointer into nonbonded_par
# nonbonded pairing (acoef, bcoef) as follows:
#   - get both atoms' atom_type_index (T_i + T_j)
#   - nonbonded_index[NTYPES*(T_i-1)+T_j] gives a pointer P_ij
#   - acoef and bcoef are ${$nonbonded_par{P_ij}}{acoef} and -{bcoef}


sub ParseTopology {
  my $name = $_[0];
  print "$name\n";
  my @return;
  if ($name =~ /\.parm7$/ ) {
    open(PARMFILE, "$_[0]") or die ("couldn't open parm file $_[0]\n");
    my $bak = $/;
    undef $/;
    my $parmfile = <PARMFILE>;
    close(PARMFILE);
    print "Read parameter file...\n";
    @return = ParseParm7($parmfile);
    $/ = $bak;
  }
  elsif ($name =~ /\.top$/ ) {
    print "Using $grompp and $gmxdump...\n";
    `cpp $_[0] > /dev/shm/temp.cpp.$$.top`;
    open(TOPFILE, "/dev/shm/temp.cpp.$$.top") or die ("couldn't open top file \n");
    my %molecule_atomncounts;
    my $molecule;
    my $natoms;
    my $directive = "";
    # write out a fake gro file for gmxdump - just needs the right number of atoms
    while (my $line = <TOPFILE>) {
      if ($line =~ /^\s*\[\s*(\S+)\s*\]\s*$/) {
        $directive = $1;
        if ($directive eq "moleculetype") {
          undef $molecule;
          while (! defined($molecule)) {
            my $line2 = <TOPFILE>;
            if ($line2 =~ /^([^;\s]+)\s+\d+/) {
              $molecule = $1;
            }
            if ($line2 =~ /^\[\s*(\S+)\s*\]$/) {
              die("failed to find molecule name in ParseToplogy\n");
            }
          }
          #print "molecule $molecule\n";
        }
      }
      if ($directive eq "atoms") {
        if ($line =~ /^\s*\d+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+/) {
          $molecule_atomncounts{$molecule} ++;
        }
      }
      if ($directive eq "molecules") {
        #print "molecule_atomncounts ",join(" ",%molecule_atomncounts),"\n";
        if ($line =~ /^([^;\s]+)\s+(\d+)/) {
          #print "line $line\n";
          if (!defined($molecule_atomncounts{$1})) {
            die("No molecule section found for $1\n");
          }
          $natoms += $2 * $molecule_atomncounts{$1} ;
        }
      }
    }
    close(TOPFILE);
    open(GROZERO,">/dev/shm/zero.$$.gro") or die ("couldn't open zero.gro\n");
    #print "$natoms atoms for grozero\n";
    print GROZERO "gro format file zero.$$.gro, garbage to feed to gmxdump\n";
    printf GROZERO "%6i\n",$natoms;
    for (my $i= 0; $i < $natoms; $i++) {
      print GROZERO "    1XXX      X    1   0.000   0.000   0.000\n";
    }
    print GROZERO " 100.00000 100.00000 100.00000\n";
    close(GROZERO);
    my $gromppout = `$grompp -f /home/fbuelen/fec/gmxdump.mdp -c /dev/shm/zero.$$.gro -p $name -o /dev/shm/gmxdump.$$.tpr -po /dev/shm/gmxdump.$$.mdp -maxwarn 100 2>&1`;
    if ($?) {print "$gromppout\n"; die();}
    my $gmxdumpout = `$gmxdump -sys -s /dev/shm/gmxdump.$$.tpr > /dev/shm/gmxdump.$$.mdp 2>&1`;
    if ($?) {print "$gmxdumpout\n"; die();}
    #`cp /dev/shm/gmxdump.$$.mdp /dev/shm/temp.$$.mdp`;
    #print `/home/fbuelen/fec/gmxdumpbugworkaround.pl /dev/shm/gmxdump.mdp`;
    #`/usr/local/gromacs/453/impi321-fftw322-gcc421/bin/grompp -f /home/fbuelen/fec/gmxdump.mdp -c /dev/shm/zero.gro -p $name -o /dev/shm/gmxdump.tpr -maxwarn 100`;
    #`/usr/local/gromacs/453/impi321-fftw322-gcc421/bin/gmxdump -sys -s /dev/shm/gmxdump.tpr > /dev/shm/gmxdump.mdp`;
    #open(GMXDUMPFILE, "/dev/shm/gmxdump.mdp") or die ("couldn't open gmxdump file /dev/shm/gmxdump.mdp\n");
    #undef $/;
    #my $gmxdumpfilein = <GMXDUMPFILE>;
    #close(GMXDUMPFILE);
    #@return = ParseGmxdump($gmxdumpfilein, "A");
    @return = ParseGmxdump("/dev/shm/gmxdump.$$.mdp", "A");
    `rm /dev/shm/gmxdump.$$.tpr`;
    `rm /dev/shm/gmxdump.$$.mdp`;
    `rm /dev/shm/zero.$$.gro`;
    `rm /dev/shm/temp.cpp.$$.top`;
  }
  else {
    die "File $name: ParseTopology accepts either amber parameter files (needs name ending .parm7) or gromacs topology files (needs name ending .top)\n";
  }
  
  return @return;
}

sub WriteTopology {
  my $AB = "";
  if (defined($_[2]) && ($_[2] eq "A")) {
    $AB = "A";
  }
  if (defined($_[2]) && ($_[2] eq "B")) {
    $AB = "B";
  }
  if ($_[0] =~ /\.parm7$/ ) {
    WriteParm7($_[0], $_[1], $AB);
  }
  elsif ($_[0] =~ /\.top$/ ) {
    WriteGTop($_[0], $_[1]);
  }
  else {
    die "WriteTopology writes either amber parameter files (needs name ending .parm7) or gromacs topology files (needs name ending .top)\n";
  }
}

sub VmdSel {
  my %args = %{$_[0]};
  my $seltext = $args{seltext};
  my $topol;
  my $coor;
  if (defined($args{coor})) {$coor = $args{coor};}
  if (defined($args{topol})) {$topol = $args{topol};}
  
  if (($seltext =~ /within/) && !defined($coor)) {
    die("VMD distance selection requires (real) coordinates to be passed\n");
  }
  my $parmfile;
  if (($topol =~/\.parm7$/) || ($topol =~/\.pdb$/) || ($topol =~/\.gro$/)) {
    $parmfile = $topol;
  }
  elsif ($topol =~/\.top$/) {
    my @collected = ParseTopology($topol);
    $parmfile = "/dev/shm/temp.$$.parm7";
    WriteTopology("/dev/shm/temp.$$.parm7", \@collected);
  }
  elsif (scalar(@{$topol})) {
    $parmfile = "/dev/shm/temp.$$.parm7";
    WriteTopology("/dev/shm/temp.$$.parm7", $topol);
  }

  if (!defined($topol)) {
    if (($coor =~/\.pdb$/) || ($coor =~/\.gro$/)) {
      $parmfile = $coor;
    }
    else {
      print "inadequate coordinates / topology supplied to VmdSel\n";
    }
  }
  
  my $vmdcmd = "";
  $vmdcmd .= "mol new $parmfile\n";
  if (defined($coor)) {$vmdcmd .= "mol addfile $coor\n";}
  $vmdcmd .= "set fp1 [open \"sel.txt\" \"w\"]\n";
  $vmdcmd .= "set sel [atomselect top \"$seltext\"]\n";
  $vmdcmd .= "puts \$fp1 [\$sel get serial]\n";
  $vmdcmd .= "close \$fp1\n";
  
  #print $vmdcmd;
  open(VMDIN,">vmdin");
  print VMDIN $vmdcmd;
  close(VMDIN);
  `cat vmdin | ~/vmdtext `;
  if ($parmfile eq "/dev/shm/temp.$$.parm7") {
    `rm $parmfile`;
  }

  my @sel = split(" ",`cat sel.txt`);
  

  if (!scalar(@sel)) {
    print "Atom selection $seltext returned no atoms ($args{topol} $args{coor})\n";
  }
  `rm vmdin`;
  `rm sel.txt`;
  return @sel;
}

sub VmdWriteSubset {
  my %args = %{$_[0]};
  my $coorout = $args{coorout};
  my $topol;
  my $coorin;
  if (defined($args{coorin})) {$coorin = $args{coorin};}
  if (defined($args{topol})) {$topol = $args{topol};}
  
  my $parmfile;
  if (($topol =~/\.parm7$/) || ($topol =~/\.pdb$/) || ($topol =~/\.gro$/)) {
    $parmfile = $topol;
  }
  elsif ($topol =~/\.top$/) {
    my @collected = ParseTopology($topol);
    $parmfile = "/dev/shm/temp.$$.parm7";
    WriteTopology("/dev/shm/temp.$$.parm7", \@collected);
  }
  elsif (scalar(@{$topol})) {
    $parmfile = "/dev/shm/temp.$$.parm7";
    WriteTopology("/dev/shm/temp.$$.parm7", $topol);
  }

  if (!defined($topol)) {
    if (($coorin =~/\.pdb$/) || ($coorin =~/\.gro$/)) {
      $parmfile = $coorin;
    }
    else {
      print "inadequate coordinates / topology supplied to VmdWriteSubset\n";
    }
  }

  my %setbeta;
  if (defined($args{setbeta})) {
    %setbeta = %{$args{setbeta}};
  }

  my @writeout;
  if (ref($args{atoms}) eq "ARRAY") {
    @writeout = @{$args{atoms}};
  }
  else {
    @writeout = VmdSel({"seltext", $args{atoms}, "topol", $parmfile, "coor", $coorin});
  }



  my $vmdcmd = "";
  if (defined($parmfile)) {
    $vmdcmd .= "mol new $parmfile\n";
    $vmdcmd .= "mol addfile $coorin\n";
  }
  else {
    $vmdcmd .= "mol new $coorin\n";
  }
  foreach my $beta (sort {$b <=> $a} keys(%setbeta)) {
    my @arr = @{$setbeta{$beta}};
    $vmdcmd .= "[atomselect top \"serial ".join(" ",@arr)."\"] set beta $beta\n";
  }
  $vmdcmd .= "[atomselect top \"serial ".join(" ",@writeout)."\"] writepdb $coorout\n";
  
  #print $vmdcmd;
  open(VMDIN,">vmdin");
  print VMDIN $vmdcmd;
  close(VMDIN);
  `cat vmdin | ~/vmdtext `;
  `rm vmdin`;
  my $pdblength = `wc -l $coorout\n`;
  if (!$pdblength) {
    print "VmdWriteSubset command did not write a pdb file ($coorout)\n";
  }
}

sub CoorFromFile {
  my $atoms = $_[0];
  my $file = $_[1];
  my $lab = $file;
  $lab =~ s/\///g;
  if (! (-e $file)) {
    die("CoorFromFile: file $file does not exist\n");
  }
  my $vmdin;
  $vmdin = "mol new $file\n";
  $vmdin .= "set fp [open \"/dev/shm/vmd.$$.$lab.tmp\" \"w\"]\n";
  $vmdin .= "puts \$fp [[atomselect top \"all\"] num]\n";
  $vmdin .= "puts \$fp [[atomselect top \"all\"] get {x y z}]\n";
  $vmdin .= "close \$fp\n\n";
  open(VMDIN,">/dev/shm/vmdin.$$.$lab") or die("couldn't open /dev/shm/vmdin.$$.$lab.$lab\n");
  print VMDIN $vmdin;
  close(VMDIN);
  print "vmd\n";
  `cat /dev/shm/vmdin.$$.$lab | ~/vmdtext `;
  print "vmd done\n";

  if (! (-e "/dev/shm/vmd.$$.$lab.tmp")) {
    sleep(1);
    `cat /dev/shm/vmdin.$$.$lab | ~/vmdtext `;
    if (! (-e "/dev/shm/vmd.$$.$lab.tmp")) {
      die("/dev/shm/vmd.$$.$lab.tmp didn't show up yet\n");
    }
  }
  my $coor = `cat /dev/shm/vmd.$$.$lab.tmp`;
  $coor =~ /(\d+)/g;
  if ($1 != scalar(keys(%{$atoms}))) {
    print "CoorFromFile: file $file has different number of atoms ($1) than passed topology (",scalar(keys(%{$atoms})),")\n";
    die();
  }
  my $acount = 1;
  my $xyzcount = 1;
  while ($coor =~ /([-+]?[0-9]*\.?[0-9]+)/g) {
    if ($xyzcount == 1) { ${${$atoms}{$acount}}{x} = $1; $xyzcount = 2;}
    elsif ($xyzcount == 2) { ${${$atoms}{$acount}}{y} = $1; $xyzcount = 3;}
    elsif ($xyzcount == 3) { ${${$atoms}{$acount}}{z} = $1; $xyzcount = 1; $acount ++;}
  }
  `rm  /dev/shm/vmd.$$.$lab.tmp`;
  `rm /dev/shm/vmdin.$$.$lab`;
}

sub WritePdb {
  my $outfile = $_[0];
  my %atoms = %{$_[1]};
  open(PDBFILE, ">$outfile") or die ("couldn't open $outfile\n");
  if (defined($_[2])) {
    print PDBFILE $_[2];
  }
  my $acount = 0;
  foreach my $atom (sort {$a <=> $b} keys(%atoms)) {
    my $line = "ATOM  ";
    $acount ++;
    $line .= sprintf("%5i ",$acount);
    $line .= sprintf("%4s ", substr(${$atoms{$atom}}{name},0,4));
    $line .=  sprintf("%3s X", substr(${$atoms{$atom}}{resname},0,3));
    my $resnum = ${$atoms{$atom}}{resid};
    if ($resnum > 9999) {$resnum = 5000 + ($resnum % 5000);}
    $line .= sprintf("%4i    ",$resnum);
    $line .= sprintf(" %7.3f",${$atoms{$atom}}{x});
    $line .= sprintf(" %7.3f",${$atoms{$atom}}{y});
    $line .= sprintf(" %7.3f",${$atoms{$atom}}{z});
    $line .= "  0.00  0.00\n";
    print PDBFILE $line;
  }
  close(PDBFILE);
}


1;
