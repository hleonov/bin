#!/usr/bin/perl

#
# This program utilizes VMD to add atoms to a list of PDB structures.
#
# The only inputs you need are a 1) a file containing a list of the
# PDB files you want to guess atoms for (see variable '$PDB_list' below)
# and you will need the properly formatted create pdb script. The
# script is given below. Place this script in a file called "create_pdb_script.tcl"
#
# To execute this program: >perl Add_Atoms_to_PDBs.pl
#
# Ed O'Brien Jr. August 11, 2004, Univ. Maryland CP, Dave Thirumalai's Group
#
#

#
# Here is the .tcl script you'll need, place it AS IS in a file called
# "create_pdb_script.tcl". DO NOT CHANGE this tcl script. 
# Note: DO NOT include '=pod', '=cut' in the file
#
=pod
package require psfgen
topology top_all27_prot_lipid.inp
        alias residue HIS HSD
        alias residue HOH TIP3
        alias residue ZN ZN2
        alias atom ILE CD1 CD
        alias atom GLY OXT OT1
        alias atom LEU OXT OT1
        alias atom PHE OXT OT1
        alias atom HOH O OH2
                                                                                                                                                             
        segment A {
        pdb temp.pdb
        first none
        last none
        }
                                                                                                                                                             
        coordpdb temp.pdb A
        guesscoord
        writepdb pdb1j0p.ent.ADD
        exit
=cut

####### INPUTS: You should only need to modify the next line ##############

$PDB_list = "Final_Master_list.txt";

###########################################################################

###### You shouldn't need to modify anything after this line ##############

$TCL_script = "create_pdb_script.tcl";

open(IN,"$PDB_list");
@List = <IN>;
close IN;

$nfiles = 0;
for $name (@List)
  {
  $nfiles = $nfiles + 1;
  # Remove white spaces at begining of string
  $name =~ s/^\s+//;
  # Remove white spaces at end of string
  $name =~ s/\s+$//;
  print "$nfiles: Adding atoms to $name\n";

  # Select ONLY the protein from the PDB
  open(TCL,">select_protein.tcl");
    print TCL "mol load pdb $name\n";
    print TCL "set protein [atomselect top protein]\n";
    print TCL '$protein writepdb temp.pdb';
    print TCL "\nexit\n";
  close TCL;   

  system "rm -f temp.pdb";
  system "rm -f junk.txt";
  system "vmd -dispdev text -e select_protein.tcl > junk.txt";

  # Modify the TCL script to submit proper PDB
  open(TCL,"$TCL_script");
  @TCL2 = <TCL>;
  close TCL;

  # Print the Modified TCL script out
  $output_name = "$name" . ".ADD";
  open(TCL,">$TCL_script");
  $l = 0;
  for $line (@TCL2)
    {  
    $l = $l + 1;
    if( $l == 20 ) { print TCL "        writepdb $output_name\n" }
    else { print TCL "$line" }
    }
  close TCL;

  # Now run the TCL script by calling VMD command line
  system "rm -f junk.txt";
  system "vmd -dispdev text -e $TCL_script > junk.txt";
  open(TEST,"$output_name");
  @TEST2 = <TEST>;
  close TEST;
  if( $#TEST2 == -1 ) { print "***** Writing coords to file FAILED! ******\n" }
  print "Finished: $name\n";
  }
print "$nfiles completed\n";
