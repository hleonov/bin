#!/usr/bin/perl
# AUTHOR: Cameron Mura (07/2001). 
# LAST MODIFIED: 05/21/2004
# USAGE: 
# 	map_dataervation_to_B.pl file1.pdb file2_sites.txt
# 
# DESCRIPTION:
# This Perl script takes 2 input files: argv[0] ("file1.pdb") is the PDB file 
# and argv[1] ("file2_sites.txt") is a file tabulating dataerved residues in 
# the format:
#
# "res#,ss(super strong)|s(trong)|m(edium)|w(eak)" ...
#
# where "res#" is the residue # and "superstrong", "strong", "medium", 
# or "weak" specify how strongly that site is dataerved.
# 
# The output is a PDB file with all of the B-factors flattened to 20.00,
# except for the "superstrong", "strong", "medium", or "weak" sites, which
# are assigned B-values of 90.00, 70.00, 40.00, or 33.00, respectively. 
# NOTE: make sure all occupancies are 1.00 for any atoms with B-factors you
# want changed.
#
# This is useful for programs like GRASP, which can color a surface by the B-factor 
# field of PDB file.  Or Robert Campbell's color_b.py module for PyMOL.

$pdb_in = $ARGV[0];
$site_file = $ARGV[1];
$bfac = "";
@data = "";
$site_line = "";

my %HoH;
open (SITES, $site_file) || die "Cannot open file \"$site_file\"\n";
while (<SITES>) {
   $site_line = $_; 
   chomp ($site_line); 
   @data = "";
   @data = split(/\s+/, $site_line);
   $HoH{$data[0]}{$data[1]} = $data[2];   
   #print "HoH{$data[0]}{$data[1]} = $data[2]\n"
}

close (SITES);

open (PDB, $pdb_in) || die "Cannot open file \"$pdb_in\"\n";

while (<PDB>)
	{
 	 $line = $_;	
	 chomp ($line); 
	 $resnum = ""; 
    $bfac = "";
    #if ($line =~ m/^(ATOM|HET)\s+\d+\s+[\w\*]+\s+\w\w\w\s+\w\s+/) {
    #print $line;
    #exit;
	 if ($line =~ m/(^(ATOM|HET)\s+\d+\s+([\w\*])+\s+\w\w\w\s+(\w+)\s+)(.+)(1\.00\s+)(\d+\.\d+)(.*)/) {
	   
       #printf "1 $1\n2 $2\n3 $3\n4 $4\n5 $5\n6 $6\n7 $7\n8 $8\n";
       #print $line,"\n";
       $chain = $3;     
	    $resnum = $4; 
       #$bfac = "20.00";
		 if (exists $HoH{$resnum}{$chain}) {
          $bfac = $HoH{$resnum}{$chain};
          #if ($HoH{$resnum}{$chain} == 1) {
#              $bfac = "100.00";               
#          } elsif ($HoH{$resnum}{$chain} == 0) {
#              $bfac = "0.00";
#          }
       }
	      		
	    #if ($bfac eq "") { $bfac = "20.00";}
	    #$shit = $5 . $bfac;
   	    #print "$1$3$4$shit$7\n";
       print "$1$5$6$bfac$8\n";
	   }
	}

