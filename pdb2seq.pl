#!/usr/local/bin/perl
# Usage: perl sequence.pl <pdb file>
# takes a *.pdb file, and outputs the amino acid sequence in a one letter format,
# seperating different chains.
$LINEMAX = 50;
%aminoAcids = (ALA => A, VAL => V, PHE => F, PRO => P, MET => M, ILE => I, 
		   LEU => L, ASP => D, GLU => E, LYS => K, ARG => R, SER => S, 
		   THR => T, TYR => Y, HIS => H, CYS => C, ASN => N, GLN => Q,
		   TRP => W, GLY => G,MSE => M);

if ($#ARGV == -1) {
	print "Usage: perl sequence.pl <pdb file>\n";
   	exit(0);
}
$pdbFile = shift;
open(PDB, $pdbFile) or die ("$pdbFile not found\n");
print "Input File is: $pdbFile\n";

$outfile = $1 if ($pdbFile =~ m/^(.+)\.\w+$/);
open(OUT, ">$outfile.seq");
print "OutputFile is: $outfile.seq\n";

$prevResNum = 0;
$sequence = "";
while ($line = <PDB>) {
    #ATOM    266  N   PHE L  34      58.914 108.628   7.378  1.00 39.60 
	#\s*\d*\s*\w*\s*          N  
    if ($line =~ m/^ATOM.{13}(\w{3})\s*(\w)\s*(\d*)/) {
	my $residue = $1;
	my $chain = $2;
	my $residueNumber = $3;
	
	if (!defined $prevChain) {   #first chain
	    $prevChain = $chain ;
	    $count=1;
	    print OUT ">$outfile-chain-$chain-$residueNumber\n";
	}
	
	
	if ($prevChain ne $chain) {   #new chain
	    printSeq($sequence);
	    print OUT "\n>$outfile-chain-$chain-$residueNumber\n";
	    $prevChain = $chain;
	    $prevResNum = 0;
	    $count=1;
	    $sequence = "";
	}

		if ($prevResNum != $residueNumber) {  #middle of chain, new residue.
		    my $diff = abs($residueNumber-$prevResNum)-1 if ($prevResNum != 0);
		    if ($diff > 0) {
			for (my $j=1; $j<=$diff; $j++) {
			    $sequence .= 'X';
			}
			print "diff: $diff  $prevResNum $residueNumber\n" ;
		    }
		    $prevResNum = $residueNumber ;
		    
		    $sequence .= $aminoAcids{$residue};		    
		  
		}
		
	    } 
    }
printSeq($sequence);
close(OUT);

sub printSeq {
    @seq = split(//,$_[0]);
    $count = 1;
    for (my $i=0; $i<=$#seq; $i++) {
	if ($count <= $LINEMAX) {
	    print OUT "$seq[$i]";
	} else {
	    print OUT "\n$seq[$i]";
	    $count=1;
	}
	$count++;
    }
}
