#!/usr/bin/perl -w
use Getopt::Long;

undef(%options);
GetOptions ('i=s' => \$input, 
			'o=s' => \$output,
			'cut' => \$options{'cut'},
			'chain=s' => \$options{'chain'},
			'term=s' => \$options{'term'},
			'solvent=s' => \$options{'solvent'},
			'rename' => \$options{'rename'},
			'from_res=s' => \$options{'from_res'},
			'to_res=s' => \$options{'to_res'},
			'het2atom' => \$options{'het2atom'});
&usage if (!defined $input || !defined $output); 

&init;
my $pdb = read_input($input);
my @outarr;
# cut protein - keep only the protein in a new file
if (defined $options{'cut'}) {
	@outarr = @{cut_protein($pdb)};
# add ACE and NAC to terminus - by atom number from original system
} elsif (defined $options{'term'}) {
	my $term_pdb = read_input($options{'term'});
	@outarr = @{add_terminus($pdb, $term_pdb)};
# edit chains - rename chains to A,B,C,D - change @chains array in init().
# only works for homo-oligomers for now.
} elsif (defined $options{'chain'}) {
	@outarr = @{edit_chains($pdb,$options{'chain'})};
# create a new file of the whole system (DMPC, SOL, Na, protein)
} elsif (defined $options{'solvent'}) {
	my $sol_pdb = read_input($options{'solvent'});
	@outarr = @{add_solvent($pdb, $sol_pdb)};
} elsif (defined $options{'rename'}) {
	&usage if (!defined $options{'from_res'} || !defined $options{'to_res'});
	@outarr = @{replace_resn($pdb, $options{'from_res'}, $options{'to_res'})}; 
} elsif (defined $options{'het2atom'}) {
	@outarr = @{het2atom($pdb)};
} else {&usage}

print "input: $input\toutput: $output\n";

open(OUT,">$output");
print OUT @outarr;
close(OUT);
#==========================================================================

#het2atom{pdbarr_ref}
sub het2atom {
	my @pdb = @{$_[0]};
	for (my $i=0; $i<=$#pdb; $i++) {
		$pdb[$i] =~ s/HETATM/ATOM  /;
		push(@new, $pdb[$i]); 
	}
	return \@new;
}
#replace_resn{pdbarr_ref, from, to)
sub replace_resn {
	my @pdb = @{$_[0]};
	for (my $i=0; $i<=$#pdb; $i++) {
		$pdb[$i] =~ s/$_[1]/$_[2]/;
		push(@new, $pdb[$i]); 
	}
	return \@new;
}

# add_solvent($pdbarr_ref, solvent_arr_ref)
# Add every atom 
sub add_solvent {
	my @pdb = @{$_[0]};
	my @s_pdb = @{$_[1]};
	for (my $i=0; $i<=$#pdb; $i++) {
		if (($pdb[$i] =~ m/^ATOM\s+(\d+)/) ||
			($pdb[$i] !~ m/(END)|(TER)/)) {
			push(@new,$pdb[$i]);
		}
	}
	for (my $i=0; $i<=$#s_pdb; $i++) {
		if ($s_pdb[$i] =~ m/^ATOM\s+\S+\s+\S+\s+(\S+)/) {
			my $aa = $1;
			if (length($aa) >4) {
				$aa = substr ($aa,0,4);
			}
			if (defined $solvent{$aa}) {
				push(@new,$s_pdb[$i]);
			}
		}
	}
	push (@new,"TER\n");
	return \@new;
}

#==========================================================================
# edit_chains($pdbarr_ref)
# edit change according to ACE seperator
sub edit_chains {
	my @pdb = @{$_[0]};
	my $last_num = $_[1];
	my $c = -1;
	my $num = 0;	# how many residues we encountered.
	my $prev = "";	#name of previous residue
	for (my $i=0; $i<=$#pdb; $i++) {
		if ($pdb[$i] =~ m/^(ATOM.{13}\w{3})..(.{4})(.*)$/) {
			my $begin = $1; 
			my $residue = $2;
			my $end = $3;
			#Found a new residue!
			if ($residue ne $prev) {
				$num++;
				#first residue of new chain
				if ($num % $last_num == 1) {	
					$c++;
				} 
			#	print "$begin $residue $num $chains[$c]\n";
				$prev=$residue;
			}
			push(@new,$begin." ".$chains[$c].$residue.$end."\n");
		} else {
			push(@new, $pdb[$i]);
		}
	}
	return \@new;
	
}
#==========================================================================
# add_terminus($pdbarr_ref, $terminus_arr_ref)
# search chain ends and insert the ACE/NAC in between.
sub add_terminus {
	my @pdb = @{$_[0]};
	my @t_pdb = @{$_[1]};
	my %ace;
	my %nac;
	#store terminus in hashes by chain
	for (my $i=0; $i<=$#t_pdb; $i++) {
		if ($t_pdb[$i] =~ m/^ATOM.{13}ACE\s(.)/) {
			$ace{$1} .= $t_pdb[$i];
		} elsif ($t_pdb[$i] =~ m/^ATOM.{13}NAC\s(.)\s+\d+/) {
			$nac{$1} .= $t_pdb[$i];
		}
	}

	my $prev_ch = "AA";
	$nac{$prev_ch} = ""; #fictious nac for beginning
	my $chain;
	for (my $i=0; $i<=$#pdb; $i++) {
		if ($pdb[$i] =~ m/^ATOM.{17}(.)/) {
			my $chain = $1;
		 	#end prev chain, start new chain
			if ($prev_ch ne $chain) {
				push(@new,$nac{$prev_ch});
				push(@new,$ace{$chain});
			}
			push(@new,$pdb[$i]);
			$prev_ch = $chain;
		} elsif ($pdb[$i] !~ m/(END)|(TER)/) {
			push(@new,$pdb[$i]); #header and remarks - keep as are
		} 
	}
	push(@new,$nac{$prev_ch}); #last chain nac
	push (@new,"TER\nENDMDL\n");
	return \@new;
}

#==========================================================================
# cut_protein($pdbarr_ref)
sub cut_protein {
	my @pdb = @{$_[0]};
	my @new;
	for (my $i=0; $i<=$#pdb; $i++) {
		if ($pdb[$i] =~ m/^ATOM\s+\S+\s+\S+\s+(\S+)/) {
			my $aa = $1;
			if (defined $aa_rep{$aa}) {
				push(@new,$pdb[$i]);
			}
		} else {
			push(@new,$pdb[$i]);
		}
	}
	return \@new;
}

#==========================================================================
sub read_input {
	open(IN,$_[0]);
	my @in = <IN>;
	close(IN);
	return \@in;
}

sub init {
	# a hash of amino acids in 3-letter format. values are ges scale.
	%aa_rep = (ALA => 1.6, CYS => 2.0, ASP => -9.2, GLU => -8.2, PHE => 3.7,
              GLY => 1.0, HIS => -3.0,ILE => 3.1,  LYS => -8.8, LEU => 2.8,
         	  MET => 3.4, ASN => -4.8,PRO => -0.2, GLN => -4.1, ARG => -12.3,
           	  SER => 0.6, THR => 1.2, VAL => 2.6,  TRP => 1.9,  TYR => -0.7,
              HID => -3.0, HIP => -3.0, HIE => -3.0, LYP => -8.8);
	@chains = ('A','B','C','D','E','F','G','H','I','J','K','L');
	%solvent = (SOL => 1, 'NA+' => 1, 'Na+' => 1, Na => 1, Ca => 1, 
				'CL-' => 1, 'Cl-' => 1, Cl => 1, 
				DMP => 1, DMPC => 1, HOH => 1, TIP3 => 1, TIP => 1);
}

sub usage {
	print "Usage: perl $0 -i infile -o outfile [options]\n";
	print "Options: -cut : cut protein from a whole system\n";
	print "\t -chain <last resid>: Edit chain identifiers\n";
	print "\t -term <term_file>: add ACE-NAC from terminus file\n";
	print "\t -solvent <solvent_file>: Add the rest of the system - DMPC, SOL, Na\n";
	print "\t -resname -from_res <FRO> -to_res <TOR>\n";
	print "\t -het2atom\n";
	exit;
}

