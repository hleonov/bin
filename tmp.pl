#!/usr/bin/perl -w
# This script should eventually run all the analysis scripts possible.

$param_file = shift || die "Usage: perl $0 <param-file>\n";

undef(%param);
$vmdt = '/Applications/VMD\ 1.8.6.app/Contents/MacOS/startup.command -dispdev text';

&read_param($param_file);

&rmsd if ($param{'rmsd'} == 1);
&water if ($param{'water'} == 1);
&connectivity if ($param{'connect'} == 1);
&res_dist if ($param{'resdist'} == 1);
&dihedrals if ($param{'dihedrals'} == 1);
&hx_tilt if ($param{'tilt'} == 1);
&contacts if ($param{'contacts'} == 1);
&conductance if ($param{'conductance'} == 1);
&ss_retention if ($param{'ss_ret'} == 1);
&gen_hole_input if ($param{'hole'} == 1);
&group_density if ($param{'g_density'} == 1);
&aa_dist if($param{'aa_dist'} ==1);
&native_contacts if ($param{'native'} == 1);
&cluster if ($param{'cluster'} == 1);

sub read_param {
	open(PAR, $_[0]) || die "Cannot open $_[0]\n$!";
	while (<PAR>) {
		next if (m/^#/ || m/^\s*$/);	#skip comments and empty lines
		s/\s//g;
		@parts = split(/=/);
		next if ($#parts<1); #skip lines without =
		$param{$parts[0]} = $parts[1];
		print "param{$parts[0]} = $parts[1]\n";
	}
	close(PAR);
	my $n = length($param{'trj'});
	$param{'trjtype'} = substr($param{'trj'},$n-3,3);
}

sub native_contacts {
	my $com1 = "editconf -f $param{'em'} -o $param{'em'}.gro";
	my $com2 = "$vmdt -e ~/vmd_scripts/native_contacts.tcl -args $param{'em'}.gro native_ca.dat $param{'trjtype'} $param{'contactskip'} $param{'trj'} $param{'cont-rad'} $param{'hx_dist'} $param{'orig_pdb'} protein and type CA";
	my $com3 = "$vmdt -e ~/vmd_scripts/native_contacts.tcl -args $param{'em'}.gro native_cs.dat $param{'trjtype'} $param{'contactskip'} $param{'trj'} $param{'cont-rad'} $param{'hx_dist'} $param{'orig_pdb'} protein and noh and not backbone";
	system($com1) == 0 || die "native contacts proc: Error in turning from-em.pdb file to .gro\n";
	unless (-e "native_ca.dat") {
		system($com2) == 0 || die "native contacts proc: Error in compute contacts for CA\n";
	}
	system($com3) == 0 || die "native contacts proc: Error in compute contacts for SC\n";

}
sub aa_dist {
	my $com1 = "$vmdt -e ~/vmd_scripts/z_aa.tcl -args $param{'pr'} $param{'trjtype'} $param{'conskip'} $param{'trj'} $param{'ch_len'}";
	my $com2 = "~/Research/m2/run_hist.csh";
	system($com1) == 0 || die "aa_dist proc: VMD Error\n";
	system($com2) == 0 || die "aa_dist proc: histogram Error\n";

}

sub cluster {
	my $com = "printf \"1\n1\n\" | g_cluster -f $param{'trj'} -s $param{'mdtpr'} -n $param{'index'} -b $param{'cl_b_time'} -method gromos -o -clid -minstruct 15 -dist -sz -cl -cutoff $param{'cl_cutoff'}";
	system($com) == 0 || die "Cluster proc: g_cluster Error\n";
}

sub group_density {
	print "======= Executing make index =========\n";
	my $com = "printf \"a P | a NTM\na OW\nq\n\" | make_ndx -f $param{'mdtpr'} -n $param{'index'} -o ../system2.ndx";
	my $com2 = "printf \"1\nP_NTM\nOW\n\" | g_density -f $param{'trj'} -n ../system2.ndx -s $param{'mdtpr'} -o density.xvg -ng 3";
	system($com) == 0 || die "Group Density proc: make-index Error\n";
	system($com2) == 0 || die "Group Density proc: g_density Error\n";
}

sub dihedrals {
	my $i=1;
	while (defined $param{"dhdlist$i"}) {
		my $listfile = $param{"dhdlist$i"};
		my $name = $param{"dhdname$i"};
		my $com = "$vmdt -e ~/vmd_scripts/compute_dihedrals.tcl -args $param{'pr'} $listfile $name\_dihedrals.txt $param{'trjtype'} $param{'dhdskip'} $param{'trj'}";
		system($com) == 0 || die "Residue dihedrals Error: Problem with computing dihedrals for $name\n";
		my $com2 = "perl ~/bin/good_smooth.pl $name\_dihedrals.txt $listfile $name\_dihedrals.sm$param{'smooth_num'}.dat $param{'smooth_num'}";
		system($com2) || print "Residue dihedrals Error: Problem with smooth for $name\n";
		$i++;
	}
}

sub res_dist {
	my $i=1;
	while (defined $param{"dist$i"}) {
			$param{"dist$i"} =~ s/\;/ /g;
#			print $param{"dist$i"}."\n";
			my $name = $param{"dname$i"};
			$com_d = "$vmdt -e ~/vmd_scripts/residue_dist.tcl -args $param{'pr'} $param{'trjtype'} $param{'distskip'} $param{'trj'} $name\_dist.txt $name\_min_dist.txt ".$param{"dist$i"};
#			print $com_d."\n";
			system($com_d) == 0 || die "Residue dist Error: Problem with computing residue distances for $name\n";
			my $com2 = "perl ~/bin/good_smooth.pl $name\_min_dist.txt $name\_min_dist.sm$param{'smooth_num'}.dat $param{'smooth_num'}";
			system($com2) || print "Residue dist Error: Problem with smooth for $name\n";
			$i++;
	}
}

sub connectivity {
	my $com1 = "$vmdt -e ~/vmd_scripts/water_connectivity.tcl -args $param{'pr'} $param{'trj'} connect_$param{'connect_rad'} $param{'res'} $param{'conskip'} $param{'connect_rad'} $param{'trjtype'}";
	system($com1) == 0 || die "Water proc: Connectivity Error\n";	
}

sub water {
	my $com2 = "$vmdt -e ~/vmd_scripts/analyze_waters_fragment.tcl -args $param{'pr'} $param{'trj'} $param{'trjtype'} results $param{'slices'} $param{'range'} $param{'cyl'} $param{'frames'}";
 	my $com3 = "$vmdt -e ~/vmd_scripts/water.tcl -args $param{'slices'} $param{'pr'} $param{'trj'} $param{'conskip'} $param{'trjtype'} $param{'cyl'}";
	system($com2) == 0 || die "Water proc: Excel histogram Error\n";	
	system($com3) == 0 || die "Water proc: Cumulative histogram Error\n";	
}

sub gen_hole_input {
	my $com1 = "printf \"1\n\" | trjconv -f $param{'trj'} -o trj_dt100.pdb -dt 100 -s $param{'mdtpr'}";
	system($com1)  == 0 || die "hole Error: cannot exec trjconv\n";	
	
}
sub hx_tilt {
	my $com = "$vmdt -e ~/vmd_scripts/helix_tilt_trj.tcl -args $param{'pr'} $param{'trjtype'} $param{'tltskip'} $param{'trj'} tilt $param{'helices'}";
	print "command: $com";
	system($com) == 0 || die "Tilt proc: Error in computing tilt\n";
	for my $j (0 .. 3) {	
		my $com2 = "perl ~/bin/good_smooth.pl avg_tilt_sh7_h$j.dat avg_tilt_sh7_h$j.sm$param{'smooth_num'}.dat $param{'smooth_num'}";
		system($com2) || print "Helix tilt Error: Problem with smooth\n";
	
	}
}

sub contacts {
	my $com = "$vmdt -e ~/vmd_scripts/contacts.tcl -args $param{'pr'} contacts.dat $param{'trjtype'} $param{'contactskip'} $param{'trj'} $param{'cont-rad'} $param{'hx_dist'}";
	my $com2 = "perl ~/bin/good_smooth.pl contacts.dat contacts.sm$param{'smooth_num'}.dat $param{'smooth_num'}";
	print "command: $com";
	system($com) == 0 || die "Contacts proc: Error in computing contacts\n";
	system($com2) == 0 || print "Contacts proc: Error with smoothing\n";
}

sub conductance {
	$param{"filter"} =~ s/\;/ /g;
	my $com = "$vmdt -e ~/vmd_scripts/conductance.tcl -args $param{'pr'} cond.log cond.txt $param{'trj'} $param{'trjtype'} $param{'filter'}";
	print "command: $com\n";
	system($com) == 0 || die "Conductance proc: Error in computing conductance\n";
}

sub ss_retention {
	my $com = "$vmdt -e ~/vmd_scripts/ss_ret_ugly2.tcl -args $param{'pr'} hx_ret.txt $param{'trjtype'} $param{'ss-skip'} $param{'trj'} $param{'ss-type'}";
	my $com2 = "perl ~/bin/good_smooth.pl hx_ret.txt hxret.sm$param{'smooth_num'}.dat $param{'smooth_num'}";
	print "command: $com\n";
	system($com) == 0 || die "SS Retention proc: Error in computing ss retention\n";
	system($com2) == 0 || print "SS Retention proc: Error with smoothing\n";

}

sub rmsd {
	#rmsd from after-EM structure
	my $com1 = "printf \"3\n3\n\" | g_rms -s $param{'prtpr'} -f $param{'trj'} -o rmsd_md_from_em.xvg";
    #rmsd from after-PR structure
    my $com2 = "printf \"3\n3\n\" | g_rms -s $param{'mdtpr'} -f $param{'trj'} -o rmsd_md.xvg";
	
	system($com1) == 0 || die "RMSD Error: cannot exec RMSD from EM structure\n";
	system($com2) == 0 || die "RMSD Error: cannot exec RMSD from PR structure\n";
	
	#rmsd of seperate helices - first define helices 
	@hx = split (/\;|:/, $param{'helices'});
	$hx_def = "";
	for (my $i=0; $i<=$#hx; $i++) {
		$hx_def .= "r $hx[$i]\n";			#including sidechains
		$hx_def_ca .= "r $hx[$i] & 3\n";	#define only CA of each helix
		$groups[$i] = "r_$hx[$i]";
		$groups_ca[$i] = "r_$hx[$i]\_\&_C-alpha";
		$hx_out[$i] = "rmsd_h$i.xvg";		#define output file
		$hx_out_ca[$i] = "rmsd_ca_h$i.xvg";	
	}
	print "my helices: $hx_def\n";
	my $ndx_com = "printf \"$hx_def $hx_def_ca q\n\" | make_ndx -f $param{'em'} -n $param{'index'} -o system2.ndx";
	system($ndx_com) == 0 || die "RMSD Error: cannot exec make_ndx\n";
	
	for ($i=0; $i<=$#hx; $i++) {
		$h_com[$i] = "printf \"$groups[$i]\n$groups[$i]\n\" | g_rms -n system2.ndx -s $param{'mdtpr'} -f $param{'trj'} -o $hx_out[$i]";	
		$h_com_ca[$i] = "printf \"$groups_ca[$i]\n$groups_ca[$i]\n\" | g_rms -n system2.ndx -s $param{'mdtpr'} -f $param{'trj'} -o $hx_out_ca[$i]";	
		system($h_com[$i]) == 0 || die "RMSD Error: cannot exec RMSD for helix $i\n";
		system($h_com_ca[$i]) == 0 || die "RMSD Error: cannot exec RMSD for ca-helix $i\n";
	}
}
