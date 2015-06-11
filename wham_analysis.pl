#!/usr/bin/perl -w
# This script runs WHAM and creates the output graph : 
#	PMF vs. Z axis with amino acid probability density along Z (Angstroms)
# Prerequisities:
# *.npdo files, *.hist files, meta.dat.
# Usage: ../../../wham_analysis.pl meta_500.dat free-aa 500

#setup WHAM parameters
$meta = shift || die "Usage: perl $0 <meta> <out> <force>\n";
$out = shift || die "Usage: perl $0 <meta> <out> <force>\n";
$force = shift || die "Usage: perl $0 <meta> <out> <force>\n";
$min_hist = 1.8;
$max_hist = 5.2;
$num_bins = 20;
$tol = 1;
$temp = 300;
$padding = 0; 

undef(%rc); 

# Execute WHAM
`~/Programs/wham/wham $min_hist $max_hist $num_bins $tol $temp $padding $meta $out`;

# Handle free energy output file
$max = `sort -n +1 $out | tail -1 | awk '{print \$2}'`;
$div = int ($max / 10);
$max_free = $div*10 + 10;	#factor to multiply probabily density so that Y axis ia 0-1 for it.
print "factor: $max_free\n";
open(FREE, $out) || die "Error: Cannot open $out\n $!\n";
<FREE>; 	#skip title line
while (<FREE>) {
	if (m/(\S+)\s+(\S+)\s+/) {
		my $z = $1*10;	#switch to Angstrom
		my $free = $2;	
		next if ($free eq 'inf');
		$rc{$z} = $free;
	}
}
close (FREE);

# Handle XMGRACE output
open(XVG, ">$out.agr") || die "Error: Cannot open $out.xvg\n $!\n";

&format_view;
&color_map;
$all_data = &rw_data;	#read hist+free energy and set
#write settings for each data set
for $i (0 .. $#set_format) {
	print XVG $set_format[$i];
}
print XVG $all_data; 	#write data

close(XVG);

# General xmgrace settings
sub format_view {
	#Get name of variant for title
	my $pwd = `pwd`;
	$pwd =~ m/umb-wham\/(\w+)\_(\w+)\//;
	my $variant = $1." ".$2;	
	my $format_text =  "# Grace project file\n#\n".
					"\@page size 842, 594\n".
					"\@background color 0\n".
					"\@with g0\n".
					"\@    world 0, 0, 60, 20\n".
					"\@    stack world 0, 0, 0, 0\n".
					"\@    znorm 1\n".
					"\@    view 0.150000, 0.150000, 1.150000, 0.850000\n".
					"\@    title \"Umbrella Sampling of Amantadine in M2\"\n".
					"\@    title font 0\n".
					"\@    title size 1.500000\n".
					"\@    title color 1\n".
					"\@    subtitle \"$variant, k=$force \"\n".
					"\@    subtitle font 0\n".
					"\@    subtitle size 1.000000\n".
					"\@    subtitle color 1\n".
					"\@    xaxes scale Normal\n".
					"\@    yaxes scale Normal\n".
					"\@    xaxes invert off\n".
					"\@    yaxes invert off\n".
					"\@    xaxis  label \"Z-axis (A)\"\n".
					"\@    yaxis  label \"Free energy\"\n".
					"\@	   legend on\n".
					"\@    legend loctype view\n".
					"\@    legend 1.2, 0.8\n".
					"\@    legend box linewidth 1.0\n".
					"\@    legend box linestyle 1\n".
					"\@    legend font 0\n".
					"\@    legend char size 0.530000\n".
					"\@    legend color 1\n".
					"\@    legend length 2\n".
					"\@    legend vgap 1\n".
					"\@    legend hgap 0\n";
	print XVG $format_text;				
}

sub rw_data {
	my $set = 0;
	$all_data = "";
	# set array so we can read amino acids in resid order
	foreach $hist (`ls -1 *.hist`) {
		if ($hist =~ m/(\d+)/){
			my $rid = $1;
			$hist_arr[$rid] = $hist;
		}
	}
	# read amino acid hist files and set their format data. 
	for my $i (1 .. $#hist_arr) {
		set_format_data($hist_arr[$i], $set);	
		$all_data .= "\@target G0.S$set\n\@type xy\n";
		open(HIS, $hist_arr[$i]) || die "Cannot open $hist_arr[$i]\n $!";
		while (<HIS>) {
			if (m/(\S+)\s+(\S+)/) {
				my $his_z = $1;
				my $his_pr = $2*$max_free;	#multiply so it is 0-1 on Y axis
				$all_data .= "$his_z $his_pr\n";
			}
		}
		close(HIS);
		$all_data .= "\&\n";
		$set++;
	}
	$all_data .= "\@target G0.S$set\n\@type xy\n";
	foreach $z (sort keys(%rc)) {
		$all_data .= $z." ".$rc{$z}."\n";
	}
	$all_data .= "\&\n";
	pmf_format_data($set);
	return $all_data;
	
}

# Set the formatting data for each amino acid probability density
sub set_format_data {
	my $name = $_[0];
	my $set = $_[1];
	$name =~ m/res\.(\w{3})\w?_(\d+)/;
	my $resname = $1;
	my $resid = $2;
	my $color = $color_map{$resname};
	my $linestyle = 1;
	my $linewidth = 1.0;
	if ($resname =~ m/ACE/ || $resname =~ m/NAC/) {
		$linestyle = 4;
		$linewidth = 1.5;
	}
#	$resname =~ tr/([A-Z])([A-Z])$/[a-z][a-z]/;
	$set_format[$set] = "\@    s$set type xy\n".
						"\@    s$set symbol 0\n".
					 	"\@    s$set symbol size 1.000000\n".
					 	"\@    s$set symbol color $color\n".
					 	"\@    s$set symbol pattern 1\n".
					 	"\@    s$set symbol linewidth $linewidth\n".
					 	"\@    s$set line type 1\n".
						"\@    s$set line linewidth $linewidth\n".
					 	"\@    s$set line linestyle $linestyle\n".
					 	"\@    s$set line color $color\n".
					 	"\@    s$set comment \"$resname $resid\"\n".
					 	"\@    s$set legend  \"$resname $resid\"\n";

						 
}

# Set the formatting data for the PMF. 
sub pmf_format_data {
	my $set = $_[0];
	$set_format[$set] = "\@    s$set type xy\n".
						"\@    s$set symbol 1\n".
					 	"\@    s$set symbol size 1.000000\n".
					 	"\@    s$set symbol color 1\n".
					 	"\@    s$set symbol pattern 1\n".
					 	"\@    s$set symbol linewidth 1.0\n".
					 	"\@    s$set line type 1\n".
					 	"\@    s$set line linestyle 1\n".
						"\@    s$set line linewidth 1.5\n".
					 	"\@    s$set line color 1\n";
}
 
sub color_map {
	#Based on the default VMD color map. Can be changed if I define the colors in the file.
	%color_map = (ALA => 1, CYS => 5, ASP => 2, GLU => 2, PHE => 12, GLY => 5,
				  HIS => 4, ILE => 13, LYS => 4, LEU => 6, MET => 5, ASN => 8, 
				  PRO => 9, GLN => 8, ARG => 4, SER => 11, THR => 11, VAL => 10,
				  TRP => 3, TYR => 15, ACE => 4, NAC => 2);

}

