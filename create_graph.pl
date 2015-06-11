#!/usr/bin/perl -w
# 13.8.08	Hadas Leonov
# Version 1
# This script should read in a data file (or a list?) containing nxy form data.
# It will also read a file containing the text data for each graph (title, legend)
# Output: agr format graph --> put into xmgrace, turn to ps, then to pdf

# list format should include the filename and which columns from 
# this file should be plotted
# filename	Xcol colN1	colN2 ..

$listfile = shift || die "Usage: perl $0 <list> <gr_params> <out>\n";
$paramfile = shift || die "Usage: perl $0 <list> <gr_params> <out>\n";
$outfile = shift || die "Usage: perl $0 <list> <gr_params> <out>\n";

undef(%j2set);
undef(@data);
undef(%t);

#Read list file. 
#split into an array where each line i contains filename and cols
open(LIST, $listfile);
	$i=0;
	while (<LIST>) {
		@{$list[$i]} = split(/\s+/, $_);
		$i++;
	}
close(LIST);

# Read the graph data
$totalN = read_data();

#read graph settings
read_param();

make_graph();

sub make_graph {
	open(OUT, ">$outfile") || die "Cannot open output file : $outfile\n";
	 print OUT &default_settings2;
	 print  OUT &format_view;
	for (my $s=0; $s<$totalN; $s++) {	
		print OUT &set_format($s);
	}
	for (my $s=0; $s<$totalN; $s++) {
	#	print "print out data for set $s\n";
		print OUT "\@target G0.S$s\n\@type xy\n";
		print OUT $data[$s]."\&\n";
	}
	close(OUT);
}

sub read_data {
	my $set=0;
	#iterating over files
	for $i (0 .. $#list) {
		print "Reading file: $list[$i][0]\n";
		open(IN,$list[$i][0]) || die "Cannot open file $list[$i][0]\n";
		
		#one file would have the same x coordinates for all (requirement!)
		my $xcol = $list[$i][1];
		
		#reset every new set to "" and create a mapping from cols to set#
		for $j (2 .. $#{$list[$i]}) {
			$j2set{$j} = $set;
			$data[$j2set{$j}] = "\n";
			$set++;
		}
		print "Number of sets: ".($#{$list[$i]}-1)."\n";
		while (<IN>) {
			s/^\s+//g;	#remove leading spaces
			next if (m/^\#|\@/);
			@line = split(/\s+/, $_);
			for $j (2 .. $#{$list[$i]}) {
				$data[$j2set{$j}] .= $line[$xcol]."\t".$line[$list[$i][$j]]."\n";
				#print "data: $data[$j2set{$j}]\n";
	#			$set_format[$j2set{$j}] = set_format($j2set{$j});
			}
		}
		close(IN);
	}
	print "Total number of sets: $set\n\n";
	return $set;
}

sub read_param {
	open(PAR, $paramfile) || die "Cannot open $_[0]\n$!";
	print "Reading Graph parameters:\n";
	while (<PAR>) {
		next if (m/^#/ || m/^\s*$/);	#skip comments and empty lines
		chomp;
		#s/\s//g;
		@parts = split(/=/);
		$parts[0] =~ s/\s*$//g;
		next if ($#parts<1); #skip lines without =
		$t{$parts[0]} = $parts[1];
		print "t{$parts[0]} = $parts[1]\n";
	}
	close(PAR);
}

sub set_format {
	my $set = $_[0];
	my $color = $set+1;
	my $leg = "legend".$set;
	#print "Configuring set $set\n";
	my $set_format = "\@    s$set type xy\n".
					"\@    s$set symbol 0\n".
				 	"\@    s$set symbol size 1.000000\n".
				 	"\@    s$set symbol color $color\n".
				 	"\@    s$set symbol pattern 1\n".
				 	"\@    s$set symbol linewidth 1.0\n".
				 	"\@    s$set line type 1\n".
					"\@    s$set line linewidth 1.0\n".
				 	"\@    s$set line linestyle 1\n".
				 	"\@    s$set line color $color\n".
#				 	"\@    s$set comment \"$resname $resid\"\n".
				 	"\@    s$set legend  \"$t{$leg} \"\n";
	
	return $set_format;

}


sub format_view {
		$format_text =  "\@background color 0\n".
					"\@with g0\n".
					"\@    world $t{xmin}, $t{ymin}, $t{xmax}, $t{ymax}\n".
					"\@    stack world 0, 0, 0, 0\n".
					"\@    znorm 1\n".
					"\@    view 0.150000, 0.150000, 1.150000, 0.850000\n".
					"\@    title \"$t{title}\"\n".
					"\@    title font 0\n".
					"\@    title size 1.500000\n".
					"\@    title color 1\n".
					"\@    subtitle \"$t{subtitle}\"\n".
					"\@    subtitle font 0\n".
					"\@    subtitle size 1.000000\n".
					"\@    subtitle color 1\n".
					"\@    xaxes scale Normal\n".
					"\@    yaxes scale Normal\n".
					"\@    xaxes invert off\n".
					"\@    yaxes invert off\n".
					"\@    xaxis  label \"$t{xlabel} \" \n".
					"\@    xaxis  tick on\n".
					"\@    xaxis  tick major $t{xmajortick}\n".
					"\@    xaxis  tick minor ticks 1\n".
					"\@    xaxis  tick default 6\n".
					"\@    xaxis  tick place rounded true\n".
					"\@    xaxis  tick in\n".
					"\@    xaxis  ticklabel on\n".
					"\@    xaxis  ticklabel format general\n".
					"\@    xaxis  ticklabel prec 5\n".
					"\@    xaxis  ticklabel formula \"$t{xformula}\"\n".
					"\@	   yaxis on\n".
					"\@    yaxis  label \"$t{ylabel}\"\n".
					"\@    yaxis  tick on\n".
					"\@    yaxis  tick major $t{ymajortick}\n".
					"\@    yaxis  tick minor ticks 1\n".
					"\@    yaxis  ticklabel on	\n".
					"\@    yaxis  ticklabel formula \"$t{yformula}\"\n".
					"\@	   legend on\n".
					"\@    legend loctype view\n".
					"\@    legend $t{legx}, $t{legy}\n".
					"\@    legend box linewidth 1.0\n".
					"\@    legend box color 1\n".
					"\@    legend box pattern 1\n".
					"\@    legend box fill color 0\n".
					"\@    legend box linestyle 1\n".
					"\@    legend font 0\n".
					"\@    legend char size 1.000000\n".
					"\@    legend color 1\n".
					"\@    legend length 2\n".
					"\@    legend vgap 1\n".
					"\@    legend hgap 1\n";
					
	return $format_text;

}


sub default_settings {
	$default_str = "# Grace project file\n#\n".
					"\@version 50118\n".
					"\@page size 792, 612\n".
					"\@page scroll 5%\n".
					"\@page inout 5%\n".
					"\@link page off\n".
					"\@map font 0 to \"Times-Roman\", \"Times-Roman\"\n".
					"\@map font 1 to \"Times-Italic\", \"Times-Italic\"\n".
					"\@map font 2 to \"Times-Bold\", \"Times-Bold\"\n".
					"\@map font 3 to \"Times-BoldItalic\", \"Times-BoldItalic\"\n".
					"\@map font 4 to \"Helvetica\", \"Helvetica\"\n".
					"\@map font 5 to \"Helvetica-Oblique\", \"Helvetica-Oblique\"\n".
					"\@map font 6 to \"Helvetica-Bold\", \"Helvetica-Bold\"\n".
					"\@map font 7 to \"Helvetica-BoldOblique\", \"Helvetica-BoldOblique\"\n".
					"\@map font 8 to \"Courier\", \"Courier\"\n".
					"\@map font 9 to \"Courier-Oblique\", \"Courier-Oblique\"\n".
					"\@map font 10 to \"Courier-Bold\", \"Courier-Bold\"\n".
					"\@map font 11 to \"Courier-BoldOblique\", \"Courier-BoldOblique\"\n".
					"\@map font 12 to \"Symbol\", \"Symbol\"\n".
					"\@map font 13 to \"ZapfDingbats\", \"ZapfDingbats\"\n".
					"\@map color 0 to (255, 255, 255), \"white\"\n".
					"\@map color 1 to (0, 0, 0), \"black\"\n".
					"\@map color 2 to (255, 0, 0), \"red\"\n".
					"\@map color 3 to (0, 255, 0), \"green\"\n".
					"\@map color 4 to (0, 0, 255), \"blue\"\n".
					"\@map color 5 to (255, 255, 0), \"yellow\"\n".
					"\@map color 6 to (188, 143, 143), \"brown\"\n".
					"\@map color 7 to (220, 220, 220), \"grey\"\n".
					"\@map color 8 to (148, 0, 211), \"violet\"\n".
					"\@map color 9 to (0, 255, 255), \"cyan\"\n".
					"\@map color 10 to (255, 0, 255), \"magenta\"\n".
					"\@map color 11 to (255, 165, 0), \"orange\"\n".
					"\@map color 12 to (114, 33, 188), \"indigo\"\n".
					"\@map color 13 to (103, 7, 72), \"maroon\"\n".
					"\@map color 14 to (64, 224, 208), \"turquoise\"\n".
					"\@map color 15 to (0, 139, 0), \"green4\"\n".
					"\@reference date 0\n".
					"\@date wrap off\n".
					"\@date wrap year 1950\n".
					"\@default linewidth 1.0\n".
					"\@default linestyle 1\n".
					"\@default color 1\n".
					"\@default pattern 1\n".
					"\@default font 0\n".
					"\@default char size 1.000000\n".
					"\@default symbol size 1.000000\n".
					"\@default sformat \"%.8g\"\n".
					"\@background color 0\n".
					"\@page background fill on\n";
					
	return $default_str;
}
		
sub default_settings2 {
	$default_str2 = "# Grace project file\n#\n".
					"\@version 50121\n".
					"\@page size 792, 612\n".
					"\@page scroll 5%\n".
					"\@page inout 5%\n".
					"\@link page off\n".
					"\@map font 0 to \"Times-Roman\", \"Times-Roman\"\n".
					"\@map font 1 to \"Times-Italic\", \"Times-Italic\"\n".
					"\@map font 2 to \"Times-Bold\", \"Times-Bold\"\n".
					"\@map font 3 to \"Times-BoldItalic\", \"Times-BoldItalic\"\n".
					"\@map font 4 to \"Helvetica\", \"Helvetica\"\n".
					"\@map font 5 to \"Helvetica-Oblique\", \"Helvetica-Oblique\"\n".
					"\@map font 6 to \"Helvetica-Bold\", \"Helvetica-Bold\"\n".
					"\@map font 7 to \"Helvetica-BoldOblique\", \"Helvetica-BoldOblique\"\n".
					"\@map font 8 to \"Courier\", \"Courier\"\n".
					"\@map font 9 to \"Courier-Oblique\", \"Courier-Oblique\"\n".
					"\@map font 10 to \"Courier-Bold\", \"Courier-Bold\"\n".
					"\@map font 11 to \"Courier-BoldOblique\", \"Courier-BoldOblique\"\n".
					"\@map font 12 to \"Symbol\", \"Symbol\"\n".
					"\@map font 13 to \"ZapfDingbats\", \"ZapfDingbats\"\n".
					"\@map color 0 to (255, 255, 255), \"white\"\n".
					"\@map color 1 to (0, 0, 0), \"black\"\n".
					"\@map color 2 to (255, 0, 0), \"red\"\n".
					"\@map color 3 to (0, 255, 0), \"green\"\n".
					"\@map color 4 to (0, 0, 255), \"blue\"\n".
					"\@map color 5 to (255, 255, 0), \"yellow\"\n".
					"\@map color 6 to (188, 143, 143), \"brown\"\n".
					"\@map color 7 to (148, 0, 211), \"violet\"\n".
					"\@map color 15 to (220, 220, 220), \"grey\"\n".				
					"\@map color 9 to (0, 255, 255), \"cyan\"\n".
					"\@map color 10 to (255, 0, 255), \"magenta\"\n".
					"\@map color 11 to (255, 165, 0), \"orange\"\n".
					"\@map color 12 to (114, 33, 188), \"indigo\"\n".
					"\@map color 13 to (103, 7, 72), \"maroon\"\n".
					"\@map color 14 to (64, 224, 208), \"turquoise\"\n".
					"\@map color 8 to (0, 139, 0), \"green4\"\n".
					"\@reference date 0\n".
					"\@date wrap off\n".
					"\@date wrap year 1950\n".
					"\@default linewidth 1.0\n".
					"\@default linestyle 1\n".
					"\@default color 1\n".
					"\@default pattern 1\n".
					"\@default font 0\n".
					"\@default char size 1.000000\n".
					"\@default symbol size 1.000000\n".
					"\@default sformat \"%.8g\"\n".
					"\@background color 0\n".
					"\@page background fill on\n";
					
	return $default_str2;
}
		
