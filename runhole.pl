#! /usr/bin/perl -w
# Run HOLE on input PDBs. 
# There are two modes of action:
# 1. single PDB - is moved to (0,0,0) first, then HOLE.
# 2. PDB with multiple models - first extract single PDBs, then run HOLE.

# Usage example: 
#	Single mode:	perl ~/bin/runhole.pl S-WT_from_md system_from_MD.pdb 1


$file = shift || die 	"Usage: perl $0 <out_base> <coord> <num> [sph]\n";
$coord = shift || die 	"Usage: perl $0 <out_base> <coord> <num> [sph]\n";
$num = shift || die 	"Usage: perl $0 <out_base> <coord> <num> [sph]\n";
$vmd = shift;

$file_inp = $file.".inp";
$file_log = $file.".log";
$file_sph = $file.".sph";
$file_sos = $file.".sos";
$file_vmd = $file.".vmd";
$file_qpt = $file.".qpt";

#system("cp /home/hadasleo/Research/m2/md/mutants_take2/Singapor_WT/backup.inp .");

#single PDB
if ($num == 1) {
	$coord =~ s/.pdb//;
	#move to 0 0 0
	system("vmd -dispdev text -e ~/vmd_scripts/move_bottom_left_to_0_0_0.tcl -args $coord.pdb $coord"."_0_0_0.pdb");
	$coord = $coord."_0_0_0.pdb";
	write_inp($coord);
#PDB of multiple models
} elsif ($num>1) {
	$models = 	`grep MODEL $coord | wc -l` - 1;
	system("~/bin/pdb_exploder $coord $models");
	$coord =~ s/.pdb//;
	write_inp("$coord\_*");
	#system("sed 's/^coord.*/coord $coord\*/' backup.inp > $file_inp");
}
print "Running HOLE...";
system("/home/hadasleo/Programs/hole2/exe/hole < $file_inp > $file_log");
print "Done\n";

if ($num>1) {
	`rm $coord\_*`;
}

if ($vmd) {
	print "Running sph & sos triangulation for VMD output...\n";
	system("/home/hadasleo/Programs/hole2/exe/sph_process -sos -color 2.3 4 $file_sph $file_sos");
	system("/home/hadasleo/Programs/hole2/exe/sos_triangle -v < $file_sos > $file_vmd");
	print "Done\n"
}

$res_out = "$file.rad.dat";
system("perl ~/bin/read_hole_results.pl $file_log $res_out matrix");

#system("/home/hadasleo/Programs/hole2/exe/sos_triangle -v < $file_sos > $file_vmd");

#sub plot {
#	system("cp /home/hadasleo/Research/m2/md/mutants_take2/Singapor_WT/pore/plot_hole_radius backup_hole.plot");
#	system("sed 's/^set title.*/set title \"$coord\"/' backup_hole.plot > tmp");
#	system("sed 's/^splot .* using/splot \"$file_gnu\" using/' hole_plot");
#}

sub write_inp {
	open(INP, ">$file_inp") || die "Cannot write inp file $file_inp\n $! \n";
	print INP	"! first use standard spherical probe\n";
	print INP	"!MUST cards\n";
	print INP	"coord  $_[0]\n";
	print INP	"radius /home/hadasleo/Programs/hole2/rad/xplor.rad	! Use xplor vdw radii\n";
	print INP	"cvect 0.0 0.0 1.0		! channel runs in Z direction\n";
	print INP	"!cpoint 0.02 0.03 -0.3\n";
	print INP	"! now optional cards\n";
	print INP	"CONNOLLY                        ! move away from spherical probe - use the new connolly option\n";
	print INP	"ignore SOL\n";
	print INP	"dotden 10\n";
	print INP	"sample 0.25\n";
	print INP	"pltout $file_qpt                       ! output qpt format quanta plot file\n";
	print INP	"SPHPDB $file_sph\n";
	print INP	"endrad 7                                        ! to avoid having enormous ends\n";
	print INP	"shorto 0\n";
	close(INP);

}
