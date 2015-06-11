#! /usr/sbin/perl

$file = $ARGV[0];
unless($ARGV[1]){
  $dotden = 5
}
else{
  $dotden = $ARGV[1];
}

# removing the .inp if it is there

$file =~ s/\.inp//g;
$file_inp = $file.".inp";
$file_log = $file.".log";
$file_sph = $file.".sph";
$file_sos = $file.".sos";
$file_vmd = $file.".vmd";
$file_qpt = $file.".qpt";
@files = ($file_log, $file_sph, $file_sos, $file_vmd, $file_qpt);
foreach $i (@files){
  if (-e $i){
    unlink($i);
  }
}
# running hole
system("/Users/hleonov/Programs/hole2/exe/hole < $file_inp > $file_log");
# making the sos (dot surface from the sph file)
system("/Users/hleonov/Programs/hole2/exe/sph_process -sos -dotden $dotden -color 2.3 4 $file_sph $file_sos");
# triangulates dot surface for solid rendered surfaces (for VMD)
system("/Users/hleonov/Programs/hole2/exe/sos_triangle -v < $file_sos > $file_vmd");



