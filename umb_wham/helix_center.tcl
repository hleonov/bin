mol load gro center_0.0_pull_k100.gro
animate read xtc center_0.0_pull_k100.xtc skip 10 waitfor all
set num [molinfo top get numframes]
set resi_begin 1
set resi_end 27
set fid [open helix_center_bb.dat w]

for {set h 0} {$h <= 3} {incr h} {
	for {set resi $resi_begin} {$resi <= $resi_end} {incr resi} {
		set z_center_arr($h,$resi) 0
	}
}
for {set fr 0} {$fr <= $num} {incr fr} {
  for {set resi $resi_begin} {$resi <= $resi_end} {incr resi} {
  	for {set h 0} {$h <= 3} {incr h} {
		set id [expr [lindex $resi 0] + $resi_end*$h]
	#	puts "resid : $id\n";
  		set amino_sel [atomselect top "resid $id and not sidechain" frame $fr]
		set z_center_arr($h,$resi) [expr $z_center_arr($h,$resi) + [lindex [measure center $amino_sel] 2]]	
	}	
		
	#set aa [atomselect top "type CA and resid $resi" frame $fr]
	#set name [$aa get resname]
  }
}

for {set h 0} {$h <= 3} {incr h} {
	for {set resi $resi_begin} {$resi <= $resi_end} {incr resi} {
		set aa [atomselect top "type CA and resid $resi" frame 0]
		set name [$aa get resname]
		set z_center [expr $z_center_arr($h,$resi) / ($num+1)]
		puts $fid "$h $resi $name $z_center"
	}
}
#puts "$resi $id2 $id3 $id4\n"
#foreach pres [$all_amino_sel get {name resid resname chain z}] {
	#	puts "$pres"
	#}
exit;
