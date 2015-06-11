mol load gro center_0.0_pull_k100.gro
animate read xtc center_0.0_pull_k100.xtc skip 10
set num [molinfo top get numframes]
set resi_begin 1
set resi_end 27

for {set resi $resi_begin} {$resi <= $resi_end} {incr resi} {
	set z_center_arr($resi) 0
}
for {set fr 0} {$fr <= $num} {incr fr} {
  for {set resi $resi_begin} {$resi <= $resi_end} {incr resi} {
	#set aa [atomselect top "type CA and resid $resi" frame $fr]
	#set name [$aa get resname]
	set id2 [expr [lindex $resi 0] + $resi_end]
	set id3 [expr $id2 + $resi_end]
	set id4 [expr $id3 + $resi_end]
	set all_amino_sel [atomselect top "resid $resi $id2 $id3 $id4 and not backbone" frame $fr]
	set z_center_arr($resi) [expr $z_center_arr($resi) + [lindex [measure center $all_amino_sel] 2]]
  }
}

for {set resi $resi_begin} {$resi <= $resi_end} {incr resi} {
	set aa [atomselect top "type CA and resid $resi" frame 0]
	set name [$aa get resname]
	set z_center [expr $z_center_arr($resi) / ($num+1)]
	puts "$resi $name $z_center"
}
#puts "$resi $id2 $id3 $id4\n"
#foreach pres [$all_amino_sel get {name resid resname chain z}] {
	#	puts "$pres"
	#}
