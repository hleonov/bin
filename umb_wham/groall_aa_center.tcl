set resi_begin 1
set resi_end 27
for {set resi $resi_begin} {$resi <= $resi_end} {incr resi} {
	set z_center_arr($resi) 0
}
set fid [open center.dat w]
set half 1
for {set coord -14.5} {$coord <= 14.5} {}  {
	mol delete top
	puts "$coord"
set coord_num $coord
	set filename "center_$coord_num\_pull_k100.gro"
	mol load gro $filename
  	for {set resi $resi_begin} {$resi <= $resi_end} {incr resi} {
	set id2 [expr [lindex $resi 0] + $resi_end]
	set id3 [expr $id2 + $resi_end]
	set id4 [expr $id3 + $resi_end]
	set all_amino_sel [atomselect top "resid $resi $id2 $id3 $id4 and not backbone" ]
	set z_center_arr($resi) [expr $z_center_arr($resi) + [lindex [measure center $all_amino_sel] 2]]
  
  }
  set coord [expr $coord + 0.5]
  	 
}

for {set resi $resi_begin} {$resi <= $resi_end} {incr resi} {
	set aa [atomselect top "type CA and resid $resi" ]
	set name [$aa get resname]
	set z_center [expr $z_center_arr($resi) / 59]
	puts $fid "$resi $name $z_center"
}
exit

#puts "$resi $id2 $id3 $id4\n"
	#foreach pres [$all_amino_sel get {name resid resname chain z}] {
	#	puts "$pres"
	#}
