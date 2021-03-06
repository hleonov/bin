set infile [lindex $argv 0]
set low [lindex $argv 1]
set high [lindex $argv 2]
set out [lindex $argv 3]

mol load gro $infile
#set ca_sel [atomselect top "chain A and protein and type CA"]
#set ca_sel [atomselect top "resid 21 and type NE1 "]
set ca_sel [atomselect top "protein and not backbone and z < $high and z > $low "]

set fid [open $out w]
foreach res [$ca_sel get {name resid resname chain z}] {
	puts $fid "$res"
	set id2 [expr [lindex $res 1] + 27]
	set id3 [expr $id2 + 27]
	set id4 [expr $id3 + 27]
	set type1 [lindex $res 0]
	#set parallel [atomselect top "type CA and resid $id2 $id3 $id4"]
	set parallel [atomselect top "type $type1 and resid $id2 $id3 $id4"]
	#foreach pres [$parallel get {name resid resname chain z}] {
	#	puts $fid "$pres"
	#}
	
}
puts $fid ""
close $fid
exit
