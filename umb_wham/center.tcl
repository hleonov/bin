# output the COM (center of mass) of the amantadine in each file (rc)

for {set i -14.5} {$i<=14.5} {set i [expr $i+0.5]} {
	mol load pdb "center_$i\.pdb"
	set amn [atomselect top "resname AMA"]
	set cent [measure center $amn]
	set cx [expr [lindex $cent 0]/10]
	set cy [expr [lindex $cent 1]/10]
	set cz [expr [lindex $cent 2]/10]
	puts "$i $cx $cy $cz"
	mol delete top
}
