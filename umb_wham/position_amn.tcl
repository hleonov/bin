#This script moves amantadine to different starting positions in the channel. 
#starting from center of channel and moving outwards (ACE-n-terminal side)
#by dz Angstrom steps.
#The center of mass of AMN is being repositioned, not NAK. 
# run : vmd -dispdev text -e ./position_amn.tcl -args initial_system.pdb <dz> 1/0


#USER NOTES: 1. which residues are saved if at all
#			 2. dz in command line is in absolute value
#			 3. adjust the limit (14) and the equality, for minus or plus
# 			 4. for the translation of amantadine only in Z - use z_trans=1

mol load pdb [lindex $argv 0]
set dz [lindex $argv 1]
set z_tran [lindex $argv 2]

# First position - center of channel
set prot [atomselect top "backbone"]
set amn [atomselect top "resname AMA"]
set amn_n [atomselect top "type NAK"]
set v1 [measure center $prot]
#set v2 [lindex [$amn_n get {x y z}] 0]
set v2 [measure center $amn weight [$amn get mass]]
set vec_dist [vecsub $v1 $v2]; #vector from amn-->protein
$amn moveby $vec_dist

#check 
#measure center $prot
#$amn_n get {x y z}

#save
set all [atomselect top "all"] 
#or 
set num 0
#set all [atomselect top "all and not resid 199 170"]
 
$all writepdb center_$num.0.pdb

#set low_d [expr $dz/2]
#set high_d [expr $dz+$low_d]

proc move {dz} {
	global z_trans amn num
	set v2 [measure center $amn weight [$amn get mass]]
	set low [expr [lindex $v2 2]-0.5]
	set high [expr [lindex $v2 2]-3.5]
	#if $dz is negative then $low and $high replace eachtother
	if {$low<$high} {
		set sel [atomselect top "(protein or resname ACE NAC) and (z>=$low) and (z<=$high)"]
	} else {
		set sel [atomselect top "(protein or resname ACE NAC) and (z>=$high) and (z<=$low)"]
	}
	set v1 [measure center $sel]
	set vec_dist [vecsub $v1 $v2]
	#lreplace <list> <start_i> <end_i> value (replace avg z by our dz)
	set vec_dist [lreplace $vec_dist 2 2 $dz] 	
	
	if {$z_trans == 1} {
		$amn moveby [list 0 0 $dz]
	} else {
		$amn moveby $vec_dist
	}
	set all [atomselect top "all"] 
	#set all [atomselect top "all and not resid 189"] 
	set num [expr $num+$dz]
	$all writepdb center_$num.pdb
}

#progress by steps of dz on z-axis and adjust to pore center + save
while {$num>=-14} {
	[move [expr -1*$dz]]
}
while {$num<=14} {
	move $dz;
}
exit
