########################################################
## VMD Tcl script to convert psf data into setup_lammps 
## compatible topology file
## To use: load only the psf file into VMD and then run this script
##
## Blake Wilson 
## blake.wilson@utdallas.edu
##
##################################################

set outname "out_test.top"

set outfile [open $outname w]


#get the index arrays




set sel [atomselect top all]
# total number of atoms
set N [$sel num]
$sel delete
######## output atoms
puts "getting atoms.."
for {set i 0} {$i < $N} {incr i} {
	set A [atomselect top "index $i"]
	set index [expr $i + 1]
	set aname [$A get name]
	set aresname [$A get resname]
	set atype [$A get type]
	set amass [$A get mass]
	set acharge [$A get charge]
	puts $outfile "atom        $index    $aresname    $aname $atype      $amass $acharge    $aresname"    
    $A delete 
}

puts $outfile " " 
puts $outfile " "
# bond list - oddly there is no molinfo for the bond list,
# so use the topo tool one
set blist [topo getbondlist]
######### output bonds
puts "getting bonds.."
set nb [llength $blist]
puts "$nb bonds"
for {set i 0} {$i < $nb} {incr i} {
	set B [lindex $blist $i]
	puts $outfile "bond    [expr [lindex $B 0] + 1]    [expr [lindex $B 1] + 1]"

}
puts $outfile " "
puts $outfile " "
unset blist
######### output angles
puts "getting angles.."
# angle list
set alist [molinfo top get angles]
set na [llength [lindex $alist 0]]
#check for empty angle list
set nfl [llength [lindex [lindex $alist 0] 0] ]
if {$nfl==0} {
	set na 0
}
puts "$na angles"
for {set i 0} {$i < $na} {incr i} {
	set A [lindex [lindex $alist 0] $i]
	puts $outfile "angle    [expr [lindex $A 1] + 1]    [expr [lindex $A 2] + 1]    [expr [lindex $A 3] + 1]"


}
unset alist
puts $outfile " "
puts $outfile " "
####### output dihedrals
puts "getting dihedrals.."
# dihedral list
set dlist [molinfo top get dihedrals]
set nd [llength [lindex $dlist 0] ]
#check for empty dihedral list
set nfl [llength [lindex [lindex $dlist 0] 0] ]
if {$nfl==0} {
	set nd 0
}
puts "$nd dihedrals"
for {set i 0} {$i < $nd} {incr i} {
	set D [lindex [lindex $dlist 0] $i]
	set a1 [expr [lindex $D 1] + 1]
	set a2 [expr [lindex $D 2] + 1]
	set a3 [expr [lindex $D 3] + 1]
	set a4 [expr [lindex $D 4] + 1]
	puts $outfile "dihedral      $a1    $a2    $a3    $a4"

}
unset dlist
puts $outfile " "
puts $outfile " "

#######output impropers
puts "getting impropers.."
# improper list
set ilist [molinfo top get impropers]
set ni [llength [lindex $ilist 0]]
#check for empty improper list
set nfl [llength [lindex [lindex $ilist 0] 0] ]
if {$nfl==0} {
	set ni 0
}
puts "$ni impropers"
for {set i 0} {$i < $ni} {incr i} {
	set I [lindex [lindex $ilist 0] $i]
	set a1 [expr [lindex $I 1] + 1]
	set a2 [expr [lindex $I 2] + 1]
	set a3 [expr [lindex $I 3] + 1]
	set a4 [expr [lindex $I 4] + 1]
	puts $outfile "improper      $a1    $a2    $a3    $a4"

}
unset ilist
puts "check output: $outname"
close $outfile
