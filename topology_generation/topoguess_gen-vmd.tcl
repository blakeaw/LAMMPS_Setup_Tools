########################################################
## VMD Tcl script to generate a setup_lammps 
## compatible topology file - this script uses the 
## topotools package to guess all the bonds, anlges, etc. 
## and outputs them in a .top format suitable for setup_lammps
## To use: load the coordinate file into VMD and then run this script
##
## Blake Wilson 
## blake.wilson@utdallas.edu
##
##################################################

package require topotools

set outname "out_topoguess.top"

set outfile [open $outname w]

#use topotools package to guess the correct bonds, angles, etc.
topo guessbonds
topo guessangles
topo guessdihedrals
topo guessimpropers

#get the index arrays
set blist [topo getbondlist]
set alist [topo getanglelist]
set dlist [topo getdihedrallist]
set ilist [topo getimproperlist]

set sel [atomselect top all]
set N [$sel num]

#output atoms
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
}

puts $outfile " " 
puts $outfile " "

#output bonds
puts "getting bonds.."
set nb [llength $blist]
for {set i 0} {$i < $nb} {incr i} {
	set B [lindex $blist $i]
	puts $outfile "bond    [expr [lindex $B 0] + 1]    [expr [lindex $B 1] + 1]"

}
puts $outfile " "
puts $outfile " "
#output angles
puts "getting angles.."
set na [llength $alist]
for {set i 0} {$i < $na} {incr i} {
	set A [lindex $alist $i]
	puts $outfile "angle    [expr [lindex $A 1] + 1]    [expr [lindex $A 2] + 1]    [expr [lindex $A 3] + 1]"


}
puts $outfile " "
puts $outfile " "
#output dihedrals
puts "getting dihedrals.."
set nd [llength $dlist]
for {set i 0} {$i < $nd} {incr i} {
	set D [lindex $dlist $i]
	set a1 [expr [lindex $D 1] + 1]
	set a2 [expr [lindex $D 2] + 1]
	set a3 [expr [lindex $D 3] + 1]
	set a4 [expr [lindex $D 4] + 1]
	puts $outfile "dihedral      $a1    $a2    $a3    $a4"

}
puts $outfile " "
puts $outfile " "
#output impropers
puts "getting impropers.."
set ni [llength $ilist]
for {set i 0} {$i < $ni} {incr i} {
	set I [lindex $ilist $i]
	set a1 [expr [lindex $I 1] + 1]
	set a2 [expr [lindex $I 2] + 1]
	set a3 [expr [lindex $I 3] + 1]
	set a4 [expr [lindex $I 4] + 1]
	puts $outfile "improper      $a1    $a2    $a3    $a4"

}

puts "check output: $outname"
close $outfile
