#!/usr/bin/tclsh
# This file is part of TopoTools, a VMD package to simplify 
# manipulating bonds other topology related properties.
#
# Copyright (c) 2009 by Axel Kohlmeyer <akohlmey@cmm.chem.upenn.edu>
#

# return info about dihedrals
# we list and count only dihedrals that are entirely within the selection.
proc ::TopoTools::dihedralinfo {infotype sel {flag none}} {

    set numdihedrals 0
    array set dihedraltypes {}
    set atidxlist [$sel list]
    set dihedrallist {}

    foreach dihedral [join [molinfo [$sel molid] get dihedrals]] {
        lassign $dihedral t a b c d

        if {([lsearch -sorted -integer $atidxlist $a] >= 0)          \
                && ([lsearch -sorted -integer $atidxlist $b] >= 0)   \
                && ([lsearch -sorted -integer $atidxlist $c] >= 0)   \
                && ([lsearch -sorted -integer $atidxlist $d] >= 0) } {
            set dihedraltypes($t) 1
            incr numdihedrals
            lappend dihedrallist $dihedral
        }
    }
    switch $infotype {

        numdihedrals      { return $numdihedrals }
        numdihedraltypes  { return [array size dihedraltypes] }
        dihedraltypenames { return [array names dihedraltypes] }
        getdihedrallist   { return $dihedrallist }
        default        { return "bug! shoot the programmer?"}
    }
}

# delete all contained dihedrals of the selection.
proc ::TopoTools::cleardihedrals {sel} {
    set mol [$sel molid]
    set atidxlist [$sel list]
    set dihedrallist {}

    foreach dihedral [join [molinfo $mol get dihedrals]] {
        lassign $dihedral t a b c d

        if {([lsearch -sorted -integer $atidxlist $a] < 0)          \
                || ([lsearch -sorted -integer $atidxlist $b] < 0)   \
                || ([lsearch -sorted -integer $atidxlist $c] < 0)   \
                || ([lsearch -sorted -integer $atidxlist $d] < 0) } {
            lappend dihedrallist $dihedral
        }
    }
    molinfo $mol set dihedrals [list $dihedrallist]
}

# reset dihedrals to data in dihedrallist
proc ::TopoTools::setdihedrallist {sel dihedrallist} {

    set mol [$sel molid]
    set atidxlist [$sel list]
    set newdihedrallist {}

    # set defaults
    set t unknown; set a -1; set b -1; set c -1; set d -1

    # preserve all dihedrals definitions that are not contained in $sel
    foreach dihedral [dihedralinfo getdihedrallist $sel] {
        lassign $dihedral t a b c d

        if {([lsearch -sorted -integer $atidxlist $a] < 0)          \
                || ([lsearch -sorted -integer $atidxlist $b] < 0)   \
                || ([lsearch -sorted -integer $atidxlist $c] < 0)   \
                || ([lsearch -sorted -integer $atidxlist $d] < 0) } {
            lappend newdihedrallist $dihedral
        }
    }

    # append new ones, but only those contained in $sel
    foreach dihedral $dihedrallist {
        lassign $dihedral t a b c d

        if {([lsearch -sorted -integer $atidxlist $a] >= 0)          \
                && ([lsearch -sorted -integer $atidxlist $b] >= 0)   \
                && ([lsearch -sorted -integer $atidxlist $c] >= 0)   \
                && ([lsearch -sorted -integer $atidxlist $d] >= 0) } {
            lappend newdihedrallist $dihedral
        }
    }

    molinfo $mol set dihedrals [list $newdihedrallist]
}

# reset dihedrals to data in dihedrallist
proc ::TopoTools::retypedihedrals {sel} {

    set mol [$sel molid]
    set dihedrallist [dihedralinfo getdihedrallist $sel]
    set atomtypes [$sel get type]
    set atomindex [$sel list]
    set newdihedrals {}
    
    foreach dihedral $dihedrallist {
        lassign $dihedral type i1 i2 i3 i4

        set idx [lsearch -sorted -integer $atomindex $i1]
        set a [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex $i2]
        set b [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex $i3]
        set c [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex $i4]
        set d [lindex $atomtypes $idx]

        if { [string compare $b $c] > 0 } { 
            set t $a; set a $d; set d $t 
            set t $b; set b $c; set c $t 
        }
        set type [join [list $a $b $c $d] "-"]

        lappend newdihedrals [list $type $i1 $i2 $i3 $i4]
    }
    setdihedrallist $sel $newdihedrals
}


# define a new dihedral or change an existing one.
proc ::TopoTools::adddihedral {mol id1 id2 id3 id4 {type unknown}} {
    if {[catch {atomselect $mol "index $id1 $id2 $id3 $id4"} sel]} {
        vmdcon -error "topology adddihedral: Invalid atom indices: $sel"
        return
    }

    # canonicalize indices
    if {$id2 > $id3} {
        set t $id2 ; set id2 $id3 ; set id3 $t 
        set t $id1 ; set id1 $id4 ; set id4 $t 
    }

    set dihedrals [join [molinfo $mol get dihedrals]]
    lappend dihedrals [list $type $id1 $id2 $id3 $id4]
    molinfo $mol set dihedrals [list $dihedrals]
}

# delete a dihedral.
proc ::TopoTools::deldihedral {mol id1 id2 id3 id4 {type unknown}} {
    if {[catch {atomselect $mol "index $id1 $id2 $id3 $id4"} sel]} {
        vmdcon -error "topology deldihedral: Invalid atom indices: $sel"
        return
    }

    # canonicalize indices
    if {$id2 > $id3} {
        set t $id2 ; set id2 $id3 ; set id3 $t 
        set t $id1 ; set id1 $id4 ; set id4 $t 
    }

    set newdihedrals {}
    foreach dihedral [join [molinfo $mol get dihedrals]] {
        lassign $dihedral t a b c d
        if { ($a != $id1) || ($b != $id2) || ($c != $id3) || ($d != $id4) } {
            lappend newdihedrals $dihedral
        }
    }
    molinfo $mol set dihedrals [list $newdihedrals]
}
