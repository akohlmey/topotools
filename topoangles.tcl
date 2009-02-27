#!/usr/bin/tclsh
# This file is part of TopoTools, a VMD package to simplify 
# manipulating bonds other topology related properties.
#
# Copyright (c) 2009 by Axel Kohlmeyer <akohlmey@cmm.chem.upenn.edu>
#

# return info about angles
# we list and count only angles that are entirely within the selection.
proc ::TopoTools::angleinfo {infotype sel {flag none}} {

    set numangles 0
    array set angletypes {}
    set atidxlist [$sel list]
    set anglelist {}

    foreach angle [join [molinfo [$sel molid] get angles]] {
        lassign $angle t a b c 

        if {([lsearch -sorted -integer $atidxlist $a] >= 0)          \
                && ([lsearch -sorted -integer $atidxlist $b] >= 0)   \
                && ([lsearch -sorted -integer $atidxlist $c] >= 0) } {
            set angletypes($t) 1
            incr numangles
            lappend anglelist $angle
        }
    }
    switch $infotype {

        numangles      { return $numangles }
        numangletypes  { return [array size angletypes] }
        angletypenames { return [array names angletypes] }
        getanglelist   { return $anglelist }
        default        { return "bug! shoot the programmer?"}
    }
}

# delete all contained angles of the selection.
proc ::TopoTools::clearangles {sel} {
    set mol [$sel molid]
    set atidxlist [$sel list]
    set anglelist {}

    foreach angle [join [molinfo $mol get angles]] {
        lassign $angle t a b c 

        if {([lsearch -sorted -integer $atidxlist $a] < 0)          \
                || ([lsearch -sorted -integer $atidxlist $b] < 0)   \
                || ([lsearch -sorted -integer $atidxlist $c] < 0) } {
            lappend anglelist $angle
        }
    }
    molinfo $mol set angles [list $anglelist]
}

# reset angles to data in anglelist
proc ::TopoTools::setanglelist {sel anglelist} {

    set mol [$sel molid]
    set atidxlist [$sel list]
    set newanglelist {}

    # set defaults
    set t unknown; set a -1; set b -1; set c -1

    # preserve all angles definitions that are not contained in $sel
    foreach angle [angleinfo getanglelist $sel] {
        lassign $angle t a b c 

        if {([lsearch -sorted -integer $atidxlist $a] < 0)          \
                || ([lsearch -sorted -integer $atidxlist $b] < 0)   \
                || ([lsearch -sorted -integer $atidxlist $c] < 0) } {
            lappend newanglelist $angle
        }
    }

    # append new ones, but only those contained in $sel
    foreach angle $anglelist {
        lassign $angle t a b c 

        if {([lsearch -sorted -integer $atidxlist $a] >= 0)          \
                && ([lsearch -sorted -integer $atidxlist $b] >= 0)   \
                && ([lsearch -sorted -integer $atidxlist $c] >= 0) } {
            lappend newanglelist $angle
        }
    }

    molinfo $mol set angles [list $newanglelist]
}

# reset angles to data in anglelist
proc ::TopoTools::retypeangles {sel} {

    set mol [$sel molid]
    set anglelist [angleinfo getanglelist $sel]
    set atomtypes [$sel get type]
    set atomindex [$sel list]
    set newangles {}
    
    foreach angle $anglelist {
        lassign $angle type i1 i2 i3

        set idx [lsearch -sorted -integer $atomindex $i1]
        set a [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex $i2]
        set b [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex $i3]
        set c [lindex $atomtypes $idx]

        if { [string compare $a $c] > 0 } { set t $a; set a $c; set c $t }
        set type [join [list $a $b $c] "-"]

        lappend newangles [list $type $i1 $i2 $i3]
    }
    setanglelist $sel $newangles
}

# define a new angle or change an existing one.
proc ::TopoTools::addangle {mol id1 id2 id3 {type unknown}} {
    if {[catch {atomselect $mol "index $id1 $id2 $id3"} sel]} {
        vmdcon -error "topology addangle: Invalid atom indices: $sel"
        return
    }

    # canonicalize indices
    if {$id1 > $id3} {set t $id1 ; set id1 $id3 ; set id3 $t } 

    set angles [join [molinfo $mol get angles]]
    lappend angles [list $type $id1 $id2 $id3]
    molinfo $mol set angles [list $angles]
}

# delete an angle.
proc ::TopoTools::delangle {mol id1 id2 id3 {type unknown}} {
    if {[catch {atomselect $mol "index $id1 $id2 $id3"} sel]} {
        vmdcon -error "topology delangle: Invalid atom indices: $sel"
        return
    }

    # canonicalize indices
    if {$id1 > $id3} {set t $id1 ; set id1 $id3 ; set id3 $t } 

    set newangles {}
    foreach angle [join [molinfo $mol get angles]] {
        lassign $angle t a b c
        if { ($a != $id1) || ($b != $id2) || ($c != $id3) } {
            lappend newangles $angle
        }
    }
    molinfo $mol set angles [list $newangles]
}
