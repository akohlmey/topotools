#!/usr/bin/tclsh
# This file is part of TopoTools, a VMD package to simplify 
# manipulating bonds other topology related properties.
#
# Copyright (c) 2009 by Axel Kohlmeyer <akohlmey@cmm.chem.upenn.edu>
#

# Return info about bonds.
# we list and count only bonds that are entirely within the selection.
proc ::TopoTools::bondinfo {infotype sel {flag none}} {

    set numbonds 0
    set bidxlist {}
    array set bondtypes {}

    set aidxlist [$sel list]
    set bondlist [$sel getbonds]
    set btyplist [$sel getbondtypes]
    set bordlist [$sel getbondorders]
    
    foreach a $aidxlist bl $bondlist tl $btyplist ol $bordlist {
        foreach b $bl t $tl o $ol {
            if {($a < $b) && ([lsearch -sorted -integer $aidxlist $b] != -1)} {
                incr numbonds
                switch $flag {
                    type   {lappend bidxlist [list $a $b $t]}
                    order  {lappend bidxlist [list $a $b $o]}
                    both   {lappend bidxlist [list $a $b $t $o]}
                    lammps {lappend bidxlist [list $numbonds $a $b $t]}
                    none   {lappend bidxlist [list $a $b]}
                }
            }
            set bondtypes($t) 1
        }
    }

    switch $infotype {
        numbonds      { return $numbonds }
        numbondtypes  { return [array size bondtypes] }
        bondtypenames { return [array names bondtypes] }
        getbondlist   { return $bidxlist }
        default       { return "bug? shoot the programmer!"}
    }
}

# delete all contained bonds of the selection.
proc ::TopoTools::clearbonds {sel} {
    
    # special optimization for "all" selection.
    if {[string equal "all" [$sel text]]} {
        set nulllist {}
        for {set i 0} {$i < [$sel num]} {incr i} {
            lappend nullist {}
        }
        $sel setbonds $nullist
        return
    }

    set mol [$sel molid]
    foreach b [bondinfo getbondlist $sel none] {
        delbond $mol [lindex $b 0] [lindex $b 1]
    }
}

# reset bonds to data in bondlist
proc ::TopoTools::setbondlist {sel flag bondlist} {

    clearbonds $sel
    set nbnd [llength $bondlist]
    if {$nbnd == 0} { return }
    # set defaults
    set n 0
    set t unknown
    set o 1
    set mol [$sel molid]
    set a -1
    set b -1
    set i 0
    set fract  [expr {100.0/$nbnd}]
    set deltat 2000
    set newt   $deltat 

    # XXX: add sanity check on data format
    foreach bond $bondlist {
        incr i
        set time [clock clicks -milliseconds]
        if {$time > $newt} {
            set percent [format "%3.1f" [expr {$i*$fract}]]
            vmdcon -info "setbondlist: $percent% done."
            display update ui
            set newt [expr {$time + $deltat}]
        }
        switch $flag {
            type   {lassign $bond a b t  }
            order  {lassign $bond a b o  }
            both   {lassign $bond a b t o}
            lammps {lassign $bond n a b t}
            none   {lassign $bond a b    }
        }
        addbond $mol $a $b $t $o
    }
    return
}

# guess bonds type names from atom types.
proc ::TopoTools::retypebonds {sel} {

    set bondlist  [bondinfo getbondlist $sel none]
    set atomtypes [$sel get type]
    set atomindex [$sel list]
    set newbonds {}

    foreach bond $bondlist {
        set idx [lsearch -sorted -integer $atomindex [lindex $bond 0]]
        set a [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex [lindex $bond 1]]
        set b [lindex $atomtypes $idx]
        set type [join [list $a $b] "-"]
        lappend newbonds [list [lindex $bond 0] [lindex $bond 1] $type]
    }
    setbondlist $sel type $newbonds
}


# define a new bond or change and existing.
proc ::TopoTools::addbond {mol id1 id2 type order} {
    if {[catch {atomselect $mol "index $id1 $id2"} sel]} {
        vmdcon -error "topology addbond: Invalid atom indices: $sel"
        return
    }

    # make sure we have consistent indexing
    lassign [$sel list] id1 id2

    set bonds [$sel getbonds]
    set bords [$sel getbondorders]
    set btype [$sel getbondtypes]

    set b1 [lindex $bonds 0]
    set b2 [lindex $bonds 1]
    set bo1 [lindex $bords 0]
    set bo2 [lindex $bords 1]
    set bt1 [lindex $btype 0]
    set bt2 [lindex $btype 1]

    # handle the first atom...
    set pos [lsearch -exact -integer $b1 $id2]
    if { $pos < 0} {
        lappend b1 $id2
        lappend bo1 $order
        lappend bt1 $type
    } else {
        set bo1 [lreplace $bo1 $pos $pos $order]
        set bt1 [lreplace $bt1 $pos $pos $type]
    }

    # ...and the second one.
    set pos [lsearch -exact -integer $b2 $id1]
    if { $pos < 0} {
        lappend b2 $id1
        lappend bo2 $order
        lappend bt2 $type
    } else {
        set bo2 [lreplace $bo2 $pos $pos $order]
        set bt2 [lreplace $bt2 $pos $pos $type]
    }

    # and write the modified data back.
    $sel setbonds [list $b1 $b2]
    if {![string equal $order 1.0]} {
        $sel setbondorders [list $bo1 $bo2]
    }
    if {![string equal $type unknown]} {
        $sel setbondtypes [list $bt1 $bt2]
    }
    $sel delete
}

# delete a bond.
proc ::TopoTools::delbond {mol id1 id2 {type unknown} {order 1.0}} {
    if {[catch {atomselect $mol "index $id1 $id2"} sel]} {
        vmdcon -error "topology delbond: Invalid atom indices: $sel"
        return
    }

    # make sure we have consistent indexing
    lassign [$sel list] id1 id2

    set bonds [$sel getbonds]

    set b1 [lindex $bonds 0]
    set b2 [lindex $bonds 1]

    # handle the first atom...
    set pos [lsearch -exact -integer $b1 $id2]
    if { $pos < 0} {
        ; # bond is not completely within selection. ignore
    } else {
        set b1 [lreplace $b1 $pos $pos]
    }

    # ...and the second one.
    set pos [lsearch -exact -integer $b2 $id1]
    if { $pos < 0} {
        ; # bond is not completely within selection. ignore...
    } else {
        set b2 [lreplace $b2 $pos $pos]
    }

    # and write the modified data back.
    $sel setbonds [list $b1 $b2]
    $sel delete
}
