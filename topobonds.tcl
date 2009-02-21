#!/usr/bin/tclsh

# return info about bonds
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
            if {$a < $b} {
                incr numbonds
                switch $flag {
                    type  lappend bidxlist [list $a $b $t]
                    order lappend bidxlist [list $a $b $o]
                    both  lappend bidxlist [list $a $b $t $o]
                    none  lappend bidxlist [list $a $b]
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
        default       { return "bug! shoot the programmer?"}
    }
}

# define a new bond or change and existing.
proc ::TopoTools::addbond {mol idx1 idx2 order type} {
    if {[catch {atomselect $mol "index $idx1 $idx2"} sel]} {
        vmdcon -error "topology addbond: Invalid atom indices: $sel"
        return
    }

    # make sure we have consistent indexing
    lassign [$sel list] idx1 idx2

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
    set pos [lsearch -exact -integer $b1 $idx2]
    if { $pos < 0} {
        lappend b1 $idx2
        lappend bo1 $order
        lappend bt1 $type
    } else {
        set bo1 [lreplace $bo1 $pos $pos $order]
        set bt1 [lreplace $bt1 $pos $pos $type]
    }

    # ...and the second one.
    set pos [lsearch -exact -integer $b2 $idx1]
    if { $pos < 0} {
        lappend b2 $idx1
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
