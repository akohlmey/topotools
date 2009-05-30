#!/usr/bin/tclsh
# TopoTools, a VMD package to simplify manipulating bonds 
# other topology related properties in VMD.

# some little helper functions

# angle definition list comparison function
proc ::TopoTools::compareangles {a b} {
    if {[lindex $a 1] < [lindex $b 1]} {
        return -1
    } elseif {[lindex $a 1] > [lindex $b 1]} {
        return 1
    } else {
        if {[lindex $a 2] < [lindex $b 2]} {
            return -1
        } elseif {[lindex $a 2] > [lindex $b 2]} {
            return 1
        } else {
            if {[lindex $a 3] < [lindex $b 3]} {
                return -1
            } elseif {[lindex $a 3] > [lindex $b 3]} {
                return 1
            } else {
                if {[llength $a] == 4} {
                    return 0
                } else {
                    if {[lindex $a 3] < [lindex $b 3]} {
                        return -1
                    } elseif  {[lindex $a 3] > [lindex $b 3]} {
                        return 1
                    } else {
                        return 0
                    }
                }
            }
        }
    }
}

# sort angle/dihedral/improper list and remove duplicates
proc ::TopoTools::sortsomething {what sel} {
    switch $what {
        angle     {
            setanglelist $sel [lsort -unique -command compareangles \
                                     [angleinfo getanglelist $sel]]
        }
        dihedral  {
            setdihedrallist $sel [lsort -unique -command comparedihedrals \
                                     [dihedralinfo getdihedrallist $sel]]
        }
        improper  {
            setimproperlist $sel [lsort -unique -command compareimpropers \
                                     [improperinfo getimproperlist $sel]]
        }
    }
}

# emulate the behavior of loading a molecule through
# the regular "mol new" command.
proc ::TopoTools::adddefaultrep {mol} {
    mol color [mol default color]
    mol rep [mol default style]
    mol selection [mol default selection]
    mol material [mol default material]
    mol addrep $mol
    display resetview
}
