#!/usr/bin/tclsh
# TopoTools, a VMD package to simplify manipulating bonds 
# other topology related properties in VMD.
#
# Copyright (c) 2009 by Axel Kohlmeyer <akohlmey@gmail.com>
# $Id: topohelpers.tcl,v 1.6 2009/11/20 19:03:32 akohlmey Exp $

# some little helper functions

# compare two lists element by element.
# return 0 if they are identical, or 1 if not.
proc ::TopoTools::listcmp {a b} {
    if {[llength $a] != [llength $b]} {
        return 1
    }
    foreach aa $a bb $b {
        if {![string equal $aa $bb]} {
            return 1
        }
    }
    return 0
}

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
# the regular "mol new" command. the options $selmod
# argument allows to append an additional modified to
# the selection, e.g. 'user > 0.1' for variable number 
# particle xyz trajectories.
proc ::TopoTools::adddefaultrep {mol {selmod none}} {
    mol color [mol default color]
    mol rep [mol default style]
    if {[string equal $selmod none]} {
        mol selection [mol default selection]
    } else {
        mol selection "([mol default selection]) and $selmod"
    }        
    mol material [mol default material]
    mol addrep $mol
    display resetview
}
