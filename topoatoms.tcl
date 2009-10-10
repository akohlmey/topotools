#!/usr/bin/tclsh
# This file is part of TopoTools, a VMD package to simplify 
# manipulating bonds other topology related properties.
#
# Copyright (c) 2009 by Axel Kohlmeyer <akohlmey@gmail.com>
#

# Return info about atoms
# we list and count only bonds that are entirely within the selection.
proc ::TopoTools::atominfo {infotype sel {flag none}} {

    set atomtypes [lsort -ascii -unique [$sel get type]]

    switch $infotype {
        numatoms      { return [$sel num] }
        numatomtypes  { return [llength $atomtypes] }
        atomtypenames { return $atomtypes }
        default       { return "bug? shoot the programmer!"}
    }
}

