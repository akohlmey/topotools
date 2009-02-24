#!/usr/bin/tclsh
# TopoTools, a VMD package to simplify manipulating bonds 
# other topology related properties in VMD.
#
# Copyright (c) 2009 by Axel Kohlmeyer <akohlmey@cmm.chem.upenn.edu>
#
package provide topotools 1.0

###################################################
# main frontend command.
###################################################

namespace eval ::TopoTools:: {
    # nothing to do here yet.
}

proc ::TopoTools::usage {} {
    vmdcon -info "usage: topo <command> \[args...\] <flags>"
    vmdcon -info ""
    vmdcon -info "common flags:"
    vmdcon -info "  -molid     <num>|top    molecule id (default: 'top')"
    vmdcon -info "  -sel       <selection>  atom selection function or text"
    vmdcon -info "                          (default: 'all')"
    vmdcon -info "  -bondtype  <typename>   bond type name (default: unknown)"
    vmdcon -info "  -bondorder <bondorder>  bond order parameter (default: 1)"
    vmdcon -info ""
    vmdcon -info "commands:"
    vmdcon -info "  help                    prints this message"
    vmdcon -info ""
    vmdcon -info "  numbonds                returns the number of unique bonds."
    vmdcon -info "  numbondtypes            returns the number of bond types."
    vmdcon -info "  bondtypenames           returns the list of bond types names."
    vmdcon -info "  clearbonds              deletes all bonds. "
    vmdcon -info "  retypebonds             resets all bond types. "
    vmdcon -info ""
    vmdcon -info "  addbond <id1> <id2>     (re-)defines a single bond."
    vmdcon -info "  delbond <id1> <id2>     deletes a single bond, if it exists."
    vmdcon -info ""
    vmdcon -info "  getbondlist \[type|order|both|none\]"
    vmdcon -info "     returns a list of unique bonds, optionally"
    vmdcon -info "     including bond order and bond type."
    vmdcon -info "  setbondlist <list> \[type|order|both|none\]" 
    vmdcon -info "     resets all bonds from a list in the same"
    vmdcon -info "     format as returned by 'topo getbondlist'."
    vmdcon -info "     bond order or -type are reset to defaults if not given."
    vmdcon -info ""
    vmdcon -info "  numangles               returns the number of unique angles"
    vmdcon -info "  numangletypes           returns the number of angle types"
    vmdcon -info "  angletypenames          returns the list of bond type names"
    vmdcon -info "  clearangles             deletes all angles. "
    vmdcon -info "  retypeangles            resets all angle types. "
    vmdcon -info ""
    vmdcon -info "  addangle <id1> <id2> <id3> \[<type>\] (re-defines) a single angle."
    vmdcon -info "  delangle <id1> <id2> <id3>  (re-defines) a single angle."
    vmdcon -info ""
    vmdcon -info ""
    vmdcon -info "  getanglelist  returns the list of angle definitions"
    vmdcon -info "                in the form {type <id1> <id2> <id3>}"
    vmdcon -info "  setanglelist <list>"
    vmdcon -info "                resets angle definitions from a list in the same"
    vmdcon -info "                format as retured by 'topo getanglelist'"
    vmdcon -info ""
    vmdcon -info "NOTE: for angle, dihedral, and improper lists, the"
    vmdcon -info "      type field currently has to be always present."
    vmdcon -info ""
    return
}

# the main frontend command.
# this takes care of all sanity checks on arguments and
# then dispatches the subcommands to the corresponding
# subroutines. 
proc topo { args } {

    set molid -1
    set seltxt all
    set selmol -1
    set bondtype unknown
    set bondorder 1.0

    set cmd ""

    # process arguments
    set newargs {}
    for {set i 0} {$i < [llength $args]} {incr i} {
        set arg [lindex $args $i]

        if {[string match -?* $arg]} {
            
            set val [lindex $args [expr $i+1]]
        
            switch -- $arg {
                -molid { 
                    if {[catch {molinfo $val get name} res]} {
                        vmdcon -error "Invalid -molid argument '$val': $res"
                        return
                    }
                    set molid $val
                    if {[string equal $molid "top"]} {
                        set molid [molinfo top]
                    }
                    incr i
                }

                -sel { 
                    if {[info commands $val] != ""} {
                        if {[catch {$val text} res]} {
                            vmdcon -error "Invalid -sel argument '$val': $res"
                            return
                        }
                        set selmol [$val molid]
                    } else {
                        set res $val
                    }
                    set seltxt $res
                    incr i
                }

                -bondtype { 
                    if {[string length $val] < 1} {
                        vmdcon -error "Invalid -bondtype argument '$val'"
                        return
                    }
                    set bondtype $val
                    incr i
                }

                -bondorder { 
                    if {[string length $val] < 1} {
                        vmdcon -error "Invalid -bondorder argument '$val'"
                        return
                    }
                    set bondorder $val
                    incr i
                }

                -- break

                default {
                    vmdcon -info "default: $arg"
                }
            }
        } else {
            lappend newargs $arg
        }
    }

    if {$molid < 0} { 
        set molid $selmol
    }
    if {$molid < 0} { 
        set molid [molinfo top]
    }
    if {($selmol >= 0) && ($selmol != $molid)} {
        vmdcon -error "Molid from selection '$selmol' does not match -molid argument '$molid'"
        return
    }

    if {[catch {atomselect $molid $seltxt} sel]} {
        vmdcon -error "Problem with atom selection: $sel"
        return
    }

    # branch out to the various subcommands
    set retval ""
    if {[llength $newargs] > 0} {
        set cmd [lindex $newargs 0]
        set newargs [lrange $newargs 1 end]
    } else {
        set newargs {}
        set cmd help
    }

    switch -- $cmd {
        getbondlist   -
        bondtypenames -
        numbondtypes  -
        numbonds {
            if {[llength $newargs] < 1} {set newargs none}
            set retval [::TopoTools::bondinfo $cmd $sel $newargs]
        }

        setbondlist  {
            set flag none
            if {[llength $newargs] > 1} {
                set flag [lindex $newargs 0]
                set newargs [lrange $newargs 1 end]
            }
            if {[llength $newargs] < 1} {set newargs none}
            set retval [::TopoTools::setbondlist $sel $flag [lindex $newargs 0]]
        }

        retypebonds {
            set retval [::TopoTools::retypebonds $sel] 
        }

        clearbonds {
            set retval [::TopoTools::clearbonds $sel] 
        }

        addbond {
            if {[llength $newargs] < 2} {
                vmdcon -error "Not enough arguments for 'topo addbond'"
                ::TopoTools::usage
                return
            }
            set retval [::TopoTools::addbond $molid \
                            [lindex $newargs 0] \
                            [lindex $newargs 1] \
                            $bondtype $bondorder]
        }

        delbond {
            if {[llength $newargs] < 2} {
                vmdcon -error "Not enough arguments for 'topo addbond'"
                ::TopoTools::usage
                return
            }
            set retval [::TopoTools::delbond $molid \
                            [lindex $newargs 0] \
                            [lindex $newargs 1] \
                            $bondtype $bondorder]
        }

        getanglelist   -
        angletypenames -
        numangletypes  -
        numangles {
            if {[llength $newargs] < 1} {set newargs none}
            set retval [::TopoTools::angleinfo $cmd $sel $newargs]
        }

        setanglelist  {
            set retval [::TopoTools::setanglelist $sel [lindex $newargs 0]]
        }

        retypeangles {
            set retval [::TopoTools::retypeangles $sel] 
        }

        clearangles {
            set retval [::TopoTools::clearangles $sel] 
        }

        addangle {
            set atype unknown
            if {[llength $newargs] < 3} {
                vmdcon -error "Not enough arguments for 'topo addangle'"
                ::TopoTools::usage
                return
            }
            if {[llength $newargs] > 3} {
                set atype [lindex $newargs 3]
            }
            set retval [::TopoTools::addangle $molid \
                            [lindex $newargs 0] \
                            [lindex $newargs 1] \
                            [lindex $newargs 2] \
                            $atype]
        }

        delangle {
            set atype unknown
            if {[llength $newargs] < 3} {
                vmdcon -error "Not enough arguments for 'topo delangle'"
                ::TopoTools::usage
                return
            }
            set retval [::TopoTools::delangle $molid \
                            [lindex $newargs 0] \
                            [lindex $newargs 1] \
                            [lindex $newargs 2] ]
        }

        help -
        default {
            ::TopoTools::usage
        }
    }
    $sel delete
    return $retval
}

source [file join $env(TOPOTOOLDIR) topobonds.tcl]
source [file join $env(TOPOTOOLDIR) topoangles.tcl]

