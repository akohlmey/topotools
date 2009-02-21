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
    # we may want to write a simple GUI at some point.
}

proc ::TopoTools::usage {} {
    vmdcon -info "usage: topology <command> \[args...\]"
    vmdcon -info ""
    vmdcon -info "commands:"
    vmdcon -info ""
    vmdcon -info "  help          prints this message"
    vmdcon -info ""
    vmdcon -info "  numbonds      returns the number of unique bonds"
    vmdcon -info "  numbondtypes  returns the number of bond types"
    vmdcon -info "  getbondlist \[type|order|both\]"
    vmdcon -info "        returns a list of unique bonds, optionally"
    vmdcon -info "        including either bond order or bond type."
    vmdcon -info "  setbondlist <list> \[type|order|both\] resets bonds from a list"
    vmdcon -info "  addbond <idx1> <idx2>  defines a new bond"
    vmdcon -info "  delbond <idx1> <idx2>  deletes a bond, if it exists"
    vmdcon -info ""
    vmdcon -info "common arguments:"
    vmdcon -info "  -molid  <num>|top        molecule id"
    vmdcon -info "  -sel    <selection>      atom selection or text"
    vmdcon -info "  -bondtype  <typename>    bond type name"
    vmdcon -info "  -bondorder <bondorder>   bond order parameter"
    vmdcon -info ""
    return
}

# the main frontend command.
# this takes care of all sanity checks on arguments and
# then dispatches the subcommands to the corresponding
# subroutines. 
proc topology { args } {

    set molid top
    set seltxt all
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
                    incr i
                }

                -sel { 
                    if {[info commands $val] != ""} {
                        if {[catch {$val text} res]} {
                            vmdcon -error "Invalid -sel argument '$val': $res"
                            return
                        }
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

    if {[string equal $molid "top"]} {
        set molid [molinfo top]
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
    vmdcon -info "command is: $cmd ,  args: $newargs / $bondorder $bondtype"
    switch -- $cmd {
        getbondlist   -
        bondtypenames -
        numbondtypes  -
        numbonds {
            set retval [::TopoTools::bondinfo $cmd $sel $newargs]
        }

        addbond {
            set retval [::TopoTools::addbond $molid \
                            [lindex $newargs 0] \
                            [lindex $newargs 1] \
                            $bondorder $bondtype]
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

