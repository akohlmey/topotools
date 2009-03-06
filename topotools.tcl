#!/usr/bin/tclsh
# TopoTools, a VMD package to simplify manipulating bonds 
# other topology related properties in VMD.
#
# TODO: 
# - topotools.tcl : use namespace variables, to cache metadata for faster execution.
#                   e.g. store lastmol/lastsel and if it the same as before, reuse data 
#                   that is present, or else invalidate cache etc.
#                   need to add a generic API for that, and options to clear and avoid caching.
# - topoatoms.tcl : provide frontend to atom name/types/numbers similar to bonds/angles/...
# - topogmx.tcl   : interface to gromacs (at least for postprocessing)
# - topoamber.tcl : 
#
# Copyright (c) 2009 by Axel Kohlmeyer <akohlmey@cmm.chem.upenn.edu>
#
package provide topotools 1.0

###################################################
# main frontend command.
###################################################

namespace eval ::TopoTools:: {
    variable version 1.0; # for allowing compatibility checks in scripts depending on this 
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


proc ::TopoTools::usage {} {
    vmdcon -info "usage: topo <command> \[args...\] <flags>"
    vmdcon -info ""
    vmdcon -info "common flags:"
    vmdcon -info "  -molid     <num>|top    molecule id (default: 'top')"
    vmdcon -info "  -sel       <selection>  atom selection function or text"
    vmdcon -info "                          (default: 'all')"
    vmdcon -info "flags only applicable to 'bond' commands:"
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
    vmdcon -info "  num(angle|dihedral|improper)s       returns the number of unique (angle|dihedral|improper)s"
    vmdcon -info "  num(angle|dihedral|improper)types   returns the number of (angle|dihedral|improper) types"
    vmdcon -info "  (angle|dihedral|improper)typenames  returns the list of bond type names"
    vmdcon -info "  clear(angle|dihedral|improper)s     deletes all (angle|dihedral|improper)s. "
    vmdcon -info "  sort(angle|dihedral|improper)s      sorts the list of (angle|dihedral|improper)s"
    vmdcon -info "                                      according to atom index and removes duplicates"
    vmdcon -info "  retype(angle|dihedral|improper)s    resets all angle types. "
    vmdcon -info ""
    vmdcon -info "  addangle <id1> <id2> <id3> \[<type>\] (re-defines) a single angle."
    vmdcon -info "  delangle <id1> <id2> <id3>  (re-defines) a single angle."
    vmdcon -info "  add(dihedral|improper) <id1> <id2> <id3> <id4> \[<type>\] (re-defines) a single (dihedral|improper)."
    vmdcon -info "  del(dihedral|improper) <id1> <id2> <id3> <id4> (re-defines) a single (dihedral|improper)."
    vmdcon -info ""
    vmdcon -info ""
    vmdcon -info "  getanglelist  returns the list of angle definitions"
    vmdcon -info "                in the form {type <id1> <id2> <id3>}"
    vmdcon -info "  setanglelist <list>"
    vmdcon -info "                resets angle definitions from a list in the same"
    vmdcon -info "                format as retured by 'topo getanglelist'"
    vmdcon -info "  get(dihedral|improper)list  returns the list of (dihedral|improper) definitions"
    vmdcon -info "                in the form {type <id1> <id2> <id3> <id4>}"
    vmdcon -info "  set(dihedral|improper)list <list>"
    vmdcon -info "                resets (dihedral|improper) definitions from a list in the same"
    vmdcon -info "                format as retured by 'topo get(dihedral|improper)list'"
    vmdcon -info "NOTE: for angle, dihedral, and improper lists, the"
    vmdcon -info "      type field currently has to be always present."
    vmdcon -info ""
    vmdcon -info ""
    vmdcon -info "  readlammpsdata <filename> \[<atomstyle>\]"
    vmdcon -info "      read atom properties, bond, angle, dihedral and other related data"
    vmdcon -info "      from a LAMMPS data file. 'atomstyle' is the value given to the 'atom_style'"
    vmdcon -info "      parameter. default value is 'full'."
    vmdcon -info "      the molecule this info is being added to must have a matching number of atoms."
    vmdcon -info "      the -sel parameter is currently ignored."
    vmdcon -info ""
    vmdcon -info "  writelammpsdata <filename> \[<atomstyle>\]"
    vmdcon -info "      write atom properties, bond, angle, dihedral and other related data"
    vmdcon -info "      to a LAMMPS data file. 'atomstyle' is the value given to the 'atom_style'"
    vmdcon -info "      parameter. default value is 'full'."
    vmdcon -info "      Only data that is present is written. "
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

    # process generic arguments and remove them
    # from argument list.
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

    set retval ""
    if {[llength $newargs] > 0} {
        set cmd [lindex $newargs 0]
        set newargs [lrange $newargs 1 end]
    } else {
        set newargs {}
        set cmd help
    }

    if { ![string equal $cmd help] } {
        if {($selmol >= 0) && ($selmol != $molid)} {
            vmdcon -error "Molid from selection '$selmol' does not match -molid argument '$molid'"
            return
        }
        if {[catch {atomselect $molid $seltxt} sel]} {
            vmdcon -error "Problem with atom selection: $sel"
            return
        }
    }

    # branch out to the various subcommands
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

        sortangles {
            set retval [::TopoTools::sortsomething angle $sel] 
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

        getdihedrallist   -
        dihedraltypenames -
        numdihedraltypes  -
        numdihedrals {
            if {[llength $newargs] < 1} {set newargs none}
            set retval [::TopoTools::dihedralinfo $cmd $sel $newargs]
        }

        setdihedrallist  {
            set retval [::TopoTools::setdihedrallist $sel [lindex $newargs 0]]
        }

        retypedihedrals {
            set retval [::TopoTools::retypedihedrals $sel] 
        }

        sortdihedrals {
            set retval [::TopoTools::sortsomething dihedral $sel] 
        }

        cleardihedrals {
            set retval [::TopoTools::cleardihedrals $sel] 
        }

        adddihedral {
            set atype unknown
            if {[llength $newargs] < 4} {
                vmdcon -error "Not enough arguments for 'topo adddihedral'"
                ::TopoTools::usage
                return
            }
            if {[llength $newargs] > 4} {
                set atype [lindex $newargs 4]
            }
            set retval [::TopoTools::adddihedral $molid \
                            [lindex $newargs 0] \
                            [lindex $newargs 1] \
                            [lindex $newargs 2] \
                            [lindex $newargs 3] \
                            $atype]
        }

        deldihedral {
            set atype unknown
            if {[llength $newargs] < 4} {
                vmdcon -error "Not enough arguments for 'topo deldihedral'"
                ::TopoTools::usage
                return
            }
            set retval [::TopoTools::deldihedral $molid \
                            [lindex $newargs 0] \
                            [lindex $newargs 1] \
                            [lindex $newargs 2] \
                            [lindex $newargs 3] ]
        }

        getimproperlist   -
        impropertypenames -
        numimpropertypes  -
        numimpropers {
            if {[llength $newargs] < 1} {set newargs none}
            set retval [::TopoTools::improperinfo $cmd $sel $newargs]
        }

        setimproperlist  {
            set retval [::TopoTools::setimproperlist $sel [lindex $newargs 0]]
        }

        retypeimpropers {
            set retval [::TopoTools::retypeimpropers $sel] 
        }

        sortimpropers {
            set retval [::TopoTools::sortsomething improper $sel] 
        }

        clearimpropers {
            set retval [::TopoTools::clearimpropers $sel] 
        }

        addimproper {
            set atype unknown
            if {[llength $newargs] < 4} {
                vmdcon -error "Not enough arguments for 'topo addimproper'"
                ::TopoTools::usage
                return
            }
            if {[llength $newargs] > 4} {
                set atype [lindex $newargs 4]
            }
            set retval [::TopoTools::addimproper $molid \
                            [lindex $newargs 0] \
                            [lindex $newargs 1] \
                            [lindex $newargs 2] \
                            [lindex $newargs 3] \
                            $atype]
        }

        delimproper {
            set atype unknown
            if {[llength $newargs] < 4} {
                vmdcon -error "Not enough arguments for 'topo delimproper'"
                ::TopoTools::usage
                return
            }
            set retval [::TopoTools::delimproper $molid \
                            [lindex $newargs 0] \
                            [lindex $newargs 1] \
                            [lindex $newargs 2] \
                            [lindex $newargs 3] ]
        }

        readlammpsdata {
            set style atom
            if {[llength $newargs] < 1} {
                vmdcon -error "Not enough arguments for 'topo readlammpsdata'"
                ::TopoTools::usage
                return
            }
            set fname [lindex $newargs 0]
            if {[llength $newargs] > 1} {
                set style [lindex $newargs 1]
            }
            set retval [::TopoTools::readlammpsdata $molid $fname $style $sel]
        }

        writelammpsdata {
            set style atom
            if {[llength $newargs] < 1} {
                vmdcon -error "Not enough arguments for 'topo writelammpsdata'"
                ::TopoTools::usage
                return
            }
            set fname [lindex $newargs 0]
            if {[llength $newargs] > 1} {
                set style [lindex $newargs 1]
            }
            set retval [::TopoTools::writelammpsdata $molid $fname $style $sel]
        }

        help -
        default {
            ::TopoTools::usage
        }
    }
    if {[info exists sel]} {
        $sel delete
    }
    return $retval
}

source [file join $env(TOPOTOOLSDIR) topobonds.tcl]
source [file join $env(TOPOTOOLSDIR) topoangles.tcl]
source [file join $env(TOPOTOOLSDIR) topodihedrals.tcl]
source [file join $env(TOPOTOOLSDIR) topoimpropers.tcl]
source [file join $env(TOPOTOOLSDIR) topolammps.tcl]
