#!/usr/bin/tclsh
# This file is part of TopoTools, a VMD package to simplify
# manipulating bonds other topology related properties.
#
# Copyright (c) 2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2022,2023 by Axel Kohlmeyer <akohlmey@gmail.com>
# $Id: topolammps.tcl,v 1.47 2023/04/21 05:41:03 johns Exp $

# high level subroutines for LAMMPS support.
#
# import a LAMMPS data file.
# this behaves almost like a molfile plugin and will create a
# new molecule (but no representation) and return its molecule id.
#
# Arguments:
# filename = name of data file
# style = atomstyle
# flags = more flags. (currently not used)
proc ::TopoTools::readlammpsdata {filename style {flags none}} {
    global M_PI
    if {[catch {open $filename r} fp]} {
        vmdcon -err "readlammpsdata: problem opening data file: $fp\n"
        return -1
    }

    # parse lammps header section.
    array set lammps [readlammpsheader $fp]
    if {$lammps(atoms) < 1} {
        vmdcon -err "readlammpsdata: failed to parse lammps data header. abort."
        return -1
    }

    # create an empty molecule and timestep
    set mol -1
    if {[catch {mol new atoms $lammps(atoms)} mol]} {
        vmdcon -err "readlammpsdata: problem creating empty molecule: $mol"
        return -1
    } else {
        animate dup $mol
    }
    mol rename $mol [file tail $filename]
    set sel [atomselect $mol all]
    set boxdim {}

    if {($lammps(xy) != {}) && ($lammps(xz) != {}) && ($lammps(yz) != {})} {
        set lammps(triclinic) 1
        vmdcon -info "readlammpsdata: detected triclinic cell."
        set a [expr {$lammps(xhi) - $lammps(xlo)}]
        set ly [expr {$lammps(yhi) - $lammps(ylo)}]
        set lz [expr {$lammps(zhi) - $lammps(zlo)}]
        set b [expr {sqrt($ly*$ly + $lammps(xy)*$lammps(xy))}]
        set c [expr {sqrt($lz*$lz + $lammps(xz)*$lammps(xz)
                          + $lammps(yz)*$lammps(yz))}]
        set alpha [expr {($lammps(xy)*$lammps(xz) + $ly*$lammps(yz))/($b*$c)}]
        set beta  [expr {$lammps(xz)/$c}]
        set gamma [expr {$lammps(xy)/$b}]
        set alpha [expr {90.0 - asin($alpha)*180.0/$M_PI}]
        set beta  [expr {90.0 - asin($beta)*180.0/$M_PI}]
        set gamma [expr {90.0 - asin($gamma)*180.0/$M_PI}]

        set boxdim [list $a $b $c $alpha $beta $gamma]
        molinfo $mol set {a b c alpha beta gamma} $boxdim
        lappend boxdim $lammps(triclinic) $a $ly $lz $lammps(xy) $lammps(xz) $lammps(yz)
    } else {
        set $lammps(triclinic) 0
        set a [expr {$lammps(xhi) - $lammps(xlo)}]
        set b [expr {$lammps(yhi) - $lammps(ylo)}]
        set c [expr {$lammps(zhi) - $lammps(zlo)}]
        set boxdim [list $a $b $c 90.0 90.0 90.0]
        molinfo $mol set {a b c alpha beta gamma} $boxdim
        lappend boxdim $lammps(triclinic) $a $b $c 0.0 0.0 0.0
    }
    set atomidmap {}
    set atommasses {}
    set atomlabels {}
    set bondlabels {}
    set anglelabels {}
    set dihedrallabels {}
    set improperlabels {}

    # now loop through the file until we find a known section header.
    # then call a subroutine that parses this section. those subroutines
    # all take the current line number as last argument and return the
    # new value, so we can keep track of the current line and print more
    # useful error messages or warnings.
    set lineno $lammps(lineno)
    while {[gets $fp line] >= 0} {
        incr lineno
        if {[regexp {^\s*Atoms} $line ]} {
            # use atom style indicated by CGCMM header or use style hint comment.
            set stylehint $lammps(style)
            regexp {^\s*Atoms\s+#\s*([a-z]+)} $line x stylehint
            # for requested atom style 'auto' use the hint value instead
            if {[string equal $style auto]} { set style $stylehint }
            if {[string equal $style unknown]} {
                vmdcon -warn "readlammpsdata: automatic atom style detection requested,"
                vmdcon -warn "readlammpsdata: but no atom style hints in data file."
                vmdcon -warn "readlammpsdata: assuming atom style 'full' instead."
                set style {full}
            }
            # check for atom style consistency
            if {![string equal unknown $stylehint] && ![string equal $style $stylehint]} {
                vmdcon -warn "readlammpsdata: requested atom style '$style' is different from"
                vmdcon -warn "readlammpsdata: style hint '$stylehint' encoded into data file."
                vmdcon -warn "readlammpsdata: this may not work. Continuing with '$style'"
            }
            # if atom style is supported
            if { ![regexp {^(atomic|bond|angle|molecular|charge|full|sphere)} $style] } {
                vmdcon -err "readlammpsdata: unsupported atom style '$style'"
                return -1
            }
            set lineno [readlammpsatoms $fp $sel $style $lammps(cgcmm) $boxdim atomlabels $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Atoms section."
                return -1
            }
            # retrieve map of atomids from user field and convert back to integer.
            set atomidmap {}
            foreach id [$sel get user] {
                lappend atomidmap [expr {int($id + 0.5)}]
            }
        } elseif {[regexp {^\s*Velocities} $line ]} {
            set lineno [readlammpsvelocities $fp $sel $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Velocities section."
                return -1
            }
        } elseif {[regexp {^\s*Masses} $line ]} {
            set lineno [readlammpsmasses $fp $mol $lammps(atomtypes) $lammps(cgcmm) atommasses atomlabels $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Masses section."
                return -1
            }
        } elseif {[regexp {^\s*Bonds} $line ]} {
            if {[llength $atomidmap] < 1} {
                vmdcon -err "readlammpsdata: Atoms section must come before Bonds in data file"
                return -1
            }
            set lineno [readlammpsbonds $fp $sel $lammps(bonds) $atomidmap bondlabels $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Bonds section."
                return -1
            }
        } elseif {[regexp {^\s*Angles} $line ]} {
            if {[llength $atomidmap] < 1} {
                vmdcon -err "readlammpsdata: Atoms section must come before Angles in data file"
                return -1
            }
            set lineno [readlammpsangles $fp $sel $lammps(angles) $atomidmap anglelabels $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Angles section."
                return -1
            }
        } elseif {[regexp {^\s*Dihedrals} $line ]} {
            if {[llength $atomidmap] < 1} {
                vmdcon -err "readlammpsdata: Atoms section must come before Dihedrals in data file"
                return -1
            }
            set lineno [readlammpsdihedrals $fp $sel $lammps(dihedrals) $atomidmap dihedrallabels $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Dihedrals section."
                return -1
            }
        } elseif {[regexp {^\s*Impropers} $line ]} {
            if {[llength $atomidmap] < 1} {
                vmdcon -err "readlammpsdata: Atoms section must come before Impropers in data file"
                return -1
            }
            set lineno [readlammpsimpropers $fp $sel $lammps(impropers) $atomidmap improperlabels $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Impropers section."
                return -1
            }
        } elseif {[regexp {^\s*(Atom Type Labels)} $line ]} {
            set lineno [readlammpslabels $fp $mol Atom $lammps(atomtypes) atomlabels $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Atom Type Labels section."
                return -1
            }
        } elseif {[regexp {^\s*(Bond Type Labels)} $line ]} {
            set lineno [readlammpslabels $fp $mol Bond $lammps(bondtypes) bondlabels $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Bond Type Labels section."
                return -1
            }
        } elseif {[regexp {^\s*(Angle Type Labels)} $line ]} {
            set lineno [readlammpslabels $fp $mol Angle $lammps(angletypes) anglelabels $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Angle Type Labels section."
                return -1
            }
        } elseif {[regexp {^\s*(Dihedral Type Labels)} $line ]} {
            set lineno [readlammpslabels $fp $mol Dihedral $lammps(dihedraltypes) dihedrallabels $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Dihedral Type Labels section."
                return -1
            }
        } elseif {[regexp {^\s*(Improper Type Labels)} $line ]} {
            set lineno [readlammpslabels $fp $mol Improper $lammps(impropertypes) improperlabels $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Improper Type Labels section."
                return -1
            }
        } elseif {[regexp {^\s*(Pair Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(atomtypes) $lineno]
        } elseif {[regexp {^\s*(PairIJ Coeffs)} $line ]} {
            set skip [expr {$lammps(atomtypes)*($lammps(atomtypes)+1)/2}]
            set lineno [skiplammpslines $fp $skip $lineno]
        } elseif {[regexp {^\s*(Bond Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(bondtypes) $lineno]
        } elseif {[regexp {^\s*(Angle Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(angletypes) $lineno]
        } elseif {[regexp {^\s*(BondBond Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(angletypes) $lineno]
        } elseif {[regexp {^\s*(BondAngle Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(angletypes) $lineno]
        } elseif {[regexp {^\s*(Dihedral Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(dihedraltypes) $lineno]
        } elseif {[regexp {^\s*(MiddleBondTorsion Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(dihedraltypes) $lineno]
        } elseif {[regexp {^\s*(EndBondTorsion Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(dihedraltypes) $lineno]
        } elseif {[regexp {^\s*(AngleTorsion Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(dihedraltypes) $lineno]
        } elseif {[regexp {^\s*(AngleAngleTorsion Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(dihedraltypes) $lineno]
        } elseif {[regexp {^\s*(BondBond13 Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(dihedraltypes) $lineno]
        } elseif {[regexp {^\s*(Improper Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(impropertypes) $lineno]
        } elseif {[regexp {^\s*(AngleAngle Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(impropertypes) $lineno]
        } elseif { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines silently
        } else {
            vmdcon -err "readlammpsdata: unkown header line: $lineno : $line "
            vmdcon -err "readlammpsdata: cannot continue. "
            return -1
        }
        set lammps(lineno) $lineno
    }
    close $fp

    # apply masses. Atoms section sets a default of 1.0.
    # since the Masses section can appear before the Atoms section
    # we have to set it here after the parsing.
    if {[llength atommasses] > 0} {
        # we have type labels, but masses are indexed by numeric type
        if {([lindex $atommasses 0] == 1) && ([llength atomlabels] > 0)} {
            foreach {t m} $atommasses {d l} $atomlabels {
                set msel [atomselect $mol "type '$l'"]
                $msel set mass $m
                $msel delete
            }
        } else {
            foreach {t m} $atommasses {
                set msel [atomselect $mol "type '$t'"]
                $msel set mass $m
                $msel delete
            }
        }
    }
    mol reanalyze $mol
    variable newaddsrep
    if {$newaddsrep} {
        adddefaultrep $mol
    }
    $sel delete
    return $mol
}

# skip over a given number of non-empty, non-comment lines
proc ::TopoTools::skiplammpslines {fp num lineno} {

    while {[gets $fp line] >= 0} {
        if {$num <= 0} break
        incr lineno
        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
          incr num -1
        }
    }
    return $lineno
}

# read lammps header section from opened file
# and return as an array.
proc ::TopoTools::readlammpsheader {fp} {
    array set lammps {
        atoms 0 atomtypes 0 bonds 0 bondtypes 0 angles 0 angletypes 0
        dihedrals 0 dihedraltypes 0 impropers 0 impropertypes 0 xtrabond 0
        xlo -0.5 xhi 0.5 ylo -0.5 yhi 0.5 zlo -0.5 zhi 0.5 xy {} xz {} yz {}
        lineno 0 cgcmm 0 triclinic 0 style unknown typemass 1 typelabels 0
    }
    set x {}

    vmdcon -info "parsing LAMMPS header."

    # first header line is skipped by LAMMPS. so we put a flag to
    # detect CMM style CG data files with additional information.
    gets $fp line
    if {[string match "*CGCMM*" $line]} {
        set lammps(cgcmm) 1
        vmdcon -info "detected CGCMM style file. will try to parse additional data."
        if {[string match "*atom_style*" $line]} {
            if { [regexp {^.*atom_style\s+(atomic|bond|angle|molecular|charge|full|sphere)\s*.*}  $line x lammps(style) ] } {
                vmdcon -info "Probable atom_style: $lammps(style)"
            }
        }
    }
    if {[string match "*LABELMAP*" $line]} {
        set lammps(typelabels) 1
        vmdcon -info "detected LABELMAP style file. will try to parse additional data."
        if {[string match "*atom_style*" $line]} {
            if { [regexp {^.*atom_style\s+(atomic|bond|angle|molecular|charge|full|sphere)\s*.*}  $line x lammps(style) ] } {
                vmdcon -info "Probable atom_style: $lammps(style)"
            }
        }
    }
    set lineno 1
    set offs [tell $fp]
    set lammps(lineno) $lineno

    #
    while {[gets $fp line] >= 0} {
        incr lineno
        if { [      regexp {^\s*(\d+)\s+atoms}     $line x lammps(atoms) ] } {
        } elseif { [regexp {^\s*(\d+)\s+bonds}     $line x lammps(bonds) ] } {
        } elseif { [regexp {^\s*(\d+)\s+angles}    $line x lammps(angles)] } {
        } elseif { [regexp {^\s*(\d+)\s+dihedrals} $line x lammps(dihedrals)] } {
        } elseif { [regexp {^\s*(\d+)\s+impropers} $line x lammps(impropers)] } {
        } elseif { [regexp {^\s*(\d+)\s+atom types}     $line x lammps(atomtypes)] } {
        } elseif { [regexp {^\s*(\d+)\s+bond types}     $line x lammps(bondtypes)] } {
        } elseif { [regexp {^\s*(\d+)\s+angle types}    $line x lammps(angletypes)] } {
        } elseif { [regexp {^\s*(\d+)\s+dihedral types} $line x lammps(dihedraltypes)] } {
        } elseif { [regexp {^\s*(\d+)\s+improper types} $line x lammps(impropertypes)] } {
        } elseif { [regexp {^\s*(\d+)\s+extra bond per atom} $line x lammps(xtrabond)] } {
        } elseif { [regexp {^\s*([-[:digit:].Ee+]+)\s+([-[:digit:].Ee+]+)\s+xlo xhi} $line \
                        x lammps(xlo) lammps(xhi)] } {
        } elseif { [regexp {^\s*([-[:digit:].Ee+]+)\s+([-[:digit:].Ee+]+)\s+ylo yhi} $line \
                        x lammps(ylo) lammps(yhi)] } {
        } elseif { [regexp {^\s*([-[:digit:].Ee+]+)\s+([-[:digit:].Ee+]+)\s+zlo zhi} $line \
                        x lammps(zlo) lammps(zhi)] } {
        } elseif { [regexp {^\s*([-[:digit:].Ee+]+)\s+([-[:digit:].Ee+]+)\s+([-[:digit:].Ee+]+)\s+xy\s+xz\s+yz} $line x lammps(xy) lammps(xz) lammps(yz)] } {
        } elseif { [regexp {^\s*(\#.*|)$} $line ] } {
        } elseif {[regexp {^\s*(Atoms|Velocities|Masses|Shapes|Dipoles|Bonds|Angles|Dihedrals|Impropers|(Pair|Bond|Angle|Dihedral|Improper) Coeffs|(Atom|Bond|Angle|Dihedral|Improper) Type Labels)} $line ]} {
            seek $fp $offs start
            break
        } else {
            vmdcon -warn "readlammpsheader: skipping unkown header line: $lineno : $line "
        }
        set offs [tell $fp]
        set lammps(lineno) $lineno
    }

    return [array get lammps]
}

# parse atom section
proc ::TopoTools::readlammpsatoms {fp sel style cgcmm boxdata typemap lineno} {
    upvar $typemap atomlabels
    global M_PI
    set numatoms [$sel num]
    set atomdata {}
    set boxx 0.0
    set boxy 0.0
    set boxz 0.0
    set alpha 90.0
    set beta  90.0
    set gamma 90.0
    set triclinic 0
    set lx 0.0
    set ly 0.0
    set lz 0.0
    set xy 0.0
    set xz 0.0
    set yz 0.0

    lassign $boxdata boxx boxy boxz alpha beta gamma triclinic lx ly lz xy xz yz

    vmdcon -info "parsing LAMMPS Atoms section with style '$style'."

    set curatoms 0
    while {[gets $fp line] >= 0} {
        incr lineno

        set atomid 0
        set resid 0
        set atomname ""
        set resname ""
        set atomtype 0
        set charge 0.0
        set mass 1.0 ; #  m=0.0 in MD gets us in trouble, so use a different default.
        set radius 1.5 ; # default radius for unknown elements.
        # XXX: we could have a guess(element|radius) utility for setting this to something better.
        set x 0.0
        set y 0.0
        set z 0.0
        set xi 0
        set yi 0
        set zi 0

        if {[regexp {^\s*(\#.*|)$} $line]} {
            # skip empty, whitespace or comment lines.
        } else {
            if {$cgcmm} {
                if {[regexp {^(.*)\#\s*(\S+)(\s+(\S+))?} $line all nline atomname dummy resname]} {
                    set line $nline
                }
            } else {
                if {[regexp {^(.*)\s+\#\s+(\S+)} $line all nline resname]} {
                    set line $nline
                }
            }
            incr curatoms
            switch $style { # XXX: use regexp based parser to detect wrong formats.

                atomic {
                    if {[llength $line] >= 8} {
                        lassign $line atomid       atomtype        x y z xi yi zi
                    } else {
                        lassign $line atomid       atomtype        x y z
                    }
                }

                bond      -
                angle     -
                molecular {
                    if {[llength $line] >= 9} {
                        lassign $line atomid resid atomtype        x y z xi yi zi
                    } else {
                        lassign $line atomid resid atomtype        x y z
                    }
                }

                charge {
                    if {[llength $line] >= 9} {
                        lassign $line atomid atomtype charge x y z xi yi zi
                    } else {
                        lassign $line atomid atomtype charge x y z
                    }
                }

                sphere {
                    if {[llength $line] >= 10} {
                        lassign $line atomid atomtype radius mass x y z xi yi zi
                    } else {
                        lassign $line atomid atomtype radius mass x y z
                    }
                    # sphere has diameter and density instead of radius and mass
                    # convert them accordingly
                    set radius [expr {0.5*$radius}]
                    set mass [expr {4.0/3.0*$M_PI*$radius*$radius*$radius*$mass}]
                }

                full {
                    if {[llength $line] >= 10} {
                        lassign $line atomid resid atomtype charge x y z xi yi zi
                    } else {
                        lassign $line atomid resid atomtype charge x y z
                    }
                }

                default   {
                    # ignore this unsupported style
                }
            }

            # sanity check on input data. if x/y/z are empty the atom style
            # must have been specified wrong. sufficient to check for $z
            if { [string length $z] == 0 } {
                vmdcon -err "readlammpsatoms: not enough data for style '$style' in line $lineno."
                return -1
            }

            # XXX: no reason to discriminate by atomid. we store the
            # atomid in the user field and assign atoms sorted by atomid
            # (for faster lookup) and can retrieve that info later and
            # use it for mapping bonds, angles, etc to it.
            ################
            # if {$atomid > $numatoms} {
            #     vmdcon -err "readlammpsatoms: only atomids 1-$numatoms are supported. $lineno : $line "
            #    return -1
            # }
            if {$cgcmm} {
                if {[string length $atomname]} {
                    set atomtype $atomname ; # if we have CGCMM data use that.
                } else {
                    set atomname $atomtype
                }
            } else {
                if {[llength $atomlabels] > 0} {
                    set idx [lsearch $atomlabels $atomtype]
                    if {$idx < 0} {
                        vmdcon -err "readlammpsatoms: type $atomtype not found in atom type labelmap. $lineno : $line "
                        return -1
                    }
                    if {[expr {$idx % 2}]} {
                        # search found symbolic type
                        set atomname $atomtype
                    } else {
                        # search found numeric type, replace with symbolic type from map
                        incr idx
                        set atomtype [lindex $atomlabels $idx]
                        set atomname $atomtype
                    }
                } else {
                    set atomname $atomtype
                }
            }
            if {$triclinic} {
                lappend atomdata [list $atomid $resid $resname $atomname $atomtype $charge \
                                      [expr {$x + $xi*$lx + $yi*$xy + $zi*$xz}] \
                                      [expr {$y + $yi*$ly + $zi*$yz}] \
                                      [expr {$z + $zi*$lz}] $mass $radius ]
            } else {
                lappend atomdata [list $atomid $resid $resname $atomname $atomtype $charge \
                                      [expr {$xi*$boxx + $x}] [expr {$yi*$boxy + $y}] \
                                      [expr {$zi*$boxz + $z}] $mass $radius ]
            }
        }
        if {$curatoms >= $numatoms} break
    }
    vmdcon -info "applying atoms data. sorted by atom id."
    $sel set {user resid resname name type charge x y z mass radius} \
        [lsort -integer -index 0 $atomdata]
    return $lineno
}

# parse masses section
proc ::TopoTools::readlammpsmasses {fp mol numtypes cgcmm massmap typemap lineno} {
    vmdcon -info "parsing LAMMPS Masses section."

    upvar $massmap massdata
    upvar $typemap atomlabels
    set massdata {}
    set curtypes 0
    while {[gets $fp line] >= 0} {
        incr lineno

        set typeid 0
        set mass 0.0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curtypes
            set typename {}
            if {[regexp {^(.*)\#\s*(\S+).*} $line all nline typename]} {
                set line $nline
            }
            if {[regexp {^\s*(\d+)\s+(.*)} $line all typeid mass]} {
                if {($typeid < 1) || ($typeid > $numtypes)} {
                    vmdcon -err "readlammpsmasses: only typeids 1-$numtypes are supported. $lineno : $line "
                    return -1
                }
                # if we have a CGCMM style data file, we have strings for types.
                if {$cgcmm && ([string length $typename] > 0)} {
                    lappend massdata $typename $mass
                } else {
                    lappend massdata $typeid $mass
                }
            } else {
                if {[regexp {^\s*(\S+)\s+(\S+)\s*$} $line all typename mass]} {
                    set idx [lsearch $atomlabels $typename]
                    if {$idx < 0} {
                        vmdcon -err "readlammpsmasses: type $typename not found in atom type labelmap. $lineno : $line "
                        return -1
                    }
                    lappend massdata $typename $mass
                } else {
                    vmdcon -err "readlammpsmasses: incorect format. $lineno : $line "
                return -1
                }
            }
        }
        if {$curtypes >= $numtypes} break
    }
    return $lineno
}

# parse velocities section
proc ::TopoTools::readlammpsvelocities {fp sel lineno} {
    set numatoms [$sel num]
    set velocitydata {}

    vmdcon -info "parsing LAMMPS Velocities section."

    set curatoms 0
    set veldata {}
    while {[gets $fp line] >= 0} {
        incr lineno

        set atomid 0
        set vx 0.0
        set vy 0.0
        set vz 0.0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curatoms
            lassign $line atomid vx vy vz

            #if {$atomid > $numatoms} {
            #    vmdcon -err "readlammpsvelocities: only atomids 1-$numatoms are supported. $lineno : $line "
            #    return -1
            #}
            lappend veldata [list $atomid $vx $vy $vz]
        }
        if {$curatoms >= $numatoms} break
    }
    if { [catch {$sel set {user vx vy vz} [lsort -integer -index 0 $veldata]} errmsg] } {
        vmdcon -warn "readlammpsvelocities: problems assigning velocities. skipping..."
    }
    return $lineno
}


# parse bond section
proc ::TopoTools::readlammpsbonds {fp sel numbonds atomidmap typemap lineno} {
    upvar $typemap bondlabels
    set curbonds 0
    set bonddata {}

    vmdcon -info "parsing LAMMPS Bonds section."
    while {[gets $fp line] >= 0} {
        incr lineno

        set num 0
        set a 0
        set b 0
        set type 0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curbonds
            lassign $line num type a b ;# XXX: use regexp based parser to detect wrong format.
            # map atomid to vmd atom index.
            set aidx [lsearch -sorted -integer $atomidmap $a]
            set bidx [lsearch -sorted -integer $atomidmap $b]
            if {($aidx < 0) || ($bidx < 0)} {
                vmdcon -err "readlammpsbond: bond with non-existent atomid on line: $lineno"
                return -1
            }
            if {[llength $bondlabels] > 0} {
                set idx [lsearch $bondlabels $type]
                if {$idx < 0} {
                    vmdcon -err "readlammpsbonds: type $type not found in bond type labelmap. $lineno : $line "
                    return -1
                }
                if {[expr {$idx % 2}] == 0} {
                    # search found numeric type, replace with symbolic type from map
                    incr idx
                    set type [lindex $bondlabels $idx]
                }
            }
            lappend bonddata [list $aidx $bidx $type]
        }
        if {$curbonds >= $numbonds} break
    }
    vmdcon -info "applying bonds data."
    setbondlist $sel type $bonddata
    return $lineno
}

# parse angle section
proc ::TopoTools::readlammpsangles {fp sel numangles atomidmap typemap lineno} {
    upvar $typemap anglelabels
    set curangles 0
    set angledata {}

    vmdcon -info "parsing LAMMPS Angles section."
    while {[gets $fp line] >= 0} {
        incr lineno

        set num 0
        set a 0
        set b 0
        set c 0
        set type 0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curangles
            lassign $line num type a b c ;# XXX: use regexp based parser to detect wrong format.
            # map atomid to vmd atom index.
            set aidx [lsearch -sorted -integer $atomidmap $a]
            set bidx [lsearch -sorted -integer $atomidmap $b]
            set cidx [lsearch -sorted -integer $atomidmap $c]
            if {($aidx < 0) || ($bidx < 0) || ($cidx < 0)} {
                vmdcon -err "readlammpsangles: angle with non-existent atomid on line: $lineno"
                return -1
            }
            if {[llength $anglelabels] > 0} {
                set idx [lsearch $anglelabels $type]
                if {$idx < 0} {
                    vmdcon -err "readlammpangles: type $type not found in angle type labelmap. $lineno : $line "
                    return -1
                }
                if {[expr {$idx % 2}] == 0} {
                    # search found numeric type, replace with symbolic type from map
                    incr idx
                    set type [lindex $anglelabels $idx]
                }
            }
            lappend angledata [list $type $aidx $bidx $cidx]
        }
        if {$curangles >= $numangles} break
    }
    vmdcon -info "applying angles data."
    setanglelist $sel $angledata
    return $lineno
}

# parse dihedral section
proc ::TopoTools::readlammpsdihedrals {fp sel numdihedrals atomidmap typemap lineno} {
    upvar $typemap dihedrallabels
    set curdihedrals 0
    set dihedraldata {}

    vmdcon -info "parsing LAMMPS Dihedrals section."
    while {[gets $fp line] >= 0} {
        incr lineno

        set num 0
        set a 0
        set b 0
        set c 0
        set d 0
        set type 0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curdihedrals
            lassign $line num type a b c d ;# XXX: use regexp based parser to detect wrong format.
            # map atomid to vmd atom index.
            set aidx [lsearch -sorted -integer $atomidmap $a]
            set bidx [lsearch -sorted -integer $atomidmap $b]
            set cidx [lsearch -sorted -integer $atomidmap $c]
            set didx [lsearch -sorted -integer $atomidmap $d]
            if {($aidx < 0) || ($bidx < 0) || ($cidx < 0) || ($didx < 0)} {
                vmdcon -err "readlammpsdihedrals: dihedral with non-existent atomid on line: $lineno"
                return -1
            }
            if {[llength $dihedrallabels] > 0} {
                set idx [lsearch $dihedrallabels $type]
                if {$idx < 0} {
                    vmdcon -err "readlammpsdihedrals: type $type not found in dihedral type labelmap. $lineno : $line "
                    return -1
                }
                if {[expr {$idx % 2}] == 0} {
                    # search found numeric type, replace with symbolic type from map
                    incr idx
                    set type [lindex $dihedrallabels $idx]
                }
            }
            lappend dihedraldata [list $type $aidx $bidx $cidx $didx]
        }
        if {$curdihedrals >= $numdihedrals} break
    }
    vmdcon -info "applying dihedrals data."
    setdihedrallist $sel $dihedraldata
    return $lineno
}

# parse improper section
proc ::TopoTools::readlammpsimpropers {fp sel numimpropers atomidmap typemap lineno} {
    upvar $typemap improperlabels
    set curimpropers 0
    set improperdata {}

    vmdcon -info "parsing LAMMPS Impropers section."
    while {[gets $fp line] >= 0} {
        incr lineno

        set num 0
        set a 0
        set b 0
        set c 0
        set d 0
        set type 0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curimpropers
            lassign $line num type a b c d;# XXX: use regexp based parser to detect wrong format.
            # map atomid to vmd atom index.
            set aidx [lsearch -sorted -integer $atomidmap $a]
            set bidx [lsearch -sorted -integer $atomidmap $b]
            set cidx [lsearch -sorted -integer $atomidmap $c]
            set didx [lsearch -sorted -integer $atomidmap $d]
            if {($aidx < 0) || ($bidx < 0) || ($cidx < 0) || ($didx < 0)} {
                vmdcon -err "readlammpsimpropers: improper with non-existent atomid on line: $lineno"
                return -1
            }
            if {[llength $improperlabels] > 0} {
                set idx [lsearch $improperlabels $type]
                if {$idx < 0} {
                    vmdcon -err "readlammpsimpropers: type $type not found in improper type labelmap. $lineno : $line "
                    return -1
                }
                if {[expr {$idx % 2}] == 0} {
                    # search found numeric type, replace with symbolic type from map
                    incr idx
                    set type [lindex $improperlabels $idx]
                }
            }
            lappend improperdata [list $type $aidx $bidx $cidx $didx]
        }
        if {$curimpropers >= $numimpropers} break
    }
    vmdcon -info "applying impropers data."
    setimproperlist $sel $improperdata
    return $lineno
}

# parse type label section
proc ::TopoTools::readlammpslabels {fp mol style numtypes typemap lineno} {
    vmdcon -info "parsing LAMMPS ${style} Type Labels section."

    upvar $typemap typedata
    set typedata {}
    set curtypes 0
    while {[gets $fp line] >= 0} {
        incr lineno

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curtypes
            if {[regexp {^(.*)\s+\#.*} $line all nline]} {
                set line $nline
            }
            if {[regexp {^\s*(\d+)\s+(\S+)\s*$} $line all typeid label]} {
                if {($typeid < 1) || ($typeid > $numtypes)} {
                    vmdcon -err "readlammpslabels: only typeids 1-$numtypes are supported. $lineno : $line "
                    return -1
                }
                lappend typedata $typeid $label
            } else {
                vmdcon -err "readlammpslabels: incorect format. $lineno : $line "
                return -1
            }
        }
        if {$curtypes >= $numtypes} break
    }
    return $lineno
}


# export internal structure data to a LAMMPS data file.
# this requires that a corresponding set of information
# is already present in VMD's memory.
# Arguments:
# mol = molecule id with matching coordinate data
# filename = name of data file
# typelabels = write data file with type labels instead of numeric types
# style = atom style
# sel = selection function for the subset to be written out.
# flags = more flags. (currently not used)
proc ::TopoTools::writelammpsdata {mol filename typelabels style sel {flags none}} {
    if {[catch {open $filename w} fp]} {
        vmdcon -err "writelammpsdata: problem opening data file: $fp\n"
        return -1
    }

    # initialize system default settings
    array set lammps {
        atoms 0 atomtypes 0 bonds 0 bondtypes 0 angles 0 angletypes 0
        dihedrals 0 dihedraltypes 0 impropers 0 impropertypes 0 xtrabond 0
        xlo 0 xhi 0 ylo 0 yhi 0 zlo 0 zhi 0 xy 0 xz 0 yz 0 triclinic 0
        style unknown typemass 1 typelabels 0
    }

    # gather available system information
    set lammps(atoms)         [$sel num]
    set lammps(bonds)         [bondinfo     numbonds     $sel]
    set lammps(angles)        [angleinfo    numangles    $sel]
    set lammps(dihedrals)     [dihedralinfo numdihedrals $sel]
    set lammps(impropers)     [improperinfo numimpropers $sel]
    set lammps(atomtypes)     [llength [lsort -ascii -unique [$sel get {type}]]]
    set lammps(bondtypes)     [bondinfo     numbondtypes     $sel]
    set lammps(angletypes)    [angleinfo    numangletypes    $sel]
    set lammps(dihedraltypes) [dihedralinfo numdihedraltypes $sel]
    set lammps(impropertypes) [improperinfo numimpropertypes $sel]
    set lammps(style) $style
    set lammps(typelabels) $typelabels

    # correct system information to allow only information valid
    # for the selected atom style
    switch $style {
        atomic -
        charge {
            set lammps(bonds) 0
            set lammps(angles) 0
            set lammps(dihedrals) 0
            set lammps(impropers) 0
            set lammps(bondtypes) 0
            set lammps(angletypes) 0
            set lammps(dihedraltypes) 0
            set lammps(impropertypes) 0
        }

        sphere {
            set lammps(typemass) 0
            set lammps(bonds) 0
            set lammps(angles) 0
            set lammps(dihedrals) 0
            set lammps(impropers) 0
            set lammps(bondtypes) 0
            set lammps(angletypes) 0
            set lammps(dihedraltypes) 0
            set lammps(impropertypes) 0
        }

        bond  {
            set lammps(angles) 0
            set lammps(dihedrals) 0
            set lammps(impropers) 0
            set lammps(angletypes) 0
            set lammps(dihedraltypes) 0
            set lammps(impropertypes) 0
        }

        angle  {
            set lammps(dihedrals) 0
            set lammps(impropers) 0
            set lammps(dihedraltypes) 0
            set lammps(impropertypes) 0
        }

        molecular -
        full -
        default { ; # print all sections
        }
    }

    # initialize simulation cell dimensions from min/max search
    lassign [measure minmax $sel -withradii] min max
    lassign $min xlo ylo zlo
    lassign $max xhi yhi zhi
    lassign [molinfo $mol get {a b c alpha beta gamma}] boxx boxy boxz alpha beta gamma

    # adjust min/max settings when box info is available.
    # try to (mostly) preserve the center of the cell, by
    # deriving it from the preceeding min/max search.
    # refuse to write data file without a box.
    set small 0.0001
    if {($boxx < $small) || ($boxy < $small) || ($boxz < $small)} {
        vmdcon -err "writelammpsdata: need to have non-zero box sizes to write a data file"
        vmdcon -err "writelammpsdata: current box sizes: {$boxx $boxy $boxz}"
        return -1
    }
    set lammps(xmid) [expr {($xlo + $xhi) * 0.5}]
    set lammps(xlo)  [expr {-0.5*$boxx + $lammps(xmid)}]
    set lammps(xhi)  [expr { 0.5*$boxx + $lammps(xmid)}]
    set lammps(ymid) [expr {($ylo + $yhi) * 0.5}]
    set lammps(ylo)  [expr {-0.5*$boxy + $lammps(ymid)}]
    set lammps(yhi)  [expr { 0.5*$boxy + $lammps(ymid)}]
    set lammps(zmid) [expr {($zlo + $zhi) * 0.5}]
    set lammps(zlo)  [expr {-0.5*$boxz + $lammps(zmid)}]
    set lammps(zhi)  [expr { 0.5*$boxz + $lammps(zmid)}]

    # if angle is not set assume orthogonal.
    if {$alpha == 0.0} { set alpha 90.0 }
    if {$beta  == 0.0} { set beta  90.0 }
    if {$gamma == 0.0} { set gamma 90.0 }

    if {($alpha != 90.0) || ($beta != 90.0) || ($gamma != 90.0)} {
        set conv 0.01745329251994329576; # pi/180.0
        set lammps(triclinic) 1
        set lammps(xy) [expr {$boxy * cos($gamma*$conv)}]
        set lammps(xz) [expr {$boxz * cos($beta*$conv)}]
        set lammps(yz) 0.0
        set ly [expr {sqrt($boxy*$boxy - $lammps(xy)*$lammps(xy))}]
        if {abs($ly) > $small} {
            set lammps(yz) [expr {($boxy*$boxz*cos($alpha*$conv)
                                   - $lammps(xy)*$lammps(xz)) / $ly}]
        }
        set lz [expr {sqrt($boxz*$boxz - $lammps(xz)*$lammps(xz) - $lammps(yz)*$lammps(yz))}]
        # update y/z-box boundaries for tilt
        if {$ly > $small} {
            set lammps(ylo)  [expr {-0.5*$ly + $lammps(ymid)}]
            set lammps(yhi)  [expr { 0.5*$ly + $lammps(ymid)}]
        }
        if {$lz > $small} {
            set lammps(zlo)  [expr {-0.5*$lz + $lammps(zmid)}]
            set lammps(zhi)  [expr { 0.5*$lz + $lammps(zmid)}]
        }
    }

    # write out supported data file sections
    writelammpsheader $fp [array get lammps]

    if {$typelabels} {
        writelammpslabelmaps $fp $sel [array get lammps]
    }

    # write out hints about type to number mappings
    # for coefficient settings
    writelammpscoeffhint $fp $sel $typelabels atoms
    if {$lammps(bonds) > 0} {
        writelammpscoeffhint $fp $sel $typelabels bonds
    }
    if {$lammps(angles) > 0} {
        writelammpscoeffhint $fp $sel $typelabels angles
    }
    if {$lammps(dihedrals) > 0} {
        writelammpscoeffhint $fp $sel $typelabels dihedrals
    }
    if {$lammps(impropers) > 0} {
        writelammpscoeffhint $fp $sel $typelabels impropers
    }
    if {$lammps(typemass) > 0} {
        writelammpsmasses $fp $sel $typelabels
    }
    writelammpsatoms $fp $sel $style $typelabels
    set atomidmap  [$sel list]
    if {$lammps(bonds) > 0} {
        writelammpsbonds $fp $sel $atomidmap $typelabels
    }
    if {$lammps(angles) > 0} {
        writelammpsangles $fp $sel $atomidmap $typelabels
    }
    if {$lammps(dihedrals) > 0} {
        writelammpsdihedrals $fp $sel $atomidmap $typelabels
    }
    if {$lammps(impropers) > 0} {
        writelammpsimpropers $fp $sel $atomidmap $typelabels
    }
    close $fp
    return 0
}


# write lammps header section to open file
proc ::TopoTools::writelammpsheader {fp flags} {
    variable version
    array set lammps $flags
    # first header line is skipped.
    if {$lammps(typelabels)} {
        puts $fp "LAMMPS data file. LABELMAP style. atom_style $lammps(style) generated by VMD/TopoTools v$version on [clock format [clock seconds]]"
    } else {
        puts $fp "LAMMPS data file. CGCMM style. atom_style $lammps(style) generated by VMD/TopoTools v$version on [clock format [clock seconds]]"
    }
    foreach key {atoms bonds angles dihedrals impropers} {
        puts $fp [format " %d %s" $lammps($key) $key]
    }

    foreach key {atomtypes bondtypes  angletypes dihedraltypes impropertypes} {
        puts $fp [format " %d %s" $lammps($key) [regsub types $key " &"]]
    }

    puts $fp [format " %.6f %.6f  xlo xhi" $lammps(xlo) $lammps(xhi)]
    puts $fp [format " %.6f %.6f  ylo yhi" $lammps(ylo) $lammps(yhi)]
    puts $fp [format " %.6f %.6f  zlo zhi" $lammps(zlo) $lammps(zhi)]

    if {$lammps(triclinic)} {
        puts $fp [format " %.6f %.6f %.6f xy xz yz" $lammps(xy) $lammps(xz) $lammps(yz)]
    }

    puts $fp ""
    return
}

# write lammps type label maps to file
proc ::TopoTools::writelammpslabelmaps {fp sel flags} {
    variable version
    array set lammps $flags

    if {$lammps(atomtypes) > 0} {
        set typemap [lsort -unique -ascii [$sel get type]]
        set typeid 1

        puts $fp " Atom Type Labels\n"
        foreach type $typemap {
            puts $fp [format " %d %s" $typeid $type]
            incr typeid
        }
    }
    puts $fp ""

    if {$lammps(bonds) > 0} {
        puts $fp " Bond Type Labels\n"
        set bid 1
        foreach bt [bondinfo bondtypenames $sel type] {
            puts $fp " $bid  $bt"
            incr bid
        }
        puts $fp ""
    }
    if {$lammps(angles) > 0} {
        puts $fp " Angle Type Labels\n"
        set aid 1
        foreach at [angleinfo angletypenames $sel] {
            puts $fp " $aid  $at"
            incr aid
        }
        puts $fp ""
    }
    if {$lammps(dihedrals) > 0} {
        puts $fp " Dihedral Type Labels\n"
        set did 1
        foreach dt [dihedralinfo dihedraltypenames $sel] {
            puts $fp " $did  $dt"
            incr did
        }
        puts $fp ""
    }
    if {$lammps(impropers) > 0} {
        puts $fp " Improper Type Labels\n"
        set iid 1
        foreach it [improperinfo impropertypenames $sel] {
            puts $fp " $iid  $it"
            incr iid
        }
        puts $fp ""
    }
    return
}

# write masses section, but only if number of masses
# matches the number of atom types and if no mass is < 0.01
proc ::TopoTools::writelammpsmasses {fp sel typelabels} {

    # first run the checks and build list of masses
    set typemap  [lsort -unique -ascii [$sel get type]]
    set masslist {}
    set mol [$sel molid]
    set selstr [$sel text]
    foreach type $typemap {
        set tsel [atomselect $mol "( $selstr ) and (type '$type')"]
        set mass [lsort -unique -real [$tsel get mass]]
        $tsel delete
        if {[llength $mass] < 1} return
        if {$mass < 0.01} return
        lappend masslist [lindex $mass 0]
    }

    # we passed the test, write out what we learned.
    vmdcon -info "writing LAMMPS Masses section."

    puts $fp " Masses\n"
    set typeid 1
    foreach mass $masslist type $typemap {
        if {$typelabels} {
            puts $fp [format " %s %.6f" $type $mass]
        } else {
            puts $fp [format " %d %.6f \# %s" $typeid $mass $type]
        }
        incr typeid
    }
    puts $fp ""
    return
}

# write atoms section
proc ::TopoTools::writelammpsatoms {fp sel style typelabels} {
    global M_PI

    vmdcon -info "writing LAMMPS Atoms section in style '$style'."

    puts $fp " Atoms # $style\n"
    set typemap [lsort -unique -ascii [$sel get type]]
    set resmap  [lsort -unique -integer [$sel get resid]]
    set atomid 0
    foreach adat [$sel get {type resid charge x y z resname mass radius}] {
        lassign $adat type resid charge x y z resname mass radius
        set atomtype [lsearch -sorted -ascii $typemap $type]
        set resid    [lsearch -sorted -integer $resmap $resid]
        incr atomid
        incr atomtype
        incr resid
        switch $style {
            atomic {
                if {$typelabels} {
                    puts $fp [format "%d %s %.6f %.6f %.6f" \
                                  $atomid        $type  $x $y $z]
                } else {
                    puts $fp [format "%d %d %.6f %.6f %.6f \# %s" \
                                  $atomid        $atomtype  $x $y $z $type]
                }
            }

            bond  -
            angle -
            molecular {
                if {$typelabels} {
                    puts $fp [format "%d %d %s %.6f %.6f %.6f \# %s" \
                                  $atomid $resid $type  $x $y $z $resname]
                } else {
                    puts $fp [format "%d %d %s %.6f %.6f %.6f \# %s" \
                                  $atomid $resid $atomtype  $x $y $z $resname]
                }
            }

            charge    {
                if {$typelabels} {
                    puts $fp [format "%d %s %.6f %.6f %.6f %.6f" \
                                  $atomid $type $charge $x $y $z]
                } else {
                    puts $fp [format "%d %d %.6f %.6f %.6f %.6f \# %s" \
                                  $atomid $atomtype $charge $x $y $z $type]
                }
            }

            sphere {
                # sphere has diameter and density instead of radius and mass
                # convert them accordingly
                set mass [expr {$mass/(4.0/3.0*$M_PI*$radius*$radius*$radius)}]
                set radius [expr {2.0*$radius}]
                if {$typelabels} {
                    puts $fp [format "%d %s %.6f %.6f %.6f %.6f %.6f" \
                                  $atomid $type $radius $mass $x $y $z]
                } else {
                    puts $fp [format "%d %d %.6f %.6f %.6f %.6f %.6f \# %s" \
                                  $atomid $atomtype $radius $mass $x $y $z $type]
                }
            }

            full      {
                if {$typelabels} {
                    puts $fp [format "%d %d %s %.6f %.6f %.6f %.6f \# %s" \
                                  $atomid $resid $type $charge $x $y $z $resname]
                } else {
                    puts $fp [format "%d %d %d %.6f %.6f %.6f %.6f \# %s %s" \
                                  $atomid $resid $atomtype $charge $x $y $z $type $resname]
                }
            }

            default   {
                # ignore this unsupported style
                # XXX: add a way to flag an error. actually the test for
                #      supported lammps atom styles should be done on a
                #      much higher level, so that we don't do unneeded work.
            }
        }
    }
    puts $fp ""
    return
}

# write bond section
proc ::TopoTools::writelammpsbonds {fp sel atomidmap typelabels} {
    set bonddata  [bondinfo getbondlist   $sel type]
    set bondtypes [bondinfo bondtypenames $sel type]
    vmdcon -info "writing LAMMPS Bonds section."
    puts $fp " Bonds\n"

    set bondid 0
    foreach bdat $bonddata {
        incr bondid
        lassign $bdat a b t
        set at1  [lsearch -integer -sorted $atomidmap $a]
        set at2  [lsearch -integer -sorted $atomidmap $b]
        # go from 0-based to 1-based indexing
        incr at1; incr at2;

        if {$typelabels} {
            puts $fp [format "%d %s %d %d" $bondid $t $at1 $at2]
        } else {
            set type [lsearch -ascii $bondtypes $t]

            # go from 0-based to 1-based indexing and write out
            incr type
            puts $fp [format "%d %d %d %d" $bondid $type $at1 $at2]
        }
    }
    puts $fp ""
    return
}

# write angle section
proc ::TopoTools::writelammpsangles {fp sel atomidmap typelabels} {
    set angledata  [angleinfo getanglelist   $sel]
    set angletypes [angleinfo angletypenames $sel]
    vmdcon -info "writing LAMMPS Angles section."
    puts $fp " Angles\n"

    set angleid 0
    foreach adat $angledata {
        incr angleid
        lassign $adat t a b c
        set at1  [lsearch -integer -sorted $atomidmap $a]
        set at2  [lsearch -integer -sorted $atomidmap $b]
        set at3  [lsearch -integer -sorted $atomidmap $c]
        # go from 0-based to 1-based indexing
        incr at1; incr at2; incr at3

        if {$typelabels} {
            puts $fp [format "%d %s %d %d %d" $angleid $t $at1 $at2 $at3]
        } else {
            set type [lsearch -ascii $angletypes $t]

            # go from 0-based to 1-based indexing and write out
            incr type
            puts $fp [format "%d %d %d %d %d" $angleid $type $at1 $at2 $at3]
        }
    }
    puts $fp ""
    return
}

# write dihedral section
proc ::TopoTools::writelammpsdihedrals {fp sel atomidmap typelabels} {
    set dihedraldata  [dihedralinfo getdihedrallist   $sel]
    set dihedraltypes [dihedralinfo dihedraltypenames $sel]
    vmdcon -info "writing LAMMPS Dihedrals section."
    puts $fp " Dihedrals\n"

    set dihedralid 0
    foreach adat $dihedraldata {
        incr dihedralid
        lassign $adat t a b c d
        set at1  [lsearch -integer -sorted $atomidmap $a]
        set at2  [lsearch -integer -sorted $atomidmap $b]
        set at3  [lsearch -integer -sorted $atomidmap $c]
        set at4  [lsearch -integer -sorted $atomidmap $d]
        # go from 0-based to 1-based indexing
        incr at1; incr at2; incr at3; incr at4

        if {$typelabels} {
            puts $fp [format "%d %s %d %d %d %d" $dihedralid $t $at1 $at2 $at3 $at4]
        } else {
            set type [lsearch -ascii $dihedraltypes $t]

            # go from 0-based to 1-based indexing and write out
            incr type
            puts $fp [format "%d %d %d %d %d %d" $dihedralid $type $at1 $at2 $at3 $at4]
        }
    }
    puts $fp ""
    return
}

# write improper section
proc ::TopoTools::writelammpsimpropers {fp sel atomidmap typelabels} {
    set improperdata  [improperinfo getimproperlist   $sel]
    set impropertypes [improperinfo impropertypenames $sel]
    vmdcon -info "writing LAMMPS Impropers section."
    puts $fp " Impropers\n"

    set improperid 0
    foreach adat $improperdata {
        incr improperid
        lassign $adat t a b c d
        set at1  [lsearch -integer -sorted $atomidmap $a]
        set at2  [lsearch -integer -sorted $atomidmap $b]
        set at3  [lsearch -integer -sorted $atomidmap $c]
        set at4  [lsearch -integer -sorted $atomidmap $d]
        # go from 0-based to 1-based indexing
        incr at1; incr at2; incr at3; incr at4

        if {$typelabels} {
            puts $fp [format "%d %s %d %d %d %d" $improperid $t $at1 $at2 $at3 $at4]
        } else {
            set type [lsearch -ascii $impropertypes $t]

            # go from 0-based to 1-based indexing and write out
            incr type
            puts $fp [format "%d %d %d %d %d %d" $improperid $type $at1 $at2 $at3 $at4]
        }
    }
    puts $fp ""
    return
}

# returns 0 if lammps atom style is supported by topotools and 1 if not.
proc ::TopoTools::checklammpsstyle {style} {
    switch $style {

        auto -
        atomic -
        bond  -
        angle -
        molecular -
        charge -
        sphere -
        full {
            return 0
        }

        default {
            return 1
        }
    }
}

# write hints about type coefficient mappings
proc ::TopoTools::writelammpscoeffhint {fp sel typelabels type} {
    switch $type {
        atoms {
            puts $fp "\# Pair Coeffs\n\#"
            set aid 1
            set atlist [lsort -ascii -unique [$sel get {type}]]
            foreach at $atlist {
                if {$typelabels} {
                    puts $fp "\# $at"
                } else {
                    puts $fp "\# $aid  $at"
                }
                incr aid
            }
        }
        bonds {
            puts $fp "\# Bond Coeffs\n\#"
            set bid 1
            foreach bt [bondinfo bondtypenames $sel type] {
                if {$typelabels} {
                    puts $fp "\# $bt"
                } else {
                    puts $fp "\# $bid  $bt"
                }
                incr bid
            }
        }
        angles {
            puts $fp "\# Angle Coeffs\n\#"
            set aid 1
            foreach at [angleinfo angletypenames $sel] {
                if {$typelabels} {
                    puts $fp "\# $at"
                } else {
                    puts $fp "\# $aid  $at"
                }
                incr aid
            }
        }
        dihedrals {
            puts $fp "\# Dihedral Coeffs\n\#"
            set did 1
            foreach dt [dihedralinfo dihedraltypenames $sel] {
                if {$typelabels} {
                    puts $fp "\# $dt"
                } else {
                    puts $fp "\# $did  $dt"
                }
                incr did
            }
        }
        impropers {
            puts $fp "\# Improper Coeffs\n\#"
            set iid 1
            foreach it [improperinfo impropertypenames $sel] {
                if {$typelabels} {
                    puts $fp "\# $it"
                } else {
                    puts $fp "\# $iid  $it"
                }
                incr iid
            }
        }
        default {
            vmdcon -warn "writelammpscoeffhint: don't know how to write hints for '$type'"
            return 1
        }
    }
    puts $fp ""
    return
}
