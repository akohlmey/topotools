#!/usr/bin/tclsh
# This file is part of TopoTools, a VMD package to simplify 
# manipulating bonds other topology related properties.
#
# Copyright (c) 2009 by Axel Kohlmeyer <akohlmey@cmm.chem.upenn.edu>
#
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
    if {[catch {open $filename r} fp]} {
        vmdcon -error "readlammpsdata: problem opening data file: $fp\n"
        return -1
    }
    
    # parse lammps header section.
    array set lammps [readlammpsheader $fp]
    if {$lammps(atoms) < 1} {
        vmdcon -error "readlammpsdata: failed to parse lammps data header. aborting."
        return -1
    } 

    # create an empty molecule and timestep
    set mol -1
    if {[catch {mol new atoms $lammps(atoms)} mol]} {
        vmdcon -error "readlammpsdata: problem creating empty molecule: $mol"
        return -1
    } else {
        animate dup $mol
    }
    mol rename $mol [file tail $filename]
    set sel [atomselect $mol all]
    set a [expr {$lammps(xhi) - $lammps(xlo)}]
    set b [expr {$lammps(yhi) - $lammps(ylo)}]
    set c [expr {$lammps(zhi) - $lammps(zlo)}]
    set boxdim [list $a $b $c]
    molinfo $mol set {a b c} $boxdim

    # now loop through the file until we find a known section header.
    # then call a subroutine that parses this section. those subroutines
    # all take the current line number as last argument and return the
    # new value, so we can keep track of the current line and print more
    # useful error messages or warnings.
    set lineno $lammps(lineno)
    while {[gets $fp line] >= 0} {
        incr lineno
        if {[regexp {^\s*Atoms} $line ]} {
            set lineno [readlammpsatoms $fp $sel $style $boxdim $lineno]
            if {$lineno < 0} {
                vmdcon -error "readlammpsdata: error reading Atoms section."
                return -1
            }
        } elseif {[regexp {^\s*Velocities} $line ]} {
            set lineno [readlammpsvelocities $fp $sel $lineno]
            if {$lineno < 0} {
                vmdcon -error "readlammpsdata: error reading Velocities section."
                return -1
            }
        } elseif {[regexp {^\s*Masses} $line ]} {
            set lineno [readlammpsmasses $fp $mol $lammps(atomtypes) $lineno]
            if {$lineno < 0} {
                vmdcon -error "readlammpsdata: error reading Masses section."
                return -1
            }
        } elseif {[regexp {^\s*Bonds} $line ]} {
            set lineno [readlammpsbonds $fp $sel $lammps(bonds) $lineno]
            if {$lineno < 0} {
                vmdcon -error "readlammpsdata: error reading Bonds section."
                return -1
            }
        } elseif {[regexp {^\s*Angles} $line ]} {
            set lineno [readlammpsangles $fp $sel $lammps(angles) $lineno]
            if {$lineno < 0} {
                vmdcon -error "readlammpsdata: error reading Angles section."
                return -1
            }
        } elseif {[regexp {^\s*Dihedrals} $line ]} {
            set lineno [readlammpsdihedrals $fp $sel $lammps(dihedrals) $lineno]
            if {$lineno < 0} {
                vmdcon -error "readlammpsdata: error reading Dihedrals section."
                return -1
            }
        } elseif {[regexp {^\s*Impropers} $line ]} {
            set lineno [readlammpsimpropers $fp $sel $lammps(impropers) $lineno]
            if {$lineno < 0} {
                vmdcon -error "readlammpsdata: error reading Impropers section."
                return -1
            }
        } elseif {[regexp {^\s*((Pair|Bond|Angle|Dihedral|Improper) Coeffs)} $line ]} {
            # add code to skip silently.
        } elseif { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else { 
            vmdcon -error "readlammpsdata: unkown header line: $lineno : $line "
            vmdcon -error "readlammpsdata: cannot continue. "
            return -1
        }    
        set lammps(lineno) $lineno
    }
    close $fp
    mol reanalyze $mol
    return $mol
}


# read lammps header section from opened file
# and return as an array.
proc ::TopoTools::readlammpsheader {fp} {
    array set lammps {
        atoms 0 atomtypes 0 bonds 0 bondtypes 0 angles 0 angletypes 0 
        dihedrals 0 dihedraltypes 0 impropers 0 impropertypes 0 xtrabond 0
        xlo 0 xhi 0 ylo 0 yhi 0 zlo 0 zhi 0 xy 0 xz 0 yz 0 
        lineno 0
    }
    set x {}

    vmdcon -info "parsing LAMMPS header."

    # skip first header line
    gets $fp line
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
        } elseif { [regexp {^\s*([-[:digit:].e+]+)\s+([-[:digit:].e+]+)\s+xlo xhi} $line \
                        x lammps(xlo) lammps(xhi)] } {
        } elseif { [regexp {^\s*([-[:digit:].e+]+)\s+([-[:digit:].e+]+)\s+ylo yhi} $line \
                        x lammps(ylo) lammps(yhi)] } {
        } elseif { [regexp {^\s*([-[:digit:].e+]+)\s+([-[:digit:].e+]+)\s+zlo zhi} $line \
                        x lammps(zlo) lammps(zhi)] } {
        } elseif { [regexp {^\s*([-[:digit:].e+]+)\s+([-[:digit:].e+]+)\s+([-[:digit:].e+]+)\s+xlo xhi} $line x lammps(xy) lammps(xz) lammps(yz)] } {
        } elseif { [regexp {^\s*(\#.*|)$} $line ] } {
        } elseif {[regexp {^\s*(Atoms|Velocities|Masses|Shapes|Dipoles|Bonds|Angles|Dihedrals|Impropers|(Pair|Bond|Angle|Dihedral|Improper) Coeffs)} $line ]} {
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
proc ::TopoTools::readlammpsatoms {fp sel style boxdata lineno} {
    set numatoms [$sel num]
    set atomdata {}
    set boxx 0.0
    set boxy 0.0 
    set boxz 0.0
    lassign $boxdata boxx boxy boxz

    vmdcon -info "parsing LAMMPS Atoms section."

    set curatoms 0
    while {[gets $fp line] >= 0} {
        incr lineno

        set atomid 0
        set resid 0
        set atomtype 0
        set charge 0.0
        set mass 1.0 ; #  m=0.0 in MD gets us in trouble, so use a different default.
        set x 0.0
        set y 0.0
        set z 0.0
        set xi 0
        set yi 0
        set zi 0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curatoms
            switch $style { # XXX: use regexp based parser to detect wrong formats.

                atomic    { 
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
                        lassign $line atomid       atomtype charge x y z xi yi zi
                    } else {
                        lassign $line atomid       atomtype charge x y z
                    }
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
            if {$atomid > $numatoms} {
                vmdcon -error "readlammpsatoms: only atomids 1-$numatoms are supported. $lineno : $line "
                return -1
            }
            lappend atomdata [list $atomid $resid $atomtype $atomtype $charge [expr {$xi*$boxx + $x}] \
                                  [expr {$yi*$boxy + $y}] [expr {$zi*$boxz + $z}] $mass ]
        }
        if {$curatoms >= $numatoms} break
    }
    vmdcon -info "applying atoms data."
    $sel set {user resid name type charge x y z mass} [lsort -integer -index 0 $atomdata]
    return $lineno
}

# parse masses section
proc ::TopoTools::readlammpsmasses {fp mol numtypes lineno} {
    vmdcon -info "parsing LAMMPS Masses section."

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
            lassign $line typeid mass
            if {$typeid > $numtypes} {
                vmdcon -error "readlammpsatoms: only typeids 1-$numtypes are supported. $lineno : $line "
                return -1
            }
            set sel [atomselect $mol "type $typeid"]
            $sel set mass $mass
            $sel delete
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
    while {[gets $fp line] >= 0} {
        incr lineno

        set atomid 0
        set vx 0.0
        set vy 0.0
        set vz 0.0
        set veldata {}

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curatoms
            lassign $line atomid vx vy vz 

            if {$atomid > $numatoms} {
                vmdcon -error "readlammpsvelocities: only atomids 1-$numatoms are supported. $lineno : $line "
                return -1
            }
            lappend veldata [list $atomid $vx $vy $vz]
        }
        if {$curatoms >= $numatoms} break
    }
    $sel set {user vx vy vz} [lsort -integer -index 0 $veldata]
    return $lineno
}


# parse bond section
proc ::TopoTools::readlammpsbonds {fp sel numbonds lineno} {
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
            incr a -1 ; # we force that only atomids 1-numatoms are used (normal case)
            incr b -1 ; # so this does work to translate into VMD numbering.
            lappend bonddata [list $a $b $type]
        }
        if {$curbonds >= $numbonds} break
    }
    vmdcon -info "applying bonds data."
    setbondlist $sel type $bonddata
    return $lineno
}

# parse angle section
proc ::TopoTools::readlammpsangles {fp sel numangles lineno} {
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
            incr a -1 ; # we force that only atomids 1-numatoms are used (normal case)
            incr b -1 ; # so this does work to translate into VMD numbering.
            incr c -1
            lappend angledata [list $type $a $b $c]
        }
        if {$curangles >= $numangles} break
    }
    vmdcon -info "applying angles data."
    setanglelist $sel $angledata
    return $lineno
}

# parse dihedral section
proc ::TopoTools::readlammpsdihedrals {fp sel numdihedrals lineno} {
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
            incr a -1 ; # we force that only atomids 1-numatoms are used (normal case)
            incr b -1 ; # so this does work to translate into VMD numbering.
            incr c -1
            incr d -1
            lappend dihedraldata [list $type $a $b $c $d]
        }
        if {$curdihedrals >= $numdihedrals} break
    }
    vmdcon -info "applying dihedrals data."
    setdihedrallist $sel $dihedraldata
    return $lineno
}

# parse improper section
proc ::TopoTools::readlammpsimpropers {fp sel numimpropers lineno} {
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
            incr a -1 ; # we force that only atomids 1-numatoms are used (normal case)
            incr b -1 ; # so this does work to translate into VMD numbering.
            incr c -1
            incr d -1
            lappend improperdata [list $type $a $b $c $d]
        }
        if {$curimpropers >= $numimpropers} break
    }
    vmdcon -info "applying impropers data."
    setimproperlist $sel $improperdata
    return $lineno
}


# export internal structure data to a LAMMPS data file.
# this requires that a corresponding set of information
# is already present in VMD's memory.
# Arguments:
# mol = molecule id with matching coordinate data
# filename = name of data file
# style = atom style
# sel = selection function for the subset to be written out.
# flags = more flags. (currently not used)
proc ::TopoTools::writelammpsdata {mol filename style sel {flags none}} {
    if {[catch {open $filename w} fp]} {
        vmdcon -error "writelammpsdata: problem opening data file: $fp\n"
        return -1
    }

    # initialize system default settings
    array set lammps {
        atoms 0 atomtypes 0 bonds 0 bondtypes 0 angles 0 angletypes 0 
        dihedrals 0 dihedraltypes 0 impropers 0 impropertypes 0 xtrabond 0
        xlo 0 xhi 0 ylo 0 yhi 0 zlo 0 zhi 0 xy 0 xz 0 yz 0 
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

    # initialize simulation cell dimensions from min/max search
    lassign [measure minmax $sel -withradii] min max
    lassign $min xlo ylo zlo
    lassign $max xhi yhi zhi
    lassign [molinfo $mol get {a b c alpha beta gamma}] boxx boxy boxz alpha beta gamma
    
    # override min/max settings where box info available. 
    # try to (mostly) preserve the center of the cell, by
    # deriving it from the preceeding min/max search.
    set small 0.0001
    if {$boxx > $small} {
        set lammps(xmid) [expr {($xlo + $xhi) * 0.5}]
        set lammps(xlo)  [expr {-0.5*$boxx + $lammps(xmid)}]
        set lammps(xhi)  [expr { 0.5*$boxx + $lammps(xmid)}]
    }
    if {$boxy > $small} {
        set lammps(ymid) [expr {($ylo + $yhi) * 0.5}]
        set lammps(ylo)  [expr {-0.5*$boxy + $lammps(ymid)}]
        set lammps(yhi)  [expr { 0.5*$boxy + $lammps(ymid)}]
    }
    if {$boxz > $small} {
        set lammps(zmid) [expr {($zlo + $zhi) * 0.5}]
        set lammps(zlo)  [expr {-0.5*$boxz + $lammps(zmid)}]
        set lammps(zhi)  [expr { 0.5*$boxz + $lammps(zmid)}]
    }
    # XXX: need check for non-orthogonal cells.

    # write out supported data file sections
    writelammpsheader $fp [array get lammps]
    writelammpsatoms $fp $sel $style
    set atomidmap  [$sel get serial]
    if {$lammps(bonds) > 0} {
        writelammpsbonds $fp $sel $atomidmap
    }
    if {$lammps(angles) > 0} {
        writelammpsangles $fp $sel $atomidmap
    }
    if {$lammps(dihedrals) > 0} {
        writelammpsdihedrals $fp $sel $atomidmap
    }
    if {$lammps(impropers) > 0} {
        writelammpsimpropers $fp $sel $atomidmap
    }
    close $fp
    return 0
}


# write lammps header section to open file
proc ::TopoTools::writelammpsheader {fp flags} {
    variable version
    array set lammps $flags
    # first header line is skipped.
    puts $fp "LAMMPS data file. CGCMM style. generated by VMD/TopoTools v$version on [clock format [clock seconds]]"

    foreach key {atoms bonds angles dihedrals impropers} {
        puts $fp [format " %d %s" $lammps($key) $key]
    }

    foreach key {atomtypes bondtypes  angletypes dihedraltypes impropertypes} {
        puts $fp [format " %d %s" $lammps($key) [regsub types $key " &"]]
    }

    puts $fp [format " %.4f %.4f  xlo xhi" $lammps(xlo) $lammps(xhi)]
    puts $fp [format " %.4f %.4f  ylo yhi" $lammps(ylo) $lammps(yhi)]
    puts $fp [format " %.4f %.4f  zlo zhi" $lammps(zlo) $lammps(zhi)]

    puts $fp ""
    return
}

# write atoms section
proc ::TopoTools::writelammpsatoms {fp sel style} {

    vmdcon -info "writing LAMMPS Atoms section in style '$style'."

    puts $fp " Atoms\n"
    set typemap [lsort -unique -ascii [$sel get type]]
    set resmap  [lsort -unique -ascii [$sel get residue]]
    set atomid 0
    foreach adat [$sel get {type residue charge x y z resname}] {
        lassign $adat type residue charge x y z resname
        incr atomid
        set atomtype [expr 1 + [lsearch -sorted -ascii $typemap $type]]
        set resid    [expr 1 + [lsearch -sorted -ascii $resmap $residue]]
        switch $style {
            atomic    { 
                puts $fp [format "%d %d %.3f %.3f %.3f \# %s" \
                              $atomid        $atomtype  $x $y $z $type] 
            }
            bond  -
            angle -
            molecular { 
                puts $fp [format "%d %d %d %.3f %.3f %.3f \# %s %s" \
                              $atomid $resid $atomtype  $x $y $z $type $resname] 
            }
            charge    { 
                puts $fp [format "%d %d %.2f %.3f %.3f %.3f \# %s" \
                              $atomid $atomtype $charge $x $y $z $type] 
            }
            full      { 
                puts $fp [format "%d %d %d %.2f %.3f %.3f %.3f \# %s %s" \
                              $atomid $resid $atomtype $charge $x $y $z $type $resname] 
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

# parse atom section
proc ::TopoTools::writelammpsatoms {fp sel style} {

    vmdcon -info "writing LAMMPS Atoms section in style '$style'."

    puts $fp " Atoms\n"
    set typemap [lsort -unique -ascii [$sel get type]]
    set resmap  [lsort -unique -ascii [$sel get residue]]
    set atomid 0
    foreach adat [$sel get {type residue charge x y z resname}] {
        lassign $adat type residue charge x y z resname
        incr atomid
        set atomtype [expr 1 + [lsearch -sorted -ascii $typemap $type]]
        set resid    [expr 1 + [lsearch -sorted -ascii $resmap $residue]]
        switch $style {
            atomic    { 
                puts $fp [format "%d %d %.3f %.3f %.3f \# %s" \
                              $atomid        $atomtype  $x $y $z $type] 
            }
            bond  -
            angle -
            molecular { 
                puts $fp [format "%d %d %d %.3f %.3f %.3f \# %s %s" \
                              $atomid $resid $atomtype  $x $y $z $type $resname] 
            }
            charge    { 
                puts $fp [format "%d %d %.2f %.3f %.3f %.3f \# %s" \
                              $atomid $atomtype $charge $x $y $z $type] 
            }
            full      { 
                puts $fp [format "%d %d %d %.2f %.3f %.3f %.3f \# %s %s" \
                              $atomid $resid $atomtype $charge $x $y $z $type $resname] 
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
proc ::TopoTools::writelammpsbonds {fp sel atomidmap} {
    set bonddata  [bondinfo getbondlist   $sel type]
    set bondtypes [bondinfo bondtypenames $sel type]
    vmdcon -info "writing LAMMPS Bonds section."
    puts $fp " Bonds\n"

    set bondid 0
    foreach bdat $bonddata {
        incr bondid
        lassign $bdat a b t
        set at1 [lindex $atomidmap $a]
        set at2 [lindex $atomidmap $b]
        set type [expr 1 + [lsearch -ascii $bondtypes $t]]
   
        puts $fp [format "%d %d %d %d" $bondid $type $at1 $at2]
    }
    puts $fp ""
    return
}

# write angle section
proc ::TopoTools::writelammpsangles {fp sel atomidmap} {
    set angledata  [angleinfo getanglelist   $sel]
    set angletypes [angleinfo angletypenames $sel]
    vmdcon -info "writing LAMMPS Angles section."
    puts $fp " Angles\n"

    set angleid 0
    foreach adat $angledata {
        incr angleid
        lassign $adat t a b c
        set at1 [lindex $atomidmap $a]
        set at2 [lindex $atomidmap $b]
        set at3 [lindex $atomidmap $c]
        set type [expr 1 + [lsearch -ascii $angletypes $t]]
   
        puts $fp [format "%d %d %d %d %d" $angleid $type $at1 $at2 $at3]
    }
    puts $fp ""
    return
}

# write dihedral section
proc ::TopoTools::writelammpsdihedrals {fp sel atomidmap} {
    set dihedraldata  [dihedralinfo getdihedrallist   $sel]
    set dihedraltypes [dihedralinfo dihedraltypenames $sel]
    vmdcon -info "writing LAMMPS Dihedrals section."
    puts $fp " Dihedrals\n"

    set dihedralid 0
    foreach adat $dihedraldata {
        incr dihedralid
        lassign $adat t a b c d
        set at1 [lindex $atomidmap $a]
        set at2 [lindex $atomidmap $b]
        set at3 [lindex $atomidmap $c]
        set at4 [lindex $atomidmap $d]
        set type [expr 1 + [lsearch -ascii $dihedraltypes $t]]
   
        puts $fp [format "%d %d %d %d %d %d" $dihedralid $type $at1 $at2 $at3 $at4]
    }
    puts $fp ""
    return
}

# write improper section
proc ::TopoTools::writelammpsimpropers {fp sel atomidmap} {
    set improperdata  [improperinfo getimproperlist   $sel]
    set impropertypes [improperinfo impropertypenames $sel]
    vmdcon -info "writing LAMMPS Impropers section."
    puts $fp " Impropers\n"

    set improperid 0
    foreach adat $improperdata {
        incr improperid
        lassign $adat t a b c d
        set at1 [lindex $atomidmap $a]
        set at2 [lindex $atomidmap $b]
        set at3 [lindex $atomidmap $c]
        set at4 [lindex $atomidmap $d]
        set type [expr 1 + [lsearch -ascii $impropertypes $t]]
   
        puts $fp [format "%d %d %d %d %d %d" $improperid $type $at1 $at2 $at3 $at4]
    }
    puts $fp ""
    return
}

# returns 0 if lammps atom style is supported by topotools and 1 if not.
proc ::TopoTools::checklammpsstyle {style} {
    switch $style {

        atomic -
        bond  -
        angle -
        molecular -
        charge -
        full {
            return 0
        }

        default {
            return 1
        }
    }
}
