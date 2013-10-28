#!/usr/bin/tclsh
# This file is part of TopoTools, a VMD package to simplify 
# manipulating bonds other topology related properties.
#
# Copyright (c) 2009,2010,2011,2012 by Axel Kohlmeyer <akohlmey@gmail.com>
# $Id: topoatoms.tcl,v 1.14 2013/04/15 09:19:28 akohlmey Exp $

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

# guess missing atomic property from periodic table data. numbers are 
# taken from the corresponding lists in the molfile plugin header.
# TODO: additional guesses: element-name, mass-element, radius-element, ...
proc ::TopoTools::guessatomdata {sel what from} {
    set selstr [$sel text]

    # element names from PTE
    set elements {X H He Li Be B C N O F Ne Na Mg Al Si P
        S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As
        Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn
        Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho
        Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po
        At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md
        No Lr Rf Db Sg Bh Hs Mt Ds Rg}

    # element masses in AMU 
    set masses {0.00000 1.00794 4.00260 6.941 9.012182 10.811  
        12.0107 14.0067 15.9994 18.9984032 20.1797 22.989770
        24.3050 26.981538 28.0855 30.973761 32.065 35.453
        39.948 39.0983 40.078 44.955910 47.867 50.9415 
        51.9961 54.938049 55.845 58.9332 58.6934 63.546 
        65.409 69.723 72.64 74.92160  78.96 79.904 83.798
        85.4678 87.62 88.90585 91.224 92.90638 95.94 98.0
        101.07 102.90550 106.42 107.8682 112.411 114.818
        118.710 121.760 127.60 126.90447 131.293 132.90545
        137.327 138.9055 140.116 140.90765 144.24 145.0 
        150.36 151.964 157.25 158.92534 162.500 164.93032 
        167.259 168.93421 173.04 174.967 178.49 180.9479
        183.84 186.207 190.23 192.217 195.078 196.96655 
        200.59 204.3833 207.2 208.98038 209.0 210.0 222.0 
        223.0 226.0 227.0 232.0381 231.03588 238.02891 
        237.0 244.0 243.0 247.0 247.0 251.0 252.0 257.0
        258.0 259.0 262.0 261.0 262.0 266.0 264.0 269.0
        268.0 271.0 272.0}

    # VdW radii, ionic radii for elements that are commonly
    # ionic in typical systems. unknown elements set to 2.0.
    set radii {1.5 1.2 1.4 1.82 2.0 2.0 1.7 1.55 1.52 
        1.47 1.54 1.36 1.18 2.0 2.1 1.8 1.8 2.27 1.88 1.76
        1.37 2.0 2.0 2.0 2.0 2.0 2.0 2.0 1.63 1.4 1.39 1.07
        2.0 1.85 1.9 1.85 2.02 2.0 2.0 2.0 2.0 2.0 2.0 2.0
        2.0 2.0 1.63 1.72 1.58 1.93 2.17 2.0 2.06 1.98 2.16
        2.1 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0
        2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 1.72 1.66
        1.55 1.96 2.02 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0
        1.86 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 
        2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0}

    switch -- "$what-$from" {
        lammps-data {
            # shortcut for lammps data files
            guessatomdata $sel element mass
            guessatomdata $sel name element
            guessatomdata $sel radius element
        }

        element-mass {
            foreach a [lsort -real -unique [$sel get mass]] {
                set s [atomselect [$sel molid] "mass $a and ( $selstr )"]
                set idx 0
                foreach m $masses {
                    # this catches most cases. 
                    # we check the few exceptions later.
                    if {[expr abs($a-$m)] < 0.65} {
                        set idx [lsearch $masses $m]
                    }
                    # this is a hydrogen or deuterium and we flag it as hydrogen.
                    if {($a > 0.0 && $a < 2.2)} {
                        set idx 1
                    }
                    # Differentiate between Bismutium and Polonium. 
                    # The normal search will detect Polonium.
                    if {($a > 207.85 && $a < 208.99)} {
                        set idx 83
                    }
                    # Differentiate between Cobalt and Nickel
                    # The normal search will detect Nickel.
                    if {($a > 56.50 && $a < 58.8133)} {
                        set idx 27
                    }
                }
                $s set element [lindex $elements $idx]
                $s delete
            }
        }

        element-name {
            foreach n [lsort -ascii -unique [$sel get name]] {
                set s [atomselect [$sel molid] "name '$n' and ( $selstr )"]
                set idx [lsearch -nocase $elements $n]
                if { $idx < 0} {
                    set n [string range $n 0 1]
                    set idx [lsearch -nocase $elements $n]
                    if {$idx < 0} {
                        set n [string range $n 0 0]
                        set idx [lsearch -nocase $elements $n]
                        if {$idx < 0} {
                            set n X
                        } else {
                            set n [lindex $elements $idx]
                        }
                    } else {
                        set n [lindex $elements $idx]
                    }
                } else {
                    set n [lindex $elements $idx]
                }
                $s set element $n
                $s delete
            }
        }

        element-type {
            foreach t [lsort -ascii -unique [$sel get type]] {
                set s [atomselect [$sel molid] "type '$t' and ( $selstr )"]
                set idx [lsearch -nocase $elements $t]
                if { $idx < 0} {
                    set t [string range $t 0 1]
                    set idx [lsearch -nocase $elements $t]
                    if {$idx < 0} {
                        set t [string range $t 0 0]
                        set idx [lsearch -nocase $elements $t]
                        if {$idx < 0} {
                            set t X
                        } else {
                            set t [lindex $elements $idx]
                        }
                    } else {
                        set t [lindex $elements $idx]
                    }
                } else {
                    set t [lindex $elements $idx]
                }
                $s set element $t
                $s delete
            }
        }

        mass-element {
            foreach e [lsort -ascii -unique [$sel get element]] {
                set s [atomselect [$sel molid] "element '$e' and ( $selstr )"]
                set idx [lsearch -nocase $elements $e]
                set m 0.0
                if {$idx >= 0} {
                    set m [lindex $masses $idx]
                }
                $s set mass $m
                $s delete
            }   
        }

        name-element {
            # name is the same as element, only we go all uppercase.
            foreach e [lsort -ascii -unique [$sel get element]] {
                set s [atomselect [$sel molid] "element '$e' and ( $selstr )"]
                $s set name [string toupper $e]
                $s delete
            }
        }

        name-type {
            $sel set name [$sel get type]
        }

        radius-element {
            foreach e [lsort -ascii -unique [$sel get element]] {
                set s [atomselect [$sel molid] "element '$e' and ( $selstr )"]
                set idx [lsearch $elements $e]
                set r 2.0
                if {$idx >= 0} {
                    set r [lindex $radii $idx]
                }
                $s set radius $r
                $s delete
            }
        }

        type-element {
            # type is the same as element, only we go all uppercase.
            foreach e [lsort -ascii -unique [$sel get element]] {
                set s [atomselect [$sel molid] "element '$e' and ( $selstr )"]
                $s set type [string toupper $e]
                $s delete
            }
        }

        type-name {
            $sel set type [$sel get name]
        }

        default {
            vmdcon -err "guessatomdata: guessing '$what' from '$from' not implemented."
            vmdcon -err "Available are: element<-mass, element<-name, mass<element "
            vmdcon -err "name<element, radius<element name<type, type<element, type<name."
            return
        }
    }
}

