#!/usr/bin/tclsh
# This file is part of TopoTools, a VMD package to simplify 
# manipulating bonds and other topology related properties.
#
# Copyright (c) 2009,2010,2011 by Axel Kohlmeyer <akohlmey@gmail.com>
# $Id: topogromacs.tcl,v 1.3 2011/02/02 21:33:29 akohlmey Exp $

# high level subroutines for supporting gromacs topology files.
#
# write a fake gromacs topology format file that can be used in combination
# with a .gro/.pdb coordinate file for generating .tpr files needed to use
# Some of the more advanced gromacs analysis tools for simulation data that
# was not generated with gromacs.
#
# IMPORTANT NOTE: this script differs from other topotools script in that
# it does not check whether fragments are fully contained in the selection.
# it will output a topology with exactly the same number of atoms as the
# selection has. in case of partially contained fragments, new molecule types
# will be created.
#
# Arguments:
# filename = name of topology file
# mol = molecule
# sel = selection
proc ::TopoTools::writegmxtop {filename mol sel {flags none}} {

    if {[catch {open $filename w} fp]} {
        vmdcon -error "writegmxtop: problem opening gromacs topology file: $fp\n"
        return -1
    }

    vmdcon -info "Generating a 'faked' gromacs topology file: $filename"
    # get a list of fragments, i.e. individual molecules
    set fragmap [lsort -integer -unique [$sel get fragment]]
    set typemap [lsort -ascii -unique [$sel get type]]
    set selstr [$sel text]

    puts $fp "; 'fake' gromacs topology generated from topotools."
    puts $fp "; WARNING| the purpose of this topology is to allow using the  |WARNING"
    puts $fp "; WARNING| analysis tools from gromacs for non gromacs data.   |WARNING"
    puts $fp "; WARNING| it cannot be used for a simulation.                 |WARNING"
    puts $fp "\n\[ defaults \]\n; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ"
    puts $fp "1 3 yes 0.5 0.5"
    puts $fp "\n\[ atomtypes \]\n; name bond_type mass charge ptype sigma epsilon"
    foreach t $typemap {
        if {[string is integer $t]} {
            puts $fp "type$t C 1.0 0.0 A 0.0 0.0"
        } else {
            puts $fp "$t C 1.0 0.0 A 0.0 0.0"
        }
    }
    puts $fp "\n\[ bondtypes \]\n; i j func b0 kb\n  C C 1 0.13 1000.0 ; totally bogus"
    ; # puts $fp "\n\[ constrainttypes \]\n;"
    puts $fp "\n\[ angletypes \]\n; i j k func th0 cth\n  C C C 1 109.500 100.0 ; totally bogus"
    puts $fp "\n\[ dihedraltypes \]\n; i j k l func coefficients\n  C C C C 1 0.0 3 10.0 ; totally bogus"
    set fraglist {}
    set fragcntr {}
    set nlold {}
    set tlold {}
    set count 0
    foreach frag $fragmap {
        set fsel [atomselect $mol "(fragment $frag) and ($selstr)"]
        set nlist [$fsel get name]
        set tlist [$fsel get type]
        if {[listcmp $nlist $nlold] || [listcmp $tlist $tlold]} {
            vmdcon -info "Found new moleculetype: fragment \#$frag natoms=[$fsel num]"
            display update ui
            if {[llength $fraglist] > [llength $fragcntr]} {
                lappend fragcntr $count
            }
            puts $fp ""
            set molname "molecule[llength $fragcntr]"
            lappend fraglist $molname
            set count 1
            set nlold $nlist
            set tlold $tlist
            puts $fp "\n\[ moleculetype \]"
            puts $fp "; Name      nrexcl\n$molname     3"
            puts $fp "\n\[ atoms \]"
            puts $fp "; nr  type  resnr residue atom cgnr charge  mass"
            set offs [lindex [$fsel get index] 0]
            set resoffs [lindex [$fsel get residue] 0]
            foreach nr [$fsel get serial] type [$fsel get type] \
                name [$fsel get name] residue [$fsel get residue] \
                resname [$fsel get resname] charge [$fsel get charge] \
                mass [$fsel get mass] {
                    # fix up some data that gromacs cannok grok
                    if {[string is integer $type]} {set type "type$type"}
                    if {[string is integer $resname]} {set resname "RES$resname"}
                    puts $fp [format "% 6d %11s % 6d %8s %6s % 6d %10.4f %10.4f"  \
                                  [expr {$nr - $offs}]  $type \
                                  [expr {$residue - $resoffs + 1}] $resname $name \
                                  [expr {$residue - $resoffs + 1}] $charge $mass ]
            }
            set list [bondinfo getbondlist $fsel none]
            if {[llength $list]} {
                puts $fp "\n\[ bonds \]\n; i  j  func"
                foreach b $list {
                    lassign $b i j
                    set i [expr {$i - $offs}]
                    set j [expr {$j - $offs}]
                    incr i; incr j
                    puts $fp "$i $j 1"
                }
            }
            set list [angleinfo getanglelist $fsel]
            if {[llength $list] > 0} {
                puts $fp "\n\[ angles \]\n; i  j  k  func"
                foreach b $list {
                    lassign $b t i j k
                    set i [expr {$i - $offs}]
                    set j [expr {$j - $offs}]
                    set k [expr {$k - $offs}]
                    incr i; incr j; incr k
                    puts $fp "$i $j $k 1"
                }
            }
            set list [dihedralinfo getdihedrallist $fsel]
            if {[llength $list] > 0} {
                puts $fp "\n\[ dihedrals \]\n; i  j  k  l  func"
                foreach b $list {
                    lassign $b t i j k l
                    set i [expr {$i - $offs}]
                    set j [expr {$j - $offs}]
                    set k [expr {$k - $offs}]
                    set l [expr {$l - $offs}]
                    incr i ; incr j; incr k ; incr l
                    puts $fp "$i $j $k $l 1"
                }
            }
        } else {
            incr count
        }
        $fsel delete
    }
    lappend fragcntr $count
    
    puts $fp "\n\[ system \]\n; Name\nvmdmolecule$mol\n"
    puts $fp "\n\[ molecules \]\n; Compound    \#mols"
    vmdcon -info "Found [llength $fraglist] moleculetypes."
    foreach name $fraglist num $fragcntr {
        vmdcon -info "$num x $name"
        puts $fp "$name    $num"
    }
    close $fp
    return
}

