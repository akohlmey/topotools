#!/usr/bin/tclsh
# TopoTools, a VMD package to simplify manipulating bonds 
# other topology related properties in VMD.

# utility commands

# merge molecules from a list of molecule ids
# to form one new "molecule", i.e. system.
proc ::TopoTools::mergemols {mids} {

    # compute total number of atoms.
    set ntotal 0
    foreach m $mids {
        if {[catch {molinfo $m get numatoms} natoms]} {
            vmdcon -error "molecule id $m does not exist."
            return -1
        } else {
            incr ntotal $natoms
        }
    }

    if {!$ntotal} {
        vmdcon -error "mergemols: combined molecule has no atoms."
        return -1
    }

    # NOTE: we clamp length of the name to avoid a buffer
    # overflow in older VMD versions.
    set newmol [string range mergedmol-[join $mids -] 0 50]
    set mol -1
    if {[catch {mol new atoms $ntotal} mol]} {
        vmdcon -error "mergemols: could not create new molecule: $mol"
        return -1
    } else {
        animate dup $mol
    }
    mol rename $mol $newmol

    # copy data over piece by piece
    set ntotal 0
    set bondlist {}
    set anglelist {}
    set dihedrallist {}
    set improperlist {}
    foreach m $mids {
        set oldsel [atomselect $m all]
        set newsel [atomselect $mol \
                        "index $ntotal to [expr $ntotal + [$oldsel num] - 1]"]

        # per atom props
        set cpylist {name type mass charge radius element x y z \
                         resname resid chain segname}
        $newsel set $cpylist [$oldsel get $cpylist]

        # assign structure data. we need to renumber indices
        set list [topo getbondlist both -molid $m]
        foreach l $list {
            lassign $l a b t o
            lappend bondlist [list [expr {$a+$ntotal}] [expr {$b+$ntotal}] $t $o]
        }

        set list [topo getanglelist -molid $m]
        foreach l $list {
            lassign $l t a b c 
            lappend anglelist [list $t [expr {$a + $ntotal}] [expr {$b + $ntotal}] \
                                    [expr {$c + $ntotal}]]
        }

        set list [topo getdihedrallist -molid $m]
        foreach l $list {
            lassign $l t a b c d
            lappend dihedrallist [list $t [expr {$a + $ntotal}] [expr {$b + $ntotal}] \
                                    [expr {$c + $ntotal}] [expr {$d + $ntotal}]]
        }
        set list [topo getimproperlist -molid $m]
        foreach l $list {
            lassign $l t a b c d
            lappend improperlist [list $t [expr {$a + $ntotal}] [expr {$b + $ntotal}] \
                                    [expr {$c + $ntotal}] [expr {$d + $ntotal}]]
        }

        incr ntotal [$oldsel num]            
        $oldsel delete
        $newsel delete
    }

    # apply structure info
    topo setbondlist both -molid $mol $bondlist
    topo setanglelist -molid $mol $anglelist
    topo setdihedrallist -molid $mol $dihedrallist
    topo setimproperlist -molid $mol $improperlist
        
    variable newaddsrep
    mol reanalyze $mol
    if {$newaddsrep} {
        adddefaultrep $mol
    }

    return $mol
}


# create a larger system by replicating the original unitcell
# arguments: molecule id of molecule to replicate
#            multiples of the cell vectors defaulting to 1
# notes:     the cell is assumed to be orthorhombic.
#
proc ::TopoTools::replicatemol {mol nx ny nz} {

    if {[string equal $mol top]} {
        set mol [molinfo top]
    }

    # build translation vectors
    set xs [expr {-($nx-1)*0.5}]
    set ys [expr {-($ny-1)*0.5}]
    set zs [expr {-($nz-1)*0.5}]
    set transvecs {}
    for {set i 0} {$i < $nx} {incr i} {
        for {set j 0} {$j < $ny} {incr j} {
            for {set k 0} {$k < $nz} {incr k} {
                lappend transvecs [list [expr {$xs + $i}] [expr {$ys + $j}] [expr {$zs + $k}]]
            }
        }
    }

    # compute total number of atoms.
    set nrepl  [llength $transvecs]
    if {!$nrepl} {
        vmdcon -error "replicatemol: no or bad nx/ny/nz replications given."
        return -1
    }
    set ntotal 0
    set natoms 0
    if {[catch {molinfo $mol get numatoms} natoms]} {
        vmdcon -error "replicatemol: molecule id $mol does not exist."
        return -1
    } else {
        set ntotal [expr {$natoms * $nrepl}]
    }
    if {!$natoms} {
        vmdcon -error "replicatemol: cannot replicate an empty molecule."
        return -1
    }

    set molname replicatedmol-$nrepl-x-$mol
    set newmol -1
    if {[catch {mol new atoms $ntotal} newmol]} {
        vmdcon -error "replicatemol: could not create new molecule: $mol"
        return -1
    } else {
        animate dup $newmol
    }
    mol rename $newmol $molname

    # copy data over piece by piece
    set ntotal 0
    set bondlist {}
    set anglelist {}
    set dihedrallist {}
    set improperlist {}

    set oldsel [atomselect $mol all]
    set obndlist [topo getbondlist both -molid $mol]
    set oanglist [topo getanglelist -molid $mol]
    set odihlist [topo getdihedrallist -molid $mol]
    set oimplist [topo getimproperlist -molid $mol]

    set box [molinfo $mol get {a b c}]
    molinfo $newmol set {a b c} [vecmul $box [list $nx $ny $nz]]

    foreach v $transvecs {
        set newsel [atomselect $newmol \
                        "index $ntotal to [expr $ntotal + [$oldsel num] - 1]"]

        # per atom props
        set cpylist {name type mass charge radius element x y z \
                         resname resid chain segname}
        $newsel set $cpylist [$oldsel get $cpylist]

        set movevec {0.0 0.0 0.0}
        if {[catch {vecmul $v $box} movevec]} {
            vmdcon -warn "failure to compute translation vector from $v: $movevec. skipping..."
            continue
        }
        $newsel moveby $movevec
        # assign structure data. we need to renumber indices
        foreach l $obndlist {
            lassign $l a b t o
            lappend bondlist [list [expr {$a+$ntotal}] [expr {$b+$ntotal}] $t $o]
        }

        foreach l $oanglist {
            lassign $l t a b c 
            lappend anglelist [list $t [expr {$a + $ntotal}] [expr {$b + $ntotal}] \
                                    [expr {$c + $ntotal}]]
        }

        foreach l $odihlist {
            lassign $l t a b c d
            lappend dihedrallist [list $t [expr {$a + $ntotal}] [expr {$b + $ntotal}] \
                                    [expr {$c + $ntotal}] [expr {$d + $ntotal}]]
        }
        foreach l $oimplist {
            lassign $l t a b c d
            lappend improperlist [list $t [expr {$a + $ntotal}] [expr {$b + $ntotal}] \
                                    [expr {$c + $ntotal}] [expr {$d + $ntotal}]]
        }

        incr ntotal [$oldsel num]            
        $newsel delete
    }
    # apply structure info
    topo setbondlist both -molid $newmol $bondlist
    topo setanglelist -molid $newmol $anglelist
    topo setdihedrallist -molid $newmol $dihedrallist
    topo setimproperlist -molid $newmol $improperlist
    
    variable newaddsrep
    mol reanalyze $newmol
    if {$newaddsrep} {
        adddefaultrep $mol
    }

    return $newmol
}

