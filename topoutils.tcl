#!/usr/bin/tclsh
# TopoTools, a VMD package to simplify manipulating bonds 
# other topology related properties in VMD.
#
# Copyright (c) 2009,2010,2011,2012,2013 by Axel Kohlmeyer <akohlmey@gmail.com>
# support for crossterms contributed by Josh Vermass <vermass2@illinois.edu>
#
# $Id: topoutils.tcl,v 1.13 2013/09/08 13:10:05 akohlmey Exp $

# utility commands

# merge molecules from a list of molecule ids
# to form one new "molecule", i.e. system.
proc ::TopoTools::mergemols {mids} {

    # compute total number of atoms and collect
    # offsets and number of atoms of each piece.
    set ntotal 0
    set offset {}
    set numlist {}
    foreach m $mids {
        if {[catch {molinfo $m get numatoms} natoms]} {
            vmdcon -err "molecule id $m does not exist."
            return -1
        } else {
            # record number of atoms and offsets for later use.
            lappend offset $ntotal
            lappend numlist $natoms
            incr ntotal $natoms
        }
    }

    if {!$ntotal} {
        vmdcon -err "mergemols: combined molecule has no atoms."
        return -1
    }

    # create new molecule to hold data.
    set mol -1
    if {[catch {mol new atoms $ntotal} mol]} {
        vmdcon -err "mergemols: could not create new molecule: $mol"
        return -1
    } else {
        animate dup $mol
    }
    mol rename $mol [string range mergedmol-[join $mids -] 0 50]

    # copy data over piece by piece
    set bondlist {}
    set anglelist {}
    set dihedrallist {}
    set improperlist {}
    set ctermlist {}
    foreach m $mids off $offset num $numlist {
        set oldsel [atomselect $m all]
        set newsel [atomselect $mol "index $off to [expr {$off+$num-1}]"]

        # per atom props
        set cpylist {name type mass charge radius element x y z \
                         resname resid chain segname}
        $newsel set $cpylist [$oldsel get $cpylist]

        # assign structure data. we need to renumber indices
        set list [topo getbondlist both -molid $m]
        foreach l $list {
            lassign $l a b t o
            lappend bondlist [list [expr {$a+$off}] [expr {$b+$off}] $t $o]
        }

        set list [topo getanglelist -molid $m]
        foreach l $list {
            lassign $l t a b c 
            lappend anglelist [list $t [expr {$a+$off}] [expr {$b+$off}] [expr {$c+$off}]]
        }

        set list [topo getdihedrallist -molid $m]
        foreach l $list {
            lassign $l t a b c d
            lappend dihedrallist [list $t [expr {$a+$off}] [expr {$b+$off}] \
                                    [expr {$c+$off}] [expr {$d+$off}]]
        }
        set list [topo getimproperlist -molid $m]
        foreach l $list {
            lassign $l t a b c d
            lappend improperlist [list $t [expr {$a + $off}] [expr {$b + $off}] \
                                    [expr {$c + $off}] [expr {$d + $off}]]
        }
        set list [topo getcrosstermlist -molid $m]
        foreach l $list {
        	lassign $l a b c d e f g h
            lappend ctermlist [list [expr {$a + $off}] [expr {$b + $off}] \
                                    [expr {$c + $off}] [expr {$d + $off}] \
                                    [expr {$e + $off}] [expr {$f + $off}] \
                                    [expr {$g + $off}] [expr {$h + $off}]]
        }
        $oldsel delete
        $newsel delete
    }

    # apply structure info
    topo setbondlist both -molid $mol $bondlist
    topo setanglelist -molid $mol $anglelist
    topo setdihedrallist -molid $mol $dihedrallist
    topo setimproperlist -molid $mol $improperlist
	topo setcrosstermlist -molid $mol $ctermlist
    # set box to be largest of the available boxes
    set amax 0.0
    set bmax 0.0
    set cmax 0.0
    foreach m $mids {
        lassign [molinfo $m get {a b c}] a b c
        if {$a > $amax} {set amax $a}
        if {$b > $bmax} {set bmax $b}
        if {$c > $cmax} {set cmax $c}
    }
    molinfo $mol set {a b c} [list $amax $bmax $cmax]

    variable newaddsrep
    mol reanalyze $mol
    if {$newaddsrep} {
        adddefaultrep $mol
    }
    return $mol
}

# build a new molecule from one or more selections
proc ::TopoTools::selections2mol {sellist} {

    # compute total number of atoms and collect
    # offsets and number of atoms of each piece.
    set ntotal 0
    set offset {}
    set numlist {}
    foreach s $sellist {
        if {[catch {$s num} natoms]} {
            vmdcon -err "selection access error: $natoms"
            return -1
        } else {
            # record number of atoms and offsets for later use.
            lappend offset $ntotal
            lappend numlist $natoms
            incr ntotal $natoms
        }
    }

    if {!$ntotal} {
        vmdcon -err "selections2mol: combined molecule has no atoms."
        return -1
    }

    # create new molecule to hold data.
    set mol -1
    if {[catch {mol new atoms $ntotal} mol]} {
        vmdcon -err "selection2mol: could not create new molecule: $mol"
        return -1
    } else {
        animate dup $mol
    }
    mol rename $mol selections2mol-[molinfo num]

    # copy data over piece by piece
    set bondlist {}
    set anglelist {}
    set dihedrallist {}
    set improperlist {}
    set ctermlist {}
    foreach sel $sellist off $offset num $numlist {
        set newsel [atomselect $mol "index $off to [expr {$off+$num-1}]"]

        # per atom props
        set cpylist {name type mass charge radius element x y z \
                         resname resid chain segname}
        $newsel set $cpylist [$sel get $cpylist]

        # get atom index map for this selection
        set atomidmap [$sel get index]

        # assign structure data. we need to renumber indices
        set list [topo getbondlist both -sel $sel]
        foreach l $list {
            lassign $l a b t o
            set anew [expr [lsearch -sorted -integer $atomidmap $a] + $off]
            set bnew [expr [lsearch -sorted -integer $atomidmap $b] + $off]
            lappend bondlist [list $anew $bnew $t $o]
        }

        set list [topo getanglelist -sel $sel]
        foreach l $list {
            lassign $l t a b c 
            set anew [expr [lsearch -sorted -integer $atomidmap $a] + $off]
            set bnew [expr [lsearch -sorted -integer $atomidmap $b] + $off]
            set cnew [expr [lsearch -sorted -integer $atomidmap $c] + $off]
            lappend anglelist [list $t $anew $bnew $cnew]
        }

        set list [topo getdihedrallist -sel $sel]
        foreach l $list {
            lassign $l t a b c d
            set anew [expr [lsearch -sorted -integer $atomidmap $a] + $off]
            set bnew [expr [lsearch -sorted -integer $atomidmap $b] + $off]
            set cnew [expr [lsearch -sorted -integer $atomidmap $c] + $off]
            set dnew [expr [lsearch -sorted -integer $atomidmap $d] + $off]
            lappend dihedrallist [list $t  $anew $bnew $cnew $dnew]
        }
        set list [topo getimproperlist -sel $sel]
        foreach l $list {
            lassign $l t a b c d
            set anew [expr [lsearch -sorted -integer $atomidmap $a] + $off]
            set bnew [expr [lsearch -sorted -integer $atomidmap $b] + $off]
            set cnew [expr [lsearch -sorted -integer $atomidmap $c] + $off]
            set dnew [expr [lsearch -sorted -integer $atomidmap $d] + $off]
            lappend improperlist [list $t  $anew $bnew $cnew $dnew]
        }
        
        set list [topo getcrosstermlist -sel $sel]
        foreach l $list {
        	lassign $l a b c d e f g h
        	set anew [expr [lsearch -sorted -integer $atomidmap $a] + $off]
            set bnew [expr [lsearch -sorted -integer $atomidmap $b] + $off]
            set cnew [expr [lsearch -sorted -integer $atomidmap $c] + $off]
            set dnew [expr [lsearch -sorted -integer $atomidmap $d] + $off]
            set enew [expr [lsearch -sorted -integer $atomidmap $e] + $off]
            set fnew [expr [lsearch -sorted -integer $atomidmap $f] + $off]
            set gnew [expr [lsearch -sorted -integer $atomidmap $g] + $off]
            set hnew [expr [lsearch -sorted -integer $atomidmap $h] + $off]
        	lappend ctermlist [list $anew $bnew $cnew $dnew $enew $fnew $gnew $hnew]
        }
        $newsel delete
    }

    # apply structure info
    topo setbondlist both -molid $mol $bondlist
    topo setanglelist -molid $mol $anglelist
    topo setdihedrallist -molid $mol $dihedrallist
    topo setimproperlist -molid $mol $improperlist
    topo setcrosstermlist -molid $mol $ctermlist
    # set box to be largest of the available boxes
    set amax 0.0
    set bmax 0.0
    set cmax 0.0
    foreach sel $sellist {
        lassign [molinfo [$sel molid] get {a b c}] a b c
        if {$a > $amax} {set amax $a}
        if {$b > $bmax} {set bmax $b}
        if {$c > $cmax} {set cmax $c}
    }
    molinfo $mol set {a b c} [list $amax $bmax $cmax]

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
        vmdcon -err "replicatemol: no or bad nx/ny/nz replications given."
        return -1
    }
    set ntotal 0
    set natoms 0
    if {[catch {molinfo $mol get numatoms} natoms]} {
        vmdcon -err "replicatemol: molecule id $mol does not exist."
        return -1
    } else {
        set ntotal [expr {$natoms * $nrepl}]
    }
    if {!$natoms} {
        vmdcon -err "replicatemol: cannot replicate an empty molecule."
        return -1
    }

    set molname replicatedmol-$nrepl-x-$mol
    set newmol -1
    if {[catch {mol new atoms $ntotal} newmol]} {
        vmdcon -err "replicatemol: could not create new molecule: $mol"
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
	set ctermlist {}
	
    set oldsel [atomselect $mol all]
    set obndlist [topo getbondlist both -molid $mol]
    set oanglist [topo getanglelist -molid $mol]
    set odihlist [topo getdihedrallist -molid $mol]
    set oimplist [topo getimproperlist -molid $mol]
    set octermlist [topo getcrosstermlist -molid $mol]

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
		foreach l $octermlist {
            lassign $l a b c d e f g h
			lappend ctermlist [list [expr {$a + $ntotal}] [expr {$b + $ntotal}] \
                                    [expr {$c + $ntotal}] [expr {$d + $ntotal}] \
                                    [expr {$e + $ntotal}] [expr {$f + $ntotal}] \
                                    [expr {$g + $ntotal}] [expr {$h + $ntotal}]]
		}
        incr ntotal [$oldsel num]            
        $newsel delete
    }
    # apply structure info
    topo setbondlist both -molid $newmol $bondlist
    topo setanglelist -molid $newmol $anglelist
    topo setdihedrallist -molid $newmol $dihedrallist
    topo setimproperlist -molid $newmol $improperlist
    topo setcrosstermlist -molid $mol $ctermlist
    
    variable newaddsrep
    mol reanalyze $newmol
    if {$newaddsrep} {
        adddefaultrep $newmol
    }

    $oldsel delete
    return $newmol
}

