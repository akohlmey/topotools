# VMD TopoTools package. Version 1.8

Copyright (c) 2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020
 by Axel Kohlmeyer <akohlmey@gmail.com>

[![DOI](https://zenodo.org/badge/13922095.svg)](https://zenodo.org/badge/latestdoi/13922095)

This package contains contributed features from:
- Josh Vermaas (fully working gromacs topology files for CHARMM)
- Konstantin W (replicatemols for non-orthogonal cells)

-------------------

## Overview

TopoTools is a plugin for [VMD](http://www.ks.uiuc.edu/Research/vmd/)
providing a collection of Tcl commands to manipulate, build, read
and write topologies (i.e. bonds, angles, dihedrals, etc.
and their corresponding properties (type, order, etc.).

## Updates

Updated and bugfix versions are also available from:
http://sites.google.com/site/akohlmey/software/topotools
and: https://github.com/akohlmey/topotools
That page also links to tutorials demonstrating the use
of the package for different purposes.

## TODO
  - improve "topo setbonds" to be more efficient for large sets
    of bonds to be added.
  - topo copybonds <fromsel> <tosel>
  - topo guessbonds <sel>  (use bondsrecalc and restore bonds
    outside the selection)
  - support for generating and replicating bonds across
    periodic boundaries (for molecular crystal models).
  - better documentation, more tutorials
  - more tools to read/write custom topology files
    (e.g. gromacs top?, amber parmtop?)
  - API to read/parse/store force field database information.

## Installation

TopoTools is written entirely in the Tcl scripting language
for use as a plugin with VMD. A version of TopoTools is already
bundled with VMD, to update it with the newer, downloaded version
unpack the TopoTools archive, which will create a directory
containing that various Tcl script files. If the directory is not
named topotools1.8, please rename it. Now, find the plugin folder
in your VMD installation where it already has a topotools1.x
directory and move the folder with the new version next to it.
If it already has the topotools1.8 folder, overwrite the files
inside with the new version. VMD should use the new version
automatically at the next start.

## Feedback

Please report any problems, bugs, suggestions, inquiries
or code contributions to the github hosted SCM page at:
https://github.com/akohlmey/topotools/issues

