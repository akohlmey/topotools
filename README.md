# VMD TopoTools package. Version 1.10

Copyright (c) 2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2022,2023,2024,2025
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

The public git repository is at https://github.com/akohlmey/topotools

TopoTools version 1.10 is the **final** release created by **me**.
I have no more plans to further develop and maintain this package.
It is therefore available for "adoption". Please contact me via email
or PM if you want to take over.

## Installation

TopoTools is written entirely in the Tcl scripting language
for use as a plugin with VMD. A version of TopoTools is already
bundled with VMD, to update it with the newer, downloaded version
unpack the TopoTools archive, which will create a directory
containing that various Tcl script files. If the directory is not
named topotools1.10, please rename it accordingly. Now, find the
plugin folder in your VMD installation where it already has
a topotools1.x directory and move the folder with the new version
next to it.  If it already has the topotools1.10 folder, overwrite
the files inside with the new version. VMD should use the new
version automatically at the next start.

