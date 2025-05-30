.SILENT:

VMFILES = pkgIndex.tcl topotools.tcl topoatoms.tcl \
        topobonds.tcl topoangles.tcl topodihedrals.tcl topoimpropers.tcl \
        topocrossterms.tcl topolammps.tcl topoutils.tcl topohelpers.tcl \
        topogromacs.tcl topovarxyz.tcl

VMVERSION = 1.10
DIR = $(PLUGINDIR)/noarch/tcl/topotools$(VMVERSION)

bins:
win32bins:
dynlibs:
staticlibs:
win32staticlibs:

distrib:
	@echo "Copying topotools $(VMVERSION) files to $(DIR)"
	mkdir -p $(DIR)
	cp $(VMFILES) $(DIR)
	cp README.md $(DIR)/README


