.SILENT:

VMFILES = pkgIndex.tcl README topotools.tcl \
	topobonds.tcl topoangles.tcl topodihedrals.tcl topoimpropers.tcl \
	topolammps.tcl

VMVERSION = 1.0
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


