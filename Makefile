.SILENT:

VMFILES = topotools.tcl topobonds.tcl topoangles.tcl pkgIndex.tcl README
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


