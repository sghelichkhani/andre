#------------------------------------------------------------------------------
#                  Central Makefile for the program TERRA
#------------------------------------------------------------------------------

# Turn on/off using direct vtk export by setting this to yes/no
USE_VTKW = yes

# Set directory information
TARGETDIR = ..
ARCHDIR   = arch
VTKWDIR   = vtkw
SRCDIR    = code
OBJDIR    = objects

# Set name executable
BIN = $(TARGETDIR)/terra

# Set information for automatic generation of documentation
TERRA_DOXYCONF   = ../Docs/SourceCode/DoxyConfig/terra.conf
TERRA_DOXYDOCDIR = ../Docs/SourceCode/terra
VTKW_DOXYCONF    = ../Docs/SourceCode/DoxyConfig/vtkw.conf
VTKW_DOXYDOCDIR  = ../Docs/SourceCode/vtkw

# Include architecture/compiler dependent settings
ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
ifneq ($(findstring doc,$(MAKECMDGOALS)),doc)
include $(ARCHDIR)/make.$(ARCH)
endif
endif

# Generate file set for dependency info
OBJ = $(OBJDIR)/driver.o $(OBJDIR)/commun.o $(OBJDIR)/commpi.o \
      $(OBJDIR)/convct.o $(OBJDIR)/adjoin.o $(OBJDIR)/opinit.o $(OBJDIR)/opcalc.o \
      $(OBJDIR)/opcoar.o $(OBJDIR)/solver.o $(OBJDIR)/utlity.o \
      $(OBJDIR)/cvinit.o $(OBJDIR)/forces.o $(OBJDIR)/energy.o \
      $(OBJDIR)/viscos.o $(OBJDIR)/plates.o $(OBJDIR)/harmon.o \
      $(OBJDIR)/topog.o  $(OBJDIR)/linpck.o $(OBJDIR)/export.o

# Add objects from the VTKW part, if this functionality is desired
ifeq "$(USE_VTKW)" "yes"
OBJ :=   $(OBJ) $(OBJDIR)/vtkw_abort.o  $(OBJDIR)/vtkw_databox.o   \
         $(OBJDIR)/vtkw_encode.o        $(OBJDIR)/vtkw_interface.o \
         $(OBJDIR)/vtkw_io.o            $(OBJDIR)/vtkw_manager.o   \
         $(OBJDIR)/vtkw_mesh.o
endif

# Default target is generation executable
$(BIN) : $(OBJDIR) $(OBJ)
	$(LK) $(LFLAGS) -o $(BIN) $(OBJ) $(LIBS)

# Generic rule for built process
$(OBJDIR)/%.o : $(SRCDIR)/%.f $(SRCDIR)/size.h $(SRCDIR)/para.h
	$(FC) $(FFLAGS) -o $@ -c $<

# Rule for building VTKW part
$(OBJDIR)/%.o : $(VTKWDIR)/%.c
	$(CC) $(VTKW_CFLAGS) -o $@ -c $<

# Targets for housekeeping
.PHONY: clean realclean doc clean-doc

$(OBJDIR) : ;
	mkdir $(OBJDIR)

clean:
	$(RM) $(OBJDIR)/*.o

realclean: clean
	$(RM) $(BIN)
	$(RM) -r $(OBJDIR)

doc:	doc-terra

doc-terra:
	doxygen $(TERRA_DOXYCONF)
	@echo " "
	@echo "-----------------------------------------------------------"
	@echo " See "$(TERRA_DOXYDOCDIR)/index.html" for documentation!"
	@echo "-----------------------------------------------------------"

doc-vtkw:
	doxygen $(VTKW_DOXYCONF)
	@echo " "
	@echo "-----------------------------------------------------------"
	@echo " See "$(VTKW_DOXYDOCDIR)/index.html" for documentation!"
	@echo "-----------------------------------------------------------"

clean-doc:
	$(RM) -r $(TERRA_DOXYDOCDIR)
	$(RM) -r $(VTKW_DOXYDOCDIR)
