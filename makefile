##############################################################################
#  Makefile for libmetreader.a
#
#    User-specified flags are in this top block
#
###############################################################################

#      This file is a component of the volcanic ash transport and dispersion model Ash3d,
#      written at the U.S. Geological Survey by Hans F. Schwaiger (hschwaiger@usgs.gov),
#      Larry G. Mastin (lgmastin@usgs.gov), and Roger P. Denlinger (roger@usgs.gov).

#      The model and its source code are products of the U.S. Federal Government and therefore
#      bear no copyright.  They may be copied, redistributed and freely incorporated 
#      into derivative products.  However as a matter of scientific courtesy we ask that
#      you credit the authors and cite published documentation of this model (below) when
#      publishing or distributing derivative products.

#      Schwaiger, H.F., Denlinger, R.P., and Mastin, L.G., 2012, Ash3d, a finite-
#         volume, conservative numerical model for ash transport and tephra deposition,
#         Journal of Geophysical Research, 117, B04204, doi:10.1029/2011JB008968. 

#      We make no guarantees, expressed or implied, as to the usefulness of the software
#      and its documentation for any purpose.  We assume no responsibility to provide
#      technical support to users of this software.

#      Sequence of commands:
#      "make"  compiles the libmetreader.a library
#      "make all" builds the library, the tools executables and copies to bin
#      "make check" runs the test cases in /tests
#      "make install" copies the contents of volcano-ash3d-metreader/bin to the install location
#                        e.g. /opt/USGS
#
#  SYSTEM specifies which compiler to use
#    Current available options are:
#      gfortran , ifort , aocc , nvhpc
#    This variable cannot be left blank
#
SYSTEM = gfortran
SYSINC = make_$(SYSTEM).inc
#
#  RUN specifies which collection of compilation flags that should be run
#    Current available options are:
#      DEBUG : includes debugging info and issues warnings
#      PROF  : includes profiling flags with some optimization
#      OPT   : includes optimizations flags for fastest runtime
#    This variable cannot be left blank

#RUN = DEBUG
#RUN = PROF
RUN = OPT

INSTALLDIR=/opt/USGS
#
# DATA FORMATS
#  For each data format you want to include in the library, set the corresponding
#  variable below to 'T'.  Set to 'F' any you do not want compiled or any unavailable
USENETCDF = T
USEGRIB   = T

# MEMORY
# If you need pointer arrays instead of allocatable arrays, set this to 'T'
USEPOINTERS = F

# DATA LOCATION
WINDROOT=/data/WindFiles
###############################################################################
#####  END OF USER SPECIFIED FLAGS  ###########################################
###############################################################################

FPPFLAGS = 
ifeq ($(USENETCDF), T)
 ncFPPFLAG = -DUSENETCDF
 ncOBJS = MetReader_NetCDF.o
 nclib = -lnetcdf -lnetcdff 
else
 ncFPPFLAG =
 ncOBJS =
 nclib =
endif
ifeq ($(USEGRIB), T)
 grbFPPFLAG = -DUSEGRIB
 grbOBJS = MetReader_GRIB.o MetReader_GRIB_index.o
 grblib = -leccodes -leccodes_f90
else
 grbFPPFLAG =
 grbOBJS =
 grblib =
endif

ifeq ($(USEPOINTERS), T)
 memFPPFLAG = -DUSEPOINTERS
else
 memFPPFLAG =
endif

# location of HoursSince and projection
USGSROOT = $(INSTALLDIR)
USGSLIBDIR = -L$(INSTALLDIR)/lib
USGSINC = -I$(INSTALLDIR)/include
USGSLIB = $(USGSINC) $(USGSLIBDIR) -lhourssince -lprojection

###############################################################################
###############################################################################

###############################################################################
# Import the compiler-specific include file.  Currently one of:
#  GNU Fortran Compiler
#  Intel Fortran Compiler
#  AMD Optimizing C/C++/Fortran Compiler (aocc)
#  Nvidia HPC Fortran Compiler (ncfortran)
include $(SYSINC)
###############################################################################

LIB = libMetReader.a

ifeq ($(USEGRIB), T)
  GRIBTOOL = bin/gen_GRIB_index
else
  GRIBTOOL =
endif

EXEC = \
 bin/MetRegrid \
 bin/MetSonde  \
 bin/MetTraj_F \
 bin/MetTraj_B \
 bin/MetCheck  \
 bin/MetProbe  \
 bin/MR_ASCII_check \
 bin/makegfsncml $(GRIBTOOL)

AUTOSCRIPTS = \
 autorun_scripts/autorun_gfs.sh \
 autorun_scripts/get_gfs.sh     \
 autorun_scripts/convert_gfs.sh \
 autorun_scripts/autorun_ecmwf.sh \
 autorun_scripts/get_ecmwf.sh     \
 autorun_scripts/convert_ecmwf.sh \
 autorun_scripts/autorun_nam.sh \
 autorun_scripts/get_nam.sh     \
 autorun_scripts/autorun_NCEP_50YearReanalysis.sh \
 autorun_scripts/get_NCEP_50YearReanalysis.sh     \
 autorun_scripts/grib2nc.sh \
 autorun_scripts/prune_windfiles.sh \
 autorun_scripts/get_gmao.sh \
 autorun_scripts/probe_volc.sh
###############################################################################
# TARGETS
###############################################################################

all: ${lib} tools
	$(info -------------------------------------)
	$(info  Before running 'make check', set the WINDROOT and USGSROOT environment variables.)
	$(info  e.g. execute: )
	$(info  export WINDROOT="${WINDROOT}")
	$(info  export USGSROOT="${USGSROOT}")
	$(info            or: source set_MR.env)
	$(info -------------------------------------)

lib: $(LIB)

tools: $(EXEC)

libMetReader.a: MetReader.F90 MetReader.o $(ncOBJS) $(grbOBJS) MetReader_Grids.o MetReader_ASCII.o makefile $(SYSINC)
	ar rcs libMetReader.a MetReader.o $(ncOBJS) $(grbOBJS) MetReader_Grids.o MetReader_ASCII.o

MetReader.o: MetReader.F90 makefile $(SYSINC)
	sh get_version.sh
	$(FC) $(FPPFLAGS) $(EXFLAGS) $(LIBS) $(USGSLIB) -c MetReader.F90
MetReader_Grids.o: MetReader_Grids.f90 MetReader.o makefile $(SYSINC)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) $(USGSLIB) -c MetReader_Grids.f90
MetReader_ASCII.o: MetReader_ASCII.f90 MetReader.o makefile $(SYSINC)
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(USGSLIB) -c MetReader_ASCII.f90

ifeq ($(USENETCDF), T)
MetReader_NetCDF.o: MetReader_NetCDF.F90 MetReader.o makefile $(SYSINC)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) $(nclib) $(USGSLIB) -c MetReader_NetCDF.F90
endif
ifeq ($(USEGRIB), T)
MetReader_GRIB_index.o: MetReader_GRIB_index.f90 makefile $(SYSINC)
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(grblib) $(USGSLIB) -c MetReader_GRIB_index.f90
MetReader_GRIB.o: MetReader_GRIB.f90 MetReader_GRIB_index.o MetReader.o makefile $(SYSINC)
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(grblib) -c MetReader_GRIB.f90
bin/gen_GRIB_index: gen_GRIB_index.f90 MetReader_GRIB_index.o makefile $(SYSINC) libMetReader.a
	mkdir -p bin
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(grblib) -c gen_GRIB_index.f90
	$(FC) $(FFLAGS) $(EXFLAGS) MetReader_GRIB_index.o gen_GRIB_index.o $(LIBS) $(grblib) -o bin/gen_GRIB_index
endif

bin/MetRegrid: tools/MetRegrid.f90 makefile $(SYSINC) libMetReader.a
	mkdir -p bin
	$(FC)  $(FFLAGS) $(EXFLAGS) tools/MetRegrid.f90 -o bin/MetRegrid -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB)

bin/MetSonde: tools/MetSonde.f90 makefile $(SYSINC) libMetReader.a
	mkdir -p bin
	$(FC) $(FFLAGS) $(EXFLAGS) -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB) -c tools/MetSonde.f90
	$(FC) $(FFLAGS) $(EXFLAGS) MetSonde.o  -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB) -o bin/MetSonde
bin/MetTraj_F: tools/MetTraj.F90 makefile $(SYSINC) libMetReader.a
	mkdir -p bin
	$(FC) $(FPPFLAGS) -DFORWARD  $(FFLAGS) $(EXFLAGS) tools/MetTraj.F90 -o bin/MetTraj_F -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB)
bin/MetTraj_B: tools/MetTraj.F90 makefile $(SYSINC) libMetReader.a
	mkdir -p bin
	$(FC) $(FPPFLAGS) -DBACKWARD $(FFLAGS) $(EXFLAGS) tools/MetTraj.F90 -o bin/MetTraj_B -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB)
bin/MetCheck: tools/MetCheck.f90 makefile $(SYSINC) libMetReader.a
	mkdir -p bin
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(nclib) $(grblib) -c tools/MetCheck.f90
	$(FC) $(FFLAGS) $(EXFLAGS) MetCheck.o -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB) -o bin/MetCheck
bin/MetProbe: tools/MetProbe.f90 makefile $(SYSINC) libMetReader.a
	mkdir -p bin
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(nclib) $(grblib) -c tools/MetProbe.f90
	$(FC) $(FFLAGS) $(EXFLAGS) MetProbe.o -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB) -o bin/MetProbe
bin/makegfsncml: tools/makegfsncml.f90 makefile $(SYSINC)
	mkdir -p bin
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(nclib) -c tools/makegfsncml.f90
	$(FC) $(FFLAGS) $(EXFLAGS) makegfsncml.o  $(LIBS) $(nclib) -o bin/makegfsncml
bin/MR_ASCII_check: tools/MR_ASCII_check.f90 makefile $(SYSINC)
	mkdir -p bin
	$(FC) $(FFLAGS) $(EXFLAGS) tools/MR_ASCII_check.f90 -o bin/MR_ASCII_check -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB)

check: $(EXEC)
	bash run_tests.sh

clean:
	rm -f *.o *__genmod.f90 *__genmod.mod
	rm -f *.mod
	rm -f lib*.a
	rm -f $(EXEC)

install: all
	install -d $(INSTALLDIR)/lib/
	install -d $(INSTALLDIR)/include/
	install -d $(INSTALLDIR)/bin/
	install -d $(INSTALLDIR)/bin/autorun_scripts
	install -d $(INSTALLDIR)/share
	install -m 644 $(LIB) $(INSTALLDIR)/lib/
	install -m 644 *.mod $(INSTALLDIR)/include/
	install -m 755 $(EXEC) $(INSTALLDIR)/bin/
	install -m 755 $(AUTOSCRIPTS) $(INSTALLDIR)/bin/autorun_scripts/
	install -m 755 tools/GMT_plot_traj.sh $(INSTALLDIR)/bin/
	install -m 644 share/volc_NOVAC.txt $(INSTALLDIR)/share/volc_NOVAC.txt

uninstall:
	rm -f $(INSTALLDIR)/lib/$(LIB)
	rm -f $(INSTALLDIR)/include/metreader.mod
	rm -f $(INSTALLDIR)/bin/MetRegrid
	rm -f $(INSTALLDIR)/bin/MetSonde
	rm -f $(INSTALLDIR)/bin/MetTraj_F
	rm -f $(INSTALLDIR)/bin/MetTraj_B
	rm -f $(INSTALLDIR)/bin/MetCheck
	rm -f $(INSTALLDIR)/bin/makegfsncml
	rm -f $(INSTALLDIR)/bin/gen_GRIB_index
	rm -f $(INSTALLDIR)/bin/MetProbe
	rm -f $(INSTALLDIR)/bin/MR_ASCII_check
	rm -f $(INSTALLDIR)/bin/autorun_scripts/autorun_gfs.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/autorun_nam.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/autorun_ecmwf.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/autorun_NCEP_50YearReanalysis.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/get_gfs.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/get_nam.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/get_ecmwf.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/get_NCEP_50YearReanalysis.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/convert_gfs.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/convert_ecmwf.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/grib2nc.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/prune_windfiles.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/get_gmao.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/probe_volc.sh
	rm -f $(INSTALLDIR)/bin/GMT_plot_traj.sh
	rm -f $(INSTALLDIR)/share/volc_NOVAC.txt

