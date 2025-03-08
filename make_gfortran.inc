###############################################################################
##########  GNU Fortran Compiler  #############################################

FCHOME = /usr
FC = $(FCHOME)/bin/gfortran
COMPINC = -I./ -I$(FCHOME)/include -I$(FCHOME)/lib64/gfortran/modules -I$(INSTALLDIR)/include -I$(FCHOME)/lib64/gfortran/modules -I$(FCHOME)/lib/x86_64-linux-gnu/fortran/gfortran-mod-15
COMPLIBS = -L./ -L$(FCHOME)/lib64 -L${INSTALLDIR}/lib -L$(FCHOME)/lib/x86_64-linux-gnu

LIBS = $(COMPLIBS) $(COMPINC)

# Debugging flags
ifeq ($(RUN), DEBUG)
    # consider including -lefence for debugging
    #  Note that -Wall   = -Waliasing -Wampersand -Wconversion -Wsurprising -Wc-binding-type
    #                      -Wintrinsics-std -Wtabs -Wintrinsic-shadow -Wline-truncation
    #                      -Wtarget-lifetime -Winteger-division -Wreal-q-constant -Wunused
    #                      -Wundefined-do-loop
    #    and     -Wextra = -Wcompare-reals, -Wunused-parameter -Wdo-subscript
    FFLAGS = -O0 -g3 -Wall -Wextra -Wcharacter-truncation -Wimplicit-interface -Wunderflow -Wuninitialized -Wmaybe-uninitialized -Wno-surprising -fimplicit-none -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace -ffpe-trap=invalid,zero,overflow -ffpe-summary=all -finit-integer=-2147483647 -finit-real=snan -fno-unsafe-math-optimizations -fsignaling-nans -frounding-math
endif

# Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g -pg -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math
endif

# Production run flags
ifeq ($(RUN), OPT)
    FFLAGS = -O3 -w -fno-math-errno -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -ftree-vectorize -funroll-loops -finline-functions 
endif

  # Preprocessing flags
FPPFLAGS = -x f95-cpp-input $(ncFPPFLAG) $(grbFPPFLAG) $(memFPPFLAG)
  # Extra flags
EXFLAGS =

###############################################################################
