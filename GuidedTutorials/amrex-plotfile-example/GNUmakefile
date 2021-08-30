# ------------------------------------------------------------------------------
# Makefile to build examples
# ------------------------------------------------------------------------------
# Note the following environment variables need to be set to in order to build:
# AMREX_INSTALL_DIR = path to AMReX installation
# MPICXX = mpicxx wrapper
#
# If any are unset, they assume the default values for compilation on Cooley
# ------------------------------------------------------------------------------

# set default values for any unset environment variables
ifeq ($(AMREX_INSTALL_DIR),)
  AMREX_INSTALL_DIR = /projects/ATPESC2019/FASTMath/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/amrex-develop-fvkknjb363aqxpt6zu2tt6oojixnlnex
endif
ifeq ($(MPICXX),)
  MPICXX = mpicxx
endif

CPPFLAGS = -Ishared -I$(AMREX_INSTALL_DIR)/include
CXXFLAGS = -O2 -std=c++11
FFLAGS = -O2
LDFLAGS = -L$(AMREX_INSTALL_DIR)/lib -Wl,-rpath,$(AMREX_INSTALL_DIR)/lib

LIBRARIES = -lamrex

LIBRARIES += -lgfortran

main.exe: main.cpp
	$(MPICXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

.PHONY: movie clean realclean pltclean

movie:
	ls -1 plt*/Header | tee movie.visit

clean:
	$(RM) *.o

realclean: clean
	$(RM) *~ *.exe

pltclean:
	$(RM) -rf plt*/ *.png *_diagnostics.txt __pycache__
