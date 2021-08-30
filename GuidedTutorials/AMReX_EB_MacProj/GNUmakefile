# ------------------------------------------------------------------------------
# Makefile to build examples
# ------------------------------------------------------------------------------
# Note the following environment variables need to be set to in order to build:
# AMREX_INSTALL_DIR = path to AMReX installation
# MPICXX = mpicxx wrapper
#
# If any are unset, they assume the default values for compilation on Cooley
# ------------------------------------------------------------------------------

DIM = 3

# set default values for any unset environment variables
ifeq ($(AMREX_INSTALL_DIR),)
  # AMREX_INSTALL_DIR = /projects/ATPESC2019/FASTMath/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/amrex-develop-fvkknjb363aqxpt6zu2tt6oojixnlnex
  AMREX_INSTALL_DIR = /projects/ATPESC2019/MathPackagesTraining/amrex/libamrex3D
endif

# These will only apply if they are not set as environment variables
# These are the paths to the AMReX and HYPRE libraries on Cooley for ATPESC
AMREX_INSTALL_DIR ?= /projects/ATPESC2019/MathPackagesTraining/amrex/libamrex3D
HYPRE_DIR         ?= /projects/ATPESC2019/MathPackagesTraining/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/hypre-develop-bjrvkgdwbwewvxhrv6wu4g6tgejcfnm4/

# Path to ParaView 5.6.1 MESA on Cooley
PARAVIEW_DIR ?= /projects/ATPESC2019/MathPackagesTraining/ParaView-5.6.1-osmesa-MPI-Linux-64bit

ifeq ($(MPICXX),)
  MPICXX = mpicxx
endif

FC = gfortran

CPPFLAGS = -Ishared -I$(AMREX_INSTALL_DIR)/include -I. -I$(HYPRE_DIR)/include
CPPFLAGS = -I$(AMREX_INSTALL_DIR)/include -I. -I$(HYPRE_DIR)/include
CXXFLAGS = -O2 -std=c++14
FFLAGS = -O2
LDFLAGS = -L$(AMREX_INSTALL_DIR)/lib -L$(HYPRE_DIR)/lib -Wl,-rpath,$(AMREX_INSTALL_DIR)/lib -Wl,-rpath,$(HYPRE_DIR)/lib

LIBRARIES = -lamrex -lHYPRE
LIBRARIES += -lgfortran
LIBRARIES += -lpthread


default: main$(DIM)d.exe

main$(DIM)d.exe: main.o MyParticleContainer.o 
	$(MPICXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

main.o: main.cpp MyParticleContainer.H 
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

MyParticleContainer.o: MyParticleContainer.cpp MyParticleContainer.H
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

.PHONY: clean realclean pltclean

clean:
	$(RM) *.o

realclean: clean
	$(RM) *~ *.exe

pltclean:
	$(RM) -rf plt*/ *.png eb*vtp

movie:
	$(PARAVIEW_DIR)/bin/pvpython paraview_vis_script.py
	ffmpeg -r 20 -i off_to_the_races.%04d.jpeg -s 1400x800 -vcodec libx264 -crf 30 -pix_fmt yuv420p off_to_the_races.mp4
