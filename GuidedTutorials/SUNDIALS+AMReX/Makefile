# ------------------------------------------------------------------------------
# Makefile to build examples
# ------------------------------------------------------------------------------
# Note the following environment variables need to be set to in order to build:
# AMREX_INSTALL_DIR = path to AMReX installation
# SUNDIALS_INSTALL_DIR = path to SUNDIALS installation
# MPICXX = mpicxx wrapper
#
# If any are unset, they assume the default values for compilation on Cooley
# ------------------------------------------------------------------------------

# set default values for any unset environment variables
ifeq ($(AMREX_INSTALL_DIR),)
#  AMREX_INSTALL_DIR = /gpfs/mira-home/balay/spack/opt/spack/linux-rhel7-haswell/gcc-8.2.0/amrex-21.07-5zqi6hw5mzgyof6xddypdjj4tdxsbk75
  AMREX_INSTALL_DIR = /lus/grand/projects/ATPESC2021/usr/MathPackages/spack/opt/spack/linux-rhel7-haswell/gcc-8.2.0/amrex-21.07-iwaygptbziyyoogyo55ktavrf7vrfoqa
endif
ifeq ($(SUNDIALS_INSTALL_DIR),)
  #SUNDIALS_INSTALL_DIR = /gpfs/mira-home/balay/spack/opt/spack/linux-rhel7-haswell/gcc-8.2.0/sundials-5.7.0-65uxu6ysj4l5ytbfpdpidbvoc727qtnn
  SUNDIALS_INSTALL_DIR = /lus/grand/projects/ATPESC2021/usr/MathPackages/spack/opt/spack/linux-rhel7-haswell/gcc-8.2.0/sundials-5.7.0-5pkorbm47t33a3ljsu2dljseiyawscvu
endif
ifeq ($(MPICXX),)
  MPICXX = mpicxx
endif

CPPFLAGS = -Ishared -I$(AMREX_INSTALL_DIR)/include -I$(SUNDIALS_INSTALL_DIR)/include
CXXFLAGS = -O2 -std=c++14
FFLAGS = -O2
LDFLAGS = -L$(AMREX_INSTALL_DIR)/lib -L$(SUNDIALS_INSTALL_DIR)/lib -L$(SUNDIALS_INSTALL_DIR)/lib64 -Wl,-rpath,$(AMREX_INSTALL_DIR)/lib,-rpath,$(SUNDIALS_INSTALL_DIR)/lib,-rpath,$(SUNDIALS_INSTALL_DIR)/lib64

LIBRARIES = -lamrex -lsundials_cvode -lsundials_arkode -lpthread

LIBRARIES += -lgfortran

default: HandsOn1.exe HandsOn2.exe HandsOn3.exe

HandsOn1.exe: HandsOn1.o shared/NVector_Multifab.o shared/Utilities.o
	$(MPICXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

HandsOn1.o: HandsOn1.cpp HandsOn1.h
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

HandsOn2.exe: HandsOn2.o shared/NVector_Multifab.o shared/Utilities.o
	$(MPICXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

HandsOn2.o: HandsOn2.cpp HandsOn2.h
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

HandsOn3.exe: HandsOn3.o shared/NVector_Multifab.o shared/Utilities.o
	$(MPICXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

HandsOn3.o: HandsOn3.cpp HandsOn3.h
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

Advection-Diffusion.exe: Advection-Diffusion.o shared/NVector_Multifab.o shared/Utilities.o
	$(MPICXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

Advection-Diffusion.o: Advection-Diffusion.cpp Advection-Diffusion.h
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

shared/NVector_Multifab.o: shared/NVector_Multifab.cpp shared/NVector_Multifab.h
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

shared/Utilities.o: shared/Utilities.cpp shared/Utilities.h
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

.PHONY: movie clean realclean pltclean

movie:
	ls -1 plt*/Header | tee movie.visit

clean:
	$(RM) *.o shared/*.o

realclean: clean
	$(RM) *~ *.exe

pltclean:
	$(RM) -rf plt*/ *.png *_diagnostics.txt __pycache__
