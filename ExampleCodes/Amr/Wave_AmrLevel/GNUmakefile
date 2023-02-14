AMREX_HOME ?= ../../../../amrex

DEBUG = FALSE

DIM = 2

USE_MPI = TRUE
USE_OMP = FALSE

USE_CUDA  = FALSE
USE_HIP   = FALSE
USE_SYCL  = FALSE

#No Fortran
BL_NO_FORT = TRUE
# No Fortran probin file
AMREX_NO_PROBINIT = TRUE

CXXSTD=c++17

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

