AMREX_HOME ?= ../../../../amrex

DEBUG = FALSE
DIM = 3
COMP = gcc
TINY_PROFILE = TRUE
USE_MPI = TRUE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_FFT = TRUE

BL_NO_FORT = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/FFT/Make.package
include Make.package


include $(AMREX_HOME)/Tools/GNUMake/Make.rules
