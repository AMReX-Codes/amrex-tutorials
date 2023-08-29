AMREX_HOME ?= ../../amrex
HEFFTE_HOME ?= ../../heffte/build

DEBUG = FALSE
DIM = 3
COMP = gcc
TINY_PROFILE = TRUE
USE_MPI = TRUE
USE_CUDA = TRUE

BL_NO_FORT = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include $(AMREX_HOME)/Src/Base/Make.package
include Make.package

VPATH_LOCATIONS += $(HEFFTE_HOME)/include
INCLUDE_LOCATIONS += $(HEFFTE_HOME)/include
LIBRARY_LOCATIONS += $(HEFFTE_HOME)/lib

libraries += -lheffte

ifeq ($(USE_CUDA),TRUE)
  libraries += -lcufft
else
  libraries += -lfftw3_mpi -lfftw3f -lfftw3
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
