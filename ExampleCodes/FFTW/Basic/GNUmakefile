# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME ?= ../../../amrex

USE_MPI      = TRUE
DEBUG        = FALSE
USE_OMP      = FALSE
COMP         = gnu
DIM          = 2

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

ifeq ($(USE_CUDA),TRUE)
  LIBRARIES += -lcufft
else
  LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3
endif
