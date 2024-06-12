# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME ?= ../../../../amrex

DEBUG        = FALSE
USE_MPI      = TRUE
USE_OMP      = FALSE
COMP         = gnu
DIM          = 3
USE_RPATH    = TRUE
USE_SUNDIALS = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ../Source/Make.package
VPATH_LOCATIONS  += ../Source
INCLUDE_LOCATIONS += ../Source

ifeq ($(USE_SUNDIALS),TRUE)
  ifeq ($(USE_CUDA),TRUE)
    SUNDIALS_ROOT ?= $(TOP)../../../../sundials/instdir_cuda
  else
    SUNDIALS_ROOT ?= $(TOP)../../../../sundials/instdir
  endif
  ifeq ($(NERSC_HOST),perlmutter)
    SUNDIALS_LIB_DIR ?= $(SUNDIALS_ROOT)/lib64
  else
    SUNDIALS_LIB_DIR ?= $(SUNDIALS_ROOT)/lib
  endif

  USE_CVODE_LIBS ?= TRUE
  USE_ARKODE_LIBS ?= TRUE

  DEFINES += -DAMREX_USE_SUNDIALS
  INCLUDE_LOCATIONS += $(SUNDIALS_ROOT)/include
  LIBRARY_LOCATIONS += $(SUNDIALS_LIB_DIR)

  LIBRARIES += -L$(SUNDIALS_LIB_DIR) -lsundials_cvode
  LIBRARIES += -L$(SUNDIALS_LIB_DIR) -lsundials_arkode
  LIBRARIES += -L$(SUNDIALS_LIB_DIR) -lsundials_nvecmanyvector
  LIBRARIES += -L$(SUNDIALS_LIB_DIR) -lsundials_core

  ifeq ($(USE_CUDA),TRUE)
    LIBRARIES += -L$(SUNDIALS_LIB_DIR) -lsundials_nveccuda
  endif

endif

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

