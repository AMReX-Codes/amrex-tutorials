# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME ?= ../../../../../amrex

DEBUG        = FALSE
USE_MPI      = TRUE
USE_OMP      = FALSE
COMP         = gnu
DIM          = 3
USE_SUNDIALS = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ../Source/Make.package
VPATH_LOCATIONS  += ../Source
INCLUDE_LOCATIONS += ../Source

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
