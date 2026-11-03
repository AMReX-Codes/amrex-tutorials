AMREX_HOME ?= ../../../amrex

DEBUG        ?= FALSE
DIM          ?= 3
COMP         ?= gnu
USE_MPI      ?= TRUE
USE_OMP      ?= FALSE
TINY_PROFILE ?= FALSE

USE_F_INTERFACES = TRUE

ifneq ($(filter TRUE,$(USE_CUDA) $(USE_HIP) $(USE_SYCL)),)
  ifeq ($(filter clean realclean,$(MAKECMDGOALS)),)
    $(error This tutorial runs on CPUs only)
  endif
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/F_Interfaces/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
