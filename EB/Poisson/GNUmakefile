AMREX_HOME ?= ../../../amrex

DEBUG        ?= FALSE
DIM          ?= 3
COMP         ?= gnu
USE_MPI      ?= TRUE
USE_OMP      ?= FALSE
USE_CUDA     ?= FALSE
USE_HIP      ?= FALSE
USE_SYCL     ?= FALSE
TINY_PROFILE ?= FALSE
BL_NO_FORT   ?= TRUE

USE_EB = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
