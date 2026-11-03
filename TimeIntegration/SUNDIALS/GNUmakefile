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

# Needs an installed SUNDIALS (>= 6), e.g., make SUNDIALS_DIR=/path/to/sundials/install
USE_SUNDIALS  = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
