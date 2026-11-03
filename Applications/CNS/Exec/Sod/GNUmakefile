AMREX_HOME ?= ../../../../../amrex

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

EBASE = CNS

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

# Problem-specific files in this directory
include ./Make.package

# CNS sources shared by all problems
include ../../Source/Make.package
VPATH_LOCATIONS   += ../../Source
INCLUDE_LOCATIONS += ../../Source

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
