# Path to the AMReX source tree. Override with `make AMREX_HOME=/path/to/amrex`.
AMREX_HOME ?= ../../../../amrex

# Build options. Each can be overridden on the command line, e.g., `make USE_CUDA=TRUE`.
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

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

# Application sources
include ../Source/Make.package
VPATH_LOCATIONS   += ../Source
INCLUDE_LOCATIONS += ../Source

# AMReX sources. Add more directories as needed, e.g., Boundary AmrCore
# LinearSolvers/MLMG EB Particle FFT.
AMREX_DIRS = Base
include $(foreach dir, $(AMREX_DIRS), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
