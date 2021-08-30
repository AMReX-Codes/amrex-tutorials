
USE_MPI  = TRUE
USE_OMP  = FALSE

COMP = gnu

DIM = 3

DEBUG = FALSE

AMREX_HOME ?= ../../..

USE_EB = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs := AmrCore 
Pdirs += Base 
Pdirs += Boundary 
Pdirs += EB
Pdirs += Particle

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
