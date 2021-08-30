#AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../../../amrex

DEBUG      = FALSE
USE_MPI    = TRUE
USE_OMP    = FALSE
USE_CUDA   = FALSE
COMP	   = gnu
DIM        = 2

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ../Source/Make.package
VPATH_LOCATIONS  += ../Source
INCLUDE_LOCATIONS += ../Source

# Pytorch directories 
PYTORCH_ROOT := ../libtorch
TORCH_LIBPATH = $(PYTORCH_ROOT)/lib

ifeq ($(USE_CUDA),TRUE)
  TORCH_LIBS = -ltorch -ltorch_cpu -lc10 -lc10_cuda -lcuda 
else
  TORCH_LIBS = -ltorch -ltorch_cpu -lc10
endif

INCLUDE_LOCATIONS += $(PYTORCH_ROOT)/include \
                     $(PYTORCH_ROOT)/include/torch/csrc/api/include
LIBRARY_LOCATIONS += $(TORCH_LIBPATH)

DEFINES += -D_GLIBCXX_USE_CXX11_ABI=1
ifeq ($(USE_CUDA),TRUE)
  LDFLAGS += -Xlinker "--no-as-needed,-rpath $(TORCH_LIBPATH) $(TORCH_LIBS)"
else
  LDFLAGS += -Wl,--no-as-needed,-rpath=$(TORCH_LIBPATH) $(TORCH_LIBS)
endif

# AMReX 
Pdirs   := Base Boundary AmrCore 
Ppack   += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
