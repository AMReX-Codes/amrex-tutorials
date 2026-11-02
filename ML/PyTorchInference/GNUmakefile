AMREX_HOME ?= ../../../amrex

DEBUG        ?= FALSE
DIM          ?= 3
COMP         ?= gnu
USE_MPI      ?= TRUE
USE_OMP      ?= FALSE
USE_CUDA     ?= FALSE
TINY_PROFILE ?= FALSE
BL_NO_FORT   ?= TRUE

# LibTorch directory (with include/ and lib/), e.g., an unpacked LibTorch download
# or the torch Python package: python3 -c "import torch, os; print(os.path.dirname(torch.__file__))"
PYTORCH_ROOT ?=
ifeq ($(PYTORCH_ROOT),)
  ifeq ($(filter clean realclean,$(MAKECMDGOALS)),)
    $(error Set PYTORCH_ROOT to the LibTorch directory)
  endif
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

TORCH_ABSPATH := $(abspath $(PYTORCH_ROOT))
SYSTEM_INCLUDE_LOCATIONS += $(TORCH_ABSPATH)/include $(TORCH_ABSPATH)/include/torch/csrc/api/include
LIBRARY_LOCATIONS += $(TORCH_ABSPATH)/lib
LIBRARIES += -Wl,-rpath,$(TORCH_ABSPATH)/lib -Wl,--no-as-needed -ltorch -ltorch_cpu -lc10
ifeq ($(USE_CUDA),TRUE)
  LIBRARIES += -ltorch_cuda -lc10_cuda
endif
LIBRARIES += -Wl,--as-needed

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
