AMREX_HOME ?= ../../..

USE_MPI  = TRUE
USE_OMP  = FALSE

COMP = gnu

DIM = 2
#DIM = 3

DEBUG = FALSE

USE_EB = TRUE
USE_PARTICLES = TRUE

USE_HYPRE = TRUE
USE_HYPRE = FALSE

include ./Make.Amr

movie:
	ffmpeg -framerate 10 -pattern_type glob -i "p*.png" -r 5 -vf scale=1024:-1 amr102.gif
