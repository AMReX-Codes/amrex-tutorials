ALL: boundary_control

SPACK_PKGS       = /home/adener/FASTMath/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5
AMREX_HOME       = ${SPACK_PKGS}/amrex-develop-fvkknjb363aqxpt6zu2tt6oojixnlnex

CFLAGS           =
FFLAGS           =
CPPFLAGS         = -I${AMREX_HOME}/include
FPPFLAGS         =
CLEANFILES       = boundary_control 
PETSC_DIR        = ${SPACK_PKGS}/petsc-develop-muikvcmmgnpc2r7socr77wu55mcmhmxi

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

LDFLAGS	        += -L${AMREX_HOME}/lib
LDLIBS          += -lamrex

boundary_control: boundary_control.o MyTest.o MyTestPlotfile.o initProb.o AMReX_BCUtil.o