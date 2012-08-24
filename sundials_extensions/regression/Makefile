###################################################################
#  Daniel R. Reynolds, SMU Mathematics
#  Copyright 2012; all rights reserved
###################################################################

# Make.machine is used for setting machine-dependent compilation
# variables: e.g. compilers, libraries, F90-to-C conversion vars
include ../Make.machine

# append problem and interface preprocessor defs to compilation defs.
CFLAGS_ = $(CFLAGS) $(SUND_DEFS)

# shortcuts for include and library paths, etc.
INCS = ${MPI_INCD} ${SUND_INCD} ${C_INCD} -I../
LIBS = ../sundials_extensions.a ${SUND_LIB} ${MPI_LIB}

# listing of all test routines
TESTS = nvector_serial_test.exe \
        nvector_parallel_test.exe \
        nvector_parallel_grid_test.exe \
        cv_analytic.exe \
        cv_analytic_sys.exe \
        cv_analytic_nonlin.exe \
        cv_VanderPol.exe \
        cv_brusselator.exe \
        cv_brusselator1D.exe \
        ida_brusselator1D_FEM.exe \
        ark_analytic.exe \
        ark_analytic_sys.exe \
        ark_analytic_nonlin.exe \
        ark_brusselator.exe \
        ark_brusselator1D.exe

# target to build all test executables
all : ${TESTS}

nvector_serial_test.exe : nvector_serial_test.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

nvector_parallel_test.exe : nvector_parallel_test.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

nvector_parallel_grid_test.exe : nvector_parallel_grid_test.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

cv_analytic.exe : cv_analytic.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

cv_analytic_sys.exe : cv_analytic_sys.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

cv_analytic_nonlin.exe : cv_analytic_nonlin.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

cv_VanderPol.exe : cv_VanderPol.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

cv_brusselator.exe : cv_brusselator.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

cv_brusselator1D.exe : cv_brusselator1D.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

ida_brusselator1D_FEM.exe : ida_brusselator1D_FEM.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

ark_analytic.exe : ark_analytic.c ark_setparams.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

ark_analytic_sys.exe : ark_analytic_sys.c ark_setparams.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

ark_analytic_nonlin.exe : ark_analytic_nonlin.c ark_setparams.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

ark_brusselator.exe : ark_brusselator.c ark_setparams.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

ark_brusselator1D.exe : ark_brusselator1D.c ark_setparams.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} -o $@

clean :
	\rm -f *.o *.exe


####### End of Makefile #######
