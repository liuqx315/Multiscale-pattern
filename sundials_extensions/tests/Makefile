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
        cv_analytic.exe

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

clean :
	\rm -f *.o *.exe


####### End of Makefile #######
