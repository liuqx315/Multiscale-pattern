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
TESTS = cv_analytic.exe \
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
        ark_brusselator1D.exe \
        ark_bruss.exe \
        ark_heat1D.exe \
        ark_robertson.exe \
        ark_robertson_root.exe

# target to build all test executables
all : ${TESTS}

%.exe : %.c 
	${CC} ${CFLAGS_} ${INCS} $< ${LIBS} ${LDFLAGS} -o $@

%.exe : %.cpp
	${CXX} ${CFLAGS_} ${INCS} $< ${LIBS} ${LDFLAGS} -o $@

%.exe : %.f90
	${F90} ${FFLAGS} $^ ${LIBS} ${LDFLAGS} -o $@

clean :
	\rm -f *.o *.exe


####### End of Makefile #######
