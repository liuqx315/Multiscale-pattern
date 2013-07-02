###################################################################
#  Daniel R. Reynolds, SMU Mathematics
#  Copyright 2013; all rights reserved
#
#  Regression directory Makefile: assumes that SUNDIALS+ARKode
#  has been built in this same source repository, with the 
#  directory structure:
#      <arkode>/regression  -- parent of this directory
#      <arkode>/sundials    -- SUNDIALS source files
#      <arkode>/install     -- SUNDIALS+ARKode install location
###################################################################

# compilers (assumes standard MPICH/OpenMPI compiler wrappers)
CC = mpicc

# SUNDIALS definitions (assumes standard installation)
SUND_DEFS = -DSUNDIALS_CASE_LOWER -DSUNDIALS_UNDERSCORE_ONE \
            -DSUNDIALS_DOUBLE_PRECISION

# compilation flags 
CFLAGS = -O2

# append problem and interface preprocessor defs to compilation defs.
CFLAGS_ = $(CFLAGS) $(SUND_DEFS)

# shortcuts for include and library paths, etc.
INCS = -I../../install/include/
LIBS = -L../../install/lib -lsundials_cvode -lsundials_nvecserial

# listing of all test routines
TESTS = cv_analytic.exe \
        cv_analytic_sys.exe \
        cv_analytic_nonlin.exe \
        cv_VanderPol.exe \
        cv_brusselator.exe \
        cv_brusselator1D.exe 

# target to build all test executables
all : ${TESTS}

%.exe : %.c 
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

clean :
	\rm -f *.o *.exe


####### End of Makefile #######
