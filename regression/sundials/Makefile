###################################################################
#  Programmer(s):  Daniel R. Reynolds @ SMU
###################################################################
#  Copyright (c) 2013, Southern Methodist University.
#  All rights reserved.
#  For details, see the LICENSE file.
###################################################################
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
INCS = -I../../install_opt/include/
LIBS = -L../../install_opt/lib -lsundials_arkode -lsundials_nvecserial -lsundials_nvecparallel -lsundials_nvecparallelgrid

# listing of all test routines
TESTS = nvector_serial_test.exe \
        nvector_parallel_test.exe \
        nvector_parallel_grid_test.exe \
        pcg_tester.exe 

# target to build all test executables
all : ${TESTS}

%.exe : %.c 
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

clean :
	\rm -f *.o *.exe


####### End of Makefile #######
