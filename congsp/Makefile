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
CC = openmpicc
CXX = openmpicxx

# SUNDIALS definitions (assumes standard installation)
SUND_DEFS = -DSUNDIALS_CASE_LOWER -DSUNDIALS_UNDERSCORE_ONE \
            -DSUNDIALS_DOUBLE_PRECISION

# compilation flags 
CFLAGS = -O2
FFLAGS = -O2
CXXFLAGS = -O2 -DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX

# append problem and interface preprocessor defs to compilation defs.
CFLAGS_ = $(CFLAGS) $(SUND_DEFS)
CXXFLAGS_ = $(CXXFLAGS) $(SUND_DEFS)

# shortcuts for include and library paths, etc.
INCS = -I../install/include/
LIBS = -L../install/lib -lsundials_arkode -lsundials_nvecparallel -lsundials_nvecserial -lm

# listing of all test routines
TESTS = ARKode_WENO.exe \
	ARKode_WENO2D.exe

# target to build all test executables
#all : ${TESTS}
all : ARKode_WENO.exe ARKode_WENO2D.exe
#%.exe : %.cpp

ARKode_WENO.exe : ARKode_WENO.cpp
	${CXX} ${CXXFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@
ARKode_WENO2D.exe : ARKode_WENO_2D.cpp
	${CXX} ${CXXFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@
clean :
	\rm -f *.o *.exe


####### End of Makefile #######
