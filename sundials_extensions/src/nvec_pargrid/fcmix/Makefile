###################################################################
#  Daniel R. Reynolds, SMU Mathematics
#  Copyright 2012; all rights reserved
###################################################################

# Make.machine is used for setting machine-dependent compilation
# variables: e.g. compilers, libraries, F90-to-C conversion vars
include ../../../Make.machine

# shortcuts for include and library paths, etc.
INCS = ${MPI_INCD} ${SUND_INCD} ${C_INCD} -I../../../include/ -I./

# objects
OBJ = fnvector_parallel_grid.o


# target to build all source files
all : ${OBJ}

%.o : %.c
	${CC} -c ${CFLAGS_} ${INCS} $^ -o $@

clean :
	\rm -f *.o



####### End of Makefile #######
