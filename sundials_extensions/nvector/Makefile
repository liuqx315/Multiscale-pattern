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
LIBS = ${SUND_LIB} ${MPI_LIB}

# objects
OBJ = nvector_parallel_ghost.o fnvector_parallel_ghost.o

# target to build all source files
all : ${OBJ}

nvector_parallel_ghost.o : nvector_parallel_ghost.c nvector_parallel_ghost.h
	${CC} -c ${CFLAGS_} ${INCS} nvector_parallel_ghost.c -o $@

fnvector_parallel_ghost.o : fnvector_parallel_ghost.c fnvector_parallel_ghost.h nvector_parallel_ghost.h
	${CC} -c ${CFLAGS_} ${INCS} fnvector_parallel_ghost.c -o $@

clean :
	\rm -f *.o



####### End of Makefile #######
