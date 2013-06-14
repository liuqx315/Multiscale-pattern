###################################################################
#  Daniel R. Reynolds, SMU Mathematics
#  Copyright 2012; all rights reserved
###################################################################

# Make.machine is used for setting machine-dependent compilation
# variables: e.g. compilers, libraries, F90-to-C conversion vars
include ../../Make.machine

# shortcuts for include and library paths, etc.
INCS = ${MPI_INCD} ${SUND_INCD} ${C_INCD} -I./ -I../ -I../../

# objects
OBJ = sundials_spfgmr.o

# target to build all source files
all : ${OBJ}

sundials_spfgmr.o : sundials_spfgmr.c sundials_spfgmr.h
	${CC} -c ${CFLAGS_} ${INCS} $< -o $@


clean :
	\rm -f *.o



####### End of Makefile #######
