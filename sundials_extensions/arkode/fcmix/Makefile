###################################################################
#  Daniel R. Reynolds, SMU Mathematics
#  Copyright 2012; all rights reserved
###################################################################

# Make.machine is used for setting machine-dependent compilation
# variables: e.g. compilers, libraries, F90-to-C conversion vars
include ../../Make.machine

# shortcuts for include and library paths, etc.
INCS = ${MPI_INCD} ${SUND_INCD} ${C_INCD} -I../ -I../../

# objects
OBJ = farkband.o \
      farkbbd.o \
      farkbp.o \
      farkdense.o \
      farkewt.o \
      farkjtimes.o \
      farklapack.o \
      farkode.o \
      farkpreco.o \
      farkroot.o


# target to build all source files
all : ${OBJ}

farkband.o : farkband.c farkode.h ../arkode_impl.h ../arkode_band.h
	${CC} -c ${CFLAGS_} ${INCS} farkband.c -o $@

farkbbd.o : farkbbd.c farkode.h farkbbd.h ../arkode_bbdpre.h ../arkode_sptfqmr.h ../arkode_spbcgs.h ../arkode_spgmr.h
	${CC} -c ${CFLAGS_} ${INCS} farkbbd.c -o $@

farkbp.o : farkbp.c farkode.h farkbp.h ../arkode_bandpre.h
	${CC} -c ${CFLAGS_} ${INCS} farkbp.c -o $@

farkdense.o : farkdense.c farkode.h ../arkode_impl.h ../arkode_dense.h
	${CC} -c ${CFLAGS_} ${INCS} farkdense.c -o $@

farkewt.o : farkewt.c farkode.h ../arkode_impl.h
	${CC} -c ${CFLAGS_} ${INCS} farkewt.c -o $@

farkjtimes.o : farkjtimes.c farkode.h ../arkode_impl.h ../arkode_spils.h
	${CC} -c ${CFLAGS_} ${INCS} farkjtimes.c -o $@

farklapack.o : farklapack.c farkode.h ../arkode_impl.h ../arkode_lapack.h
	${CC} -c ${CFLAGS_} ${INCS} farklapack.c -o $@

farkode.o : farkode.c farkode.h ../arkode_impl.h ../arkode_band.h ../arkode_dense.h ../arkode_spgmr.h ../arkode_spbcgs.h ../arkode_sptfqmr.h
	${CC} -c ${CFLAGS_} ${INCS} farkode.c -o $@

farkpreco.o : farkpreco.c farkode.h ../arkode_impl.h ../arkode_spils.h
	${CC} -c ${CFLAGS_} ${INCS} farkpreco.c -o $@

farkroot.o : farkroot.c farkode.h farkroot.h ../arkode_impl.h
	${CC} -c ${CFLAGS_} ${INCS} farkroot.c -o $@



clean :
	\rm -f *.o



####### End of Makefile #######
