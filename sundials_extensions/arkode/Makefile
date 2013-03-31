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
OBJ = arkode.o \
      arkode_band.o \
      arkode_bandpre.o \
      arkode_bbdpre.o \
      arkode_butcher.o \
      arkode_dense.o \
      arkode_direct.o \
      arkode_io.o \
      arkode_lapack.o \
      arkode_pcg.o \
      arkode_spbcgs.o \
      arkode_sptfqmr.o \
      arkode_spgmr.o \
      arkode_spils.o


# target to build all source files
all : ${OBJ}

arkode.o : arkode.c arkode_impl.h
	${CC} -c ${CFLAGS_} ${INCS} arkode.c -o $@

arkode_band.o : arkode_band.c arkode_band.h arkode_impl.h arkode_direct_impl.h
	${CC} -c ${CFLAGS_} ${INCS} arkode_band.c -o $@

arkode_bandpre.o : arkode_bandpre.c arkode_impl.h arkode_bandpre_impl.h arkode_pcg.h arkode_spils_impl.h arkode_sptfqmr.h arkode_spbcgs.h arkode_spgmr.h
	${CC} -c ${CFLAGS_} ${INCS} arkode_bandpre.c -o $@

arkode_bbdpre.o : arkode_bbdpre.c arkode_impl.h arkode_bbdpre_impl.h arkode_pcg.h arkode_spils_impl.h arkode_sptfqmr.h arkode_spbcgs.h arkode_spgmr.h
	${CC} -c ${CFLAGS_} ${INCS} arkode_bbdpre.c -o $@

arkode_butcher.o : arkode_butcher.c arkode_impl.h
	${CC} -c ${CFLAGS_} ${INCS} arkode_butcher.c -o $@

arkode_dense.o : arkode_dense.c arkode_dense.h arkode_direct_impl.h arkode_impl.h
	${CC} -c ${CFLAGS_} ${INCS} arkode_dense.c -o $@

arkode_direct.o : arkode_direct.c arkode_impl.h arkode_direct_impl.h
	${CC} -c ${CFLAGS_} ${INCS} arkode_direct.c -o $@

arkode_io.o : arkode_io.c arkode_impl.h
	${CC} -c ${CFLAGS_} ${INCS} arkode_io.c -o $@

arkode_lapack.o : arkode_lapack.c arkode_impl.h arkode_direct_impl.h arkode_lapack.h
	${CC} -c ${CFLAGS_} ${INCS} arkode_lapack.c -o $@

arkode_pcg.o : arkode_pcg.c arkode_pcg.h arkode_spils_impl.h arkode_impl.h
	${CC} -c ${CFLAGS_} ${INCS} arkode_pcg.c -o $@

arkode_spbcgs.o : arkode_spbcgs.c arkode_spbcgs.h arkode_spils_impl.h arkode_impl.h
	${CC} -c ${CFLAGS_} ${INCS} arkode_spbcgs.c -o $@

arkode_sptfqmr.o : arkode_sptfqmr.c arkode_sptfqmr.h arkode_spils_impl.h arkode_impl.h
	${CC} -c ${CFLAGS_} ${INCS} arkode_sptfqmr.c -o $@

arkode_spgmr.o : arkode_spgmr.c arkode_spils.h
	${CC} -c ${CFLAGS_} ${INCS} arkode_spgmr.c -o $@

arkode_spils.o: arkode_spils.c
	${CC} -c ${CFLAGS_} ${INCS} arkode_spils.c -o $@


clean :
	\rm -f *.o



####### End of Makefile #######
