###################################################################
#  Daniel R. Reynolds, SMU Mathematics
#  Copyright 2012; all rights reserved
###################################################################

# Make.machine is used for setting machine-dependent compilation
# variables: e.g. compilers, libraries, F90-to-C conversion vars
include ../../Make.machine

# append problem and interface preprocessor defs to compilation defs.
CFLAGS_ = $(CFLAGS) $(SUND_DEFS)
CXXFLAGS_ = $(CXXFLAGS) $(SUND_DEFS)

# shortcuts for include and library paths, etc.
INCS = ${MPI_INCD} ${SUND_INCD} ${C_INCD} -I../../include/
LIBS = ../../sundials_extensions.a ${SUND_LIB} ${MPI_LIB}

# listing of all test routines
TESTS = ark_analytic.exe \
        ark_analytic_sys.exe \
        ark_analytic_nonlin.exe \
        ark_analytic_nonlin_back.exe \
        ark_brusselator.exe \
        ark_brusselator1D.exe \
        ark_heat1D.exe \
        ark_heat2D.exe \
        ark_robertson.exe \
        ark_robertson_root.exe \
        ark_set_check.exe \
        ark_bruss.exe \
	ark_emep.exe \
        ark_hires.exe \
        ark_medakzo.exe \
        ark_orego.exe \
        ark_pollu.exe \
        ark_ringmod.exe \
        ark_rober.exe \
	ark_vdpol.exe \
	ark_vdpolm.exe

# target to build all test executables
all : ${TESTS}

%.exe : %.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

%.exe : %.cpp init_from_file.c
	${CXX} ${CXXFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_bruss.exe : ark_bruss.f90
	${F90} ${FFLAGS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_emep.exe : ark_emep.f90 emep.f
	${F90} ${FFLAGS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_hires.exe : ark_hires.f90 hires.f
	${F90} ${FFLAGS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_medakzo.exe : ark_medakzo.f90 medakzo.f
	${F90} ${FFLAGS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_orego.exe : ark_orego.f90 orego.f
	${F90} ${FFLAGS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_pollu.exe : ark_pollu.f90 pollu.f
	${F90} ${FFLAGS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_ringmod.exe : ark_ringmod.f90 ringmod.f
	${F90} ${FFLAGS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_rober.exe : ark_rober.f90 rober.f
	${F90} ${FFLAGS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_vdpol.exe : ark_vdpol.f90 vdpol.f
	${F90} ${FFLAGS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_vdpolm.exe : ark_vdpolm.f90 vdpolm.f
	${F90} ${FFLAGS} $^ ${LIBS} ${LDFLAGS} -o $@

clean :
	\rm -f *.o *.exe


####### End of Makefile #######
