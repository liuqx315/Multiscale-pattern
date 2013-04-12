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
TESTS = nvector_serial_test.exe \
        nvector_parallel_test.exe \
        nvector_parallel_grid_test.exe \
        pcg_tester.exe \
        cv_analytic.exe \
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
        ark_heat1D.exe \
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

nvector_serial_test.exe : nvector_serial_test.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

nvector_parallel_test.exe : nvector_parallel_test.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

nvector_parallel_grid_test.exe : nvector_parallel_grid_test.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

pcg_tester.exe : pcg_tester.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

cv_analytic.exe : cv_analytic.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

cv_analytic_sys.exe : cv_analytic_sys.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

cv_analytic_nonlin.exe : cv_analytic_nonlin.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

cv_VanderPol.exe : cv_VanderPol.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

cv_brusselator.exe : cv_brusselator.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

cv_brusselator1D.exe : cv_brusselator1D.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

ida_brusselator1D_FEM.exe : ida_brusselator1D_FEM.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_analytic.exe : ark_analytic.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_analytic_sys.exe : ark_analytic_sys.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_analytic_nonlin.exe : ark_analytic_nonlin.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_brusselator.exe : ark_brusselator.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_brusselator1D.exe : ark_brusselator1D.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_heat1D.exe : ark_heat1D.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_robertson.exe : ark_robertson.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_robertson_root.exe : ark_robertson_root.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_set_check.exe : ark_set_check.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

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
