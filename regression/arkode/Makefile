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
CXX = mpicxx
F90 = mpif90

# SUNDIALS definitions (assumes standard installation)
SUND_DEFS = -DSUNDIALS_CASE_LOWER -DSUNDIALS_UNDERSCORE_ONE \
            -DSUNDIALS_DOUBLE_PRECISION

# compilation flags 
CFLAGS = -O2
FFLAGS = -O2
CXXFLAGS = -O2 -DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX
# CFLAGS = -O0 -g
# FFLAGS = -O0 -g
# CXXFLAGS = -O0 -g -DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX
OMPFLAGS = -fopenmp

# append problem and interface preprocessor defs to compilation defs.
CFLAGS_ = $(CFLAGS) $(SUND_DEFS)
CXXFLAGS_ = $(CXXFLAGS) $(SUND_DEFS)

# shortcuts for include and library paths, etc.
INCS = -I../../install_opt/include/
LIBS = ../../install_opt/lib/libsundials_arkode.a \
       ../../install_opt/lib/libsundials_nvecparallel.a \
       ../../install_opt/lib/libsundials_nvecserial.a \
       -lm
FLIBS = ../../install_opt/lib/libsundials_farkode.a \
        ../../install_opt/lib/libsundials_fnvecserial.a \
        ../../install_opt/lib/libsundials_arkode.a \
        ../../install_opt/lib/libsundials_nvecserial.a
KLULIBS = /usr/local/klu-1.2.0/gnu/lib/libklu.a \
          /usr/local/klu-1.2.0/gnu/lib/libcolamd.a \
          /usr/local/klu-1.2.0/gnu/lib/libamd.a \
          /usr/local/klu-1.2.0/gnu/lib/libbtf.a
SLULIBS = -L/usr/local/superlu-mt-2.1/gnu_pthread/lib -lsuperlu_mt \
          -L/usr/lib/libblas -lblas
# INCS = -I../../install_dbg/include/
# LIBS = ../../install_dbg/lib/libsundials_arkode.a \
#        ../../install_dbg/lib/libsundials_nvecparallel.a \
#        ../../install_dbg/lib/libsundials_nvecserial.a \
#        -lm
# FLIBS = ../../install_dbg/lib/libsundials_farkode.a \
#         ../../install_dbg/lib/libsundials_fnvecserial.a \
#         ../../install_dbg/lib/libsundials_arkode.a \
#         ../../install_dbg/lib/libsundials_nvecserial.a
# KLULIBS = /usr/local/klu-1.2.0/gnu_debug/lib/libklu.a \
#           /usr/local/klu-1.2.0/gnu_debug/lib/libcolamd.a \
#           /usr/local/klu-1.2.0/gnu_debug/lib/libamd.a \
#           /usr/local/klu-1.2.0/gnu_debug/lib/libbtf.a
# SLULIBS = -L/usr/local/superlu-mt-2.1/gnu_pthread_debug/lib -lsuperlu_mt \
#           -L/usr/lib/libblas -lblas

# listing of all test routines
TESTS = ark_analytic.exe \
        ark_analytic_sys.exe \
        ark_analytic_nonlin.exe \
        ark_analytic_nonlin_back.exe \
        ark_brusselator.exe \
        ark_brusselator_Ma.exe \
        ark_brusselator_Mb.exe \
        ark_brusselator_Mc.exe \
        ark_brusselator_Md.exe \
        ark_brusselator_Me.exe \
        ark_brusselator1D.exe \
        ark_brusselator1D_FEM.exe \
        ark_heat1D.exe \
        ark_heat1D_adapt.exe \
        ark_heat2D.exe \
        ark_robertson.exe \
        ark_robertson_root.exe \
        ark_set_check.exe \
        ark_bruss.exe \
        ark_hires.exe \
        ark_medakzo.exe \
        ark_orego.exe \
        ark_pollu.exe \
        ark_ringmod.exe \
        ark_rober.exe \
	ark_vdpol.exe \
	ark_vdpolm.exe

KLUTESTS = ark_brusselator1D_klu.exe \
           ark_brusselator1D_FEM_klu.exe \
           ark_robertson_klu.exe

SLUTESTS = ark_brusselator1D_slu.exe \
           ark_brusselator1D_FEM_slu.exe \
           ark_robertson_slu.exe

OMPTESTS = ark_brusselator1D_omp.exe

PTTESTS = ark_brusselator1D_pt_slu.exe


# target to build all test executables
all : ${TESTS}

# target to build all SuperLU_MT test executables
klu : ${KLUTESTS}

# target to build all KLU test executables
slu : ${SLUTESTS}

# target to build all OpenMP test executables
omp : ${OMPTESTS}

# target to build all Pthreads test executables
pt : ${PTTESTS}


# specific build rules for C programs linking with KLU
ark_robertson_klu.exe : ark_robertson_klu.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${KLULIBS} ${LDFLAGS} -o $@

ark_brusselator1D_klu.exe : ark_brusselator1D_klu.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${KLULIBS} ${LDFLAGS} -o $@

ark_brusselator1D_FEM_klu.exe : ark_brusselator1D_FEM_klu.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${KLULIBS} ${LDFLAGS} -o $@

# specific build rules for C programs linking with SuperLU_MT
ark_robertson_slu.exe : ark_robertson_slu.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${SLULIBS} ${LDFLAGS} -o $@

ark_brusselator1D_slu.exe : ark_brusselator1D_slu.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${SLULIBS} ${LDFLAGS} -o $@

ark_brusselator1D_FEM_slu.exe : ark_brusselator1D_FEM_slu.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${SLULIBS} ${LDFLAGS} -o $@

ark_brusselator1D_pt_slu.exe : ark_brusselator1D_pt_slu.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ../../install_opt/lib/libsundials_nvecpthreads.a ${SLULIBS} ${LDFLAGS} -o $@

# specific build rules for C programs using OpenMP
ark_brusselator1D_omp.exe : ark_brusselator1D_omp.c init_from_file.c
	${CC} ${OMPFLAGS} ${CFLAGS_} ${INCS} $^ ${LIBS} ../../install_opt/lib/libsundials_nvecopenmp.a ${LDFLAGS} -o $@


# specific build rules for regression testing with different preprocessor directives
ark_brusselator_Ma.exe : ark_brusselator_M.c init_from_file.c   # dense-dense
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

ark_brusselator_Mb.exe : ark_brusselator_M.c init_from_file.c   # spgmr-spgmr
	${CC} ${CFLAGS_} ${INCS} $^ -DUSE_ITERATIVE -DUSE_SPGMR -DMASS_USE_ITERATIVE -DMASS_USE_SPGMR ${LIBS} ${LDFLAGS} -o $@

ark_brusselator_Mc.exe : ark_brusselator_M.c init_from_file.c   # spbcg-pcg
	${CC} ${CFLAGS_} ${INCS} $^ -DUSE_ITERATIVE -DUSE_SPBCG -DMASS_USE_ITERATIVE -DMASS_USE_PCG ${LIBS} ${LDFLAGS} -o $@

ark_brusselator_Md.exe : ark_brusselator_M.c init_from_file.c   # sptfqmr-spbcg
	${CC} ${CFLAGS_} ${INCS} $^ -DUSE_ITERATIVE -DUSE_SPTFQMR -DMASS_USE_ITERATIVE -DMASS_USE_SPBCG ${LIBS} ${LDFLAGS} -o $@

ark_brusselator_Me.exe : ark_brusselator_M.c init_from_file.c   # spgmr-dense
	${CC} ${CFLAGS_} ${INCS} $^ -DUSE_ITERATIVE -DUSE_SPGMR ${LIBS} ${LDFLAGS} -o $@

# specific build rules for Fortran programs
ark_bruss.exe : ark_bruss.f90
	${F90} ${FFLAGS} $^ ${FLIBS} ${LDFLAGS} -o $@

ark_hires.exe : ark_hires.f90 hires.f
	${F90} ${FFLAGS} $^ ${FLIBS} ${LDFLAGS} -o $@

ark_medakzo.exe : ark_medakzo.f90 medakzo.f
	${F90} ${FFLAGS} $^ ${FLIBS} ${LDFLAGS} -o $@

ark_orego.exe : ark_orego.f90 orego.f
	${F90} ${FFLAGS} $^ ${FLIBS} ${LDFLAGS} -o $@

ark_pollu.exe : ark_pollu.f90 pollu.f
	${F90} ${FFLAGS} $^ ${FLIBS} ${LDFLAGS} -o $@

ark_ringmod.exe : ark_ringmod.f90 ringmod.f
	${F90} ${FFLAGS} $^ ${FLIBS} ${LDFLAGS} -o $@

ark_rober.exe : ark_rober.f90 rober.f
	${F90} ${FFLAGS} $^ ${FLIBS} ${LDFLAGS} -o $@

ark_vdpol.exe : ark_vdpol.f90 vdpol.f
	${F90} ${FFLAGS} $^ ${FLIBS} ${LDFLAGS} -o $@

ark_vdpolm.exe : ark_vdpolm.f90 vdpolm.f
	${F90} ${FFLAGS} $^ ${FLIBS} ${LDFLAGS} -o $@

# general build rules for C programs
%.exe : %.c init_from_file.c
	${CC} ${CFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

# general build rules for C++ programs
%.exe : %.cpp init_from_file.c
	${CXX} ${CXXFLAGS_} ${INCS} $^ ${LIBS} ${LDFLAGS} -o $@

clean :
	\rm -f *.o *.exe diags_*.txt bruss_*.txt heat*.txt output.txt


####### End of Makefile #######
