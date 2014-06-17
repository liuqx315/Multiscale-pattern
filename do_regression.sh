#!/bin/bash
# ----------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ SMU
# ----------------------------------------------------------------
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
# ----------------------------------------------------------------


# update the ARKode_dev source code
hg pull
hg up

# create build directory
mkdir build_opt
cd build_opt

# configure SUNDIALS
cmake \
-DKLU_ENABLE:BOOL="1" \
-DOPENMP_ENABLE:BOOL="1" \
-DPTHREAD_ENABLE:BOOL="1" \
-DCXX_ENABLE:BOOL="1" \
-DLAPACK_ENABLE:BOOL="1" \
-DF90_ENABLE:BOOL="1" \
-DMPI_ENABLE:BOOL="1" \
-DFCMIX_ENABLE:BOOL="1" \
-DEXAMPLES_INSTALL_PATH:STRING="$PWD/../examples_opt" \
-DCMAKE_INSTALL_PREFIX:PATH="$PWD/../install_opt" \
-DCMAKE_C_COMPILER:FILEPATH="/usr/bin/mpicc" \
-DCMAKE_CXX_COMPILER:FILEPATH="/usr/bin/mpicxx" \
-DCMAKE_Fortran_COMPILER:FILEPATH="/usr/bin/mpif90" \
-DMPI_RUN_COMMAND:STRING="mpirun" \
-DMPI_MPICC:FILEPATH="/usr/bin/mpicc" \
-DMPI_MPICXX:FILEPATH="/usr/bin/mpicxx" \
-DMPI_MPIF77:FILEPATH="/usr/bin/mpif77" \
-DMPI_MPIF90:FILEPATH="/usr/bin/mpif90" \
-DCMAKE_C_FLAGS:STRING="-Wall -ansi -pedantic -O2" \
-DCMAKE_C_FLAGS_DEBUG:STRING="-g" \
-DCMAKE_C_FLAGS_RELEASE:STRING="-O3 -DNDEBUG" \
-DCMAKE_CXX_FLAGS:STRING="-Wall -O2" \
-DCMAKE_CXX_FLAGS_DEBUG:STRING="-g" \
-DCMAKE_CXX_FLAGS_RELEASE:STRING="-O3 -DNDEBUG" \
-DCMAKE_Fortran_FLAGS:STRING="-Wall -O2" \
-DCMAKE_Fortran_FLAGS_DEBUG:STRING="-g" \
-DCMAKE_Fortran_FLAGS_RELEASE:STRING="-O3" \
-DOpenMP_C_FLAGS:STRING="-fopenmp" \
-DOpenMP_CXX_FLAGS:STRING="-fopenmp" \
-DSUNDIALS_RT_LIBRARY:FILEPATH="-L/usr/lib/x86_64-linux-gnu -lrt" \
-DAMD_LIBRARY:FILEPATH="-L/usr/local/klu-1.2.0/gnu/lib -lamd" \
-DKLU_INCLUDE_DIR:STRING="/usr/local/klu-1.2.0/gnu/include" \
-DBTF_LIBRARY:FILEPATH="-L/usr/local/klu-1.2.0/gnu/lib -lbtf" \
-DKLU_LIBRARY:FILEPATH="-L/usr/local/klu-1.2.0/gnu/lib -lklu" \
-DCOLAMD_LIBRARY:FILEPATH="-L/usr/local/klu-1.2.0/gnu/lib -lcolamd" \
-DLAPACK_LIBRARIES:STRING="-L/usr/lib -llapack -L/usr/lib/libblas -lblas" \
../sundials

# build and install sundials
make
make install

# build regression tests
#   arkode
cd regression/arkode
make all klu
cd -
#   cvode
cd regression/cvode
make 
cd -
#   ida
cd regression/ida
make 
cd -
#   sundials
cd regression/sundials
make 
cd -


# run arkode regression tests
cd regression/arkode
./regression_runner.py -q -l
cd -
