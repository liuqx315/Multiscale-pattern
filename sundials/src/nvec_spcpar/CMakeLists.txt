# ---------------------------------------------------------------
# $Revision: 1.3 $
# $Date: 2009-02-17 02:58:48 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the parallel NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_SPCPAR\n\")")

IF(MPI_MPICC)
  # use MPI_MPICC as the compiler
  SET(CMAKE_C_COMPILER ${MPI_MPICC})
ELSE(MPI_MPICC)
  # add MPI_INCLUDE_PATH to include directories
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(MPI_MPICC)

# Add variable nvecspcparallel_SOURCES with the sources for the NVECSPCPARALLEL lib
SET(nvecspcparallel_SOURCES nvector_spcparallel.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECSPCPARALLEL library
SET(shared_SOURCES sundials_math.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/src/sundials/ shared_SOURCES)

# Add variable nvecspcparallel_HEADERS with the exported NVECSPCPARALLEL header files
SET(nvecspcparallel_HEADERS nvector_spcparallel.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nvecspcparallel_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECSPCPARALLEL library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECSPCPARALLEL library
IF(BUILD_STATIC_LIBS)
  ADD_LIBRARY(sundials_nvecspcparallel_static STATIC ${nvecspcparallel_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecspcparallel_static
    PROPERTIES OUTPUT_NAME sundials_nvecspcparallel CLEAN_DIRECT_OUTPUT 1)
  INSTALL(TARGETS sundials_nvecspcparallel_static DESTINATION lib)
ENDIF(BUILD_STATIC_LIBS)

# Rules for building and installing the shared library:
#  - Add the build target for the NVECSPCPARALLEL library
#  - Set the library name and make sure it is not deleted
#  - Set VERSION and SOVERSION for shared libraries
#  - Install the NVECSPCPARALLEL library
IF(BUILD_SHARED_LIBS)
  ADD_LIBRARY(sundials_nvecspcparallel_shared SHARED ${nvecspcparallel_SOURCES} ${shared_SOURCES})
  SET_TARGET_PROPERTIES(sundials_nvecspcparallel_shared
    PROPERTIES OUTPUT_NAME sundials_nvecspcparallel CLEAN_DIRECT_OUTPUT 1)
  SET_TARGET_PROPERTIES(sundials_nvecspcparallel_shared 
    PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
  INSTALL(TARGETS sundials_nvecspcparallel_shared DESTINATION lib)
ENDIF(BUILD_SHARED_LIBS)

# Install the NVECSPCPARALLEL header files
INSTALL(FILES ${nvecspcparallel_HEADERS} DESTINATION include/nvector)

# If FCMIX is enabled and MPI-F77 works, build and install the FNVECSPCPARALLEL library
IF(FCMIX_ENABLE AND MPIF_FOUND)
  SET(fnvecspcparallel_SOURCES fnvector_spcparallel.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fnvecspcparallel_static STATIC ${fnvecspcparallel_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecspcparallel_static
      PROPERTIES OUTPUT_NAME sundials_fnvecspcparallel CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fnvecspcparallel_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fnvecspcparallel_shared SHARED ${fnvecspcparallel_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecspcparallel_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecspcparallel CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fnvecspcparallel_shared
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    INSTALL(TARGETS sundials_fnvecspcparallel_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_ENABLE AND MPIF_FOUND)
#
MESSAGE(STATUS "Added NVECTOR_SPCPARALLEL module")
