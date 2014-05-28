/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds and Radu Serban @LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 */

#ifndef _NVECTOR_SPCPARALLEL_H
#define _NVECTOR_SPCPARALLEL_H

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <mpi.h>

#include <sundials/sundials_nvector.h>

/*
 * -----------------------------------------------------------------
 * PART I: SPCPARALLEL implementation of N_Vector               
 * -----------------------------------------------------------------
 */

/* define MPI data types */

#if defined(SUNDIALS_SINGLE_PRECISION)
#define SPVEC_REAL_MPI_TYPE MPI_FLOAT

#elif defined(SUNDIALS_DOUBLE_PRECISION)
#define SPVEC_REAL_MPI_TYPE MPI_DOUBLE

#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SPVEC_REAL_MPI_TYPE MPI_LONG_DOUBLE

#endif

#define SPVEC_INTEGER_MPI_TYPE MPI_LONG

/* The SPCPARALLEL implementation of the N_Vector 'content' structure 
 * contains the global length of the vector, local valid lengths and 
 * ghost dimensions in x, y and z-directions of the vector, a pointer 
 * to an array of real components, and the MPI communicator 
 */

struct _N_VectorContent_SpcParallel {
  int Ngrp;              /* number of variable groups               */
  int *Nspc;             /* number of species in each group         */
  long int Nx;           /* local x-mesh vector length              */
  long int Ny;           /* local y-mesh vector length              */
  long int Nz;           /* local z-mesh vector length              */
  long int NGx;          /* x-width of ghost boundary               */
  long int NGy;          /* y-width of ghost boundary               */
  long int NGz;          /* z-width of ghost boundary               */
  long int *n1;          /* local vector lengths for each group     */
  long int n;            /* local vector length                     */
  long int N;            /* global vector length                    */
  realtype *data;        /* local data array                        */
  realtype **gdata;      /* pointers in data at start of group data */
  booleantype own_data;  /* flag for ownership of data              */
  MPI_Comm comm;         /* pointer to MPI communicator             */
};

typedef struct _N_VectorContent_SpcParallel *N_VectorContent_SpcParallel;

/*
 * -----------------------------------------------------------------
 * PART II: Macros                                              
 * --------------------------------------------------------------
 * Notes
 *
 * When looping over the components of an N_Vector v, it is
 * more efficient to first obtain the component array via
 * v_data = NV_DATA_MHD(v) and then access v_data[i] within
 * the loop than it is to use NV_Ith_MHD(v,i) within the
 * loop.
 *
 * -----------------------------------------------------------------
 */

/* Vector content */

#define SPV_CONTENT(v) ( (N_VectorContent_SpcParallel)(v->content) )

/* Number of groups */

#define SPV_NGROUPS(v) (SPV_CONTENT(v)->Ngrp)

/* Number of species in group ig */

#define SPV_NSPECIES(v,ig) (SPV_CONTENT(v)->Nspc[ig])

/* Local grid dimensions */

#define SPV_XLENGTH(v) (SPV_CONTENT(v)->Nx)
#define SPV_YLENGTH(v) (SPV_CONTENT(v)->Ny)
#define SPV_ZLENGTH(v) (SPV_CONTENT(v)->Nz)

#define SPV_XGHOST(v) (SPV_CONTENT(v)->NGx)
#define SPV_YGHOST(v) (SPV_CONTENT(v)->NGy)
#define SPV_ZGHOST(v) (SPV_CONTENT(v)->NGz)

#define SPV_XGLENGTH(v) (SPV_XLENGTH(v)+2*SPV_XGHOST(v))
#define SPV_YGLENGTH(v) (SPV_YLENGTH(v)+2*SPV_YGHOST(v))
#define SPV_ZGLENGTH(v) (SPV_ZLENGTH(v)+2*SPV_ZGHOST(v))

/* Local data dimension for group ig */

#define SPV_LOCLENGTH1(v,ig) (SPV_CONTENT(v)->n1[ig])

/* Local data dimension */

#define SPV_LOCLENGTH(v)  (SPV_CONTENT(v)->n)

/* Global data dimension */

#define SPV_GLOBLENGTH(v) (SPV_CONTENT(v)->N)

/* Data access */

#define SPV_DATA(v) (SPV_CONTENT(v)->data)

#define SPV_GDATA(v,ig) (SPV_CONTENT(v)->gdata[ig])

#define SPV_Ith(v,ig,i) (SPV_DATA(v)[i])

#define SPV_GIth(v,ig,i) (SPV_GDATA(v,ig)[i])

/* MPI communicator */

#define SPV_COMM(v) (SPV_CONTENT(v)->comm)

/* Data ownership flag */

#define SPV_OWN_DATA(v) (SPV_CONTENT(v)->own_data)

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by SPCPARALLEL
 * 
 * CONSTRUCTORS:
 *    N_VNew_SpcParallel
 *    N_VNewEmpty_SpcParallel
 *    N_VMake_SpcParallel
 *    N_VLoad_SpcParallel
 *    N_VCloneVectorArray_SpcParallel        
 *    N_VCloneVectorArrayEmpty_SpcParallel   
 * DESTRUCTORS:
 *    N_VDestroy_SpcParallel              
 *    N_VDestroyVectorArray_SpcParallel    
 * OTHER:
 *    N_VPrint_SpcParallel
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : N_VNew_SpcParallel
 * -----------------------------------------------------------------
 * This function creates and allocates memory for an 
 * SPCPARALLEL vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNew_SpcParallel(MPI_Comm comm, int Ngrp, int *Nspc, 
					    long int Nx,  long int Ny,  long int Nz, 
					    long int NGx, long int NGy, long int NGz);
/*
 * -----------------------------------------------------------------
 * Function : N_VNewEmpty_SpcParallel
 * -----------------------------------------------------------------
 * This function creates a new SPCPARALLEL N_Vector with an empty
 * (NULL) data array.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_SpcParallel(MPI_Comm comm,  int Ngrp, int *Nspc,
						 long int Nx,  long int Ny,  long int Nz, 
						 long int NGx, long int NGy, long int NGz);

/*
 * -----------------------------------------------------------------
 * Function : N_VMake_SpcParallel
 * -----------------------------------------------------------------
 * This function creates an SPCPARALLEL vector and attaches to it
 * the user-supplied data (which must be a contiguous array)
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VMake_SpcParallel(MPI_Comm comm, int Ngrp, int *Nspc, 
					     long int Nx,  long int Ny,  long int Nz, 
					     long int NGx, long int NGy, long int NGz,
					     realtype *data);

/*
 * -----------------------------------------------------------------
 * Function : N_VLoad_SpcParallel
 * -----------------------------------------------------------------
 * This function copies into a given SPCPARALLEL N_Vector the data
 * for one species. The user must indicate the group index and the
 * species index within that group.
 * The data must be a realtype array in the correct order.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VLoad_SpcParallel(N_Vector v, int igrp, int ispc, realtype *data);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArray_SpcParallel
 * -----------------------------------------------------------------
 * This function creates an array of SPCPARALLEL vectors by cloning
 * from a given vector w
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArray_SpcParallel(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArrayEmpty_SpcParallel
 * -----------------------------------------------------------------
 * This function creates an array of SPCPARALLEL vectors with empty
 * (NULL) data array by cloning from a given vector w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArrayEmpty_SpcParallel(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VDestroyVectorArray_SpcParallel
 * -----------------------------------------------------------------
 * This function frees an array of N_Vector created with 
 * N_VCloneVectorArray_SpcParallel.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VDestroyVectorArray_SpcParallel(N_Vector *vs, int count);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrint_SpcParallel
 * -----------------------------------------------------------------
 * This function prints the content of an SPCPARALLEL vector 
 * to a file or to stdout (if fname=="")
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrint_SpcParallel(N_Vector v, char *fname);

/*
 * -----------------------------------------------------------------
 * SPCPARALLEL implementations of the vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector    N_VCloneEmpty_SpcParallel(N_Vector);
SUNDIALS_EXPORT N_Vector    N_VClone_SpcParallel(N_Vector);
SUNDIALS_EXPORT void        N_VDestroy_SpcParallel(N_Vector);
SUNDIALS_EXPORT void        N_VSpace_SpcParallel(N_Vector, long int *, long int *);
SUNDIALS_EXPORT realtype   *N_VGetArrayPointer_SpcParallel(N_Vector);
SUNDIALS_EXPORT void        N_VSetArrayPointer_SpcParallel(realtype *, N_Vector);
SUNDIALS_EXPORT void        N_VLinearSum_SpcParallel(realtype, N_Vector, realtype, N_Vector, N_Vector);
SUNDIALS_EXPORT void        N_VConst_SpcParallel(realtype, N_Vector);
SUNDIALS_EXPORT void        N_VProd_SpcParallel(N_Vector, N_Vector, N_Vector);
SUNDIALS_EXPORT void        N_VDiv_SpcParallel(N_Vector, N_Vector, N_Vector);
SUNDIALS_EXPORT void        N_VScale_SpcParallel(realtype, N_Vector, N_Vector);
SUNDIALS_EXPORT void        N_VAbs_SpcParallel(N_Vector, N_Vector);
SUNDIALS_EXPORT void        N_VInv_SpcParallel(N_Vector, N_Vector);
SUNDIALS_EXPORT void        N_VAddConst_SpcParallel(N_Vector, realtype, N_Vector);
SUNDIALS_EXPORT realtype    N_VDotProd_SpcParallel(N_Vector, N_Vector);
SUNDIALS_EXPORT realtype    N_VMaxNorm_SpcParallel(N_Vector);
SUNDIALS_EXPORT realtype    N_VWrmsNorm_SpcParallel(N_Vector, N_Vector);
SUNDIALS_EXPORT realtype    N_VWrmsNormMask_SpcParallel(N_Vector, N_Vector, N_Vector);
SUNDIALS_EXPORT realtype    N_VMin_SpcParallel(N_Vector);
SUNDIALS_EXPORT realtype    N_VWL2Norm_SpcParallel(N_Vector, N_Vector);
SUNDIALS_EXPORT realtype    N_VL1Norm_SpcParallel(N_Vector);
SUNDIALS_EXPORT void        N_VCompare_SpcParallel(realtype, N_Vector, N_Vector);
SUNDIALS_EXPORT booleantype N_VInvTest_SpcParallel(N_Vector, N_Vector);
SUNDIALS_EXPORT booleantype N_VConstrMask_SpcParallel(N_Vector, N_Vector, N_Vector);   
SUNDIALS_EXPORT realtype    N_VMinQuotient_SpcParallel(N_Vector, N_Vector);



#ifdef __cplusplus
}
#endif

#endif

