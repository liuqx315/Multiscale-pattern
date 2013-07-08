#!/usr/bin/env python
# script to perform regression tests on ARKODE solvers.
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import sys
import time
import numpy as np
import arkode_tools as ark

#### Utility functions ####

def check_tests(testlist,nsttol,ovtol):
    """ This routine takes in a string containing a set of      """
    """ executable names, a tolerance on the minimum number of  """
    """ steps that must be run to consider the test 'completed' """
    """ and a tolerance on the allowable oversolve (larger is   """
    """ better).  It then runs the desired tests, and checks    """
    """ whether the tests pass.                                 """
    import shlex
    import os
    import subprocess
    iret = 0;
    for i in range(len(testlist)):
        tret = 0
        [nst,ast,nfe,nfi,lset,nfi_lset,nJe,nnewt,ncf,nef,merr,rerr,ov,rt] = ark.run_test(testlist[i],0);
        # check for nst >= nsttol (in case something fails at initialization)
        if (nst < nsttol):
            tret = 1;
            sys.stdout.write("\n  %s \033[91m failure (too few steps: %i < %i) \033[94m [%.2g s]\033[0m" % (testlist[i], nst, nsttol, rt))
        # check for oversolve >= ovtol (fits within allowable error)
        if ((ov < ovtol) or (ov != ov)):
            tret = 1;
            sys.stdout.write("\n  %s \033[91m failure (too much error: %g < %g)\033[94m [%.2g s]\033[0m" % (testlist[i], ov, ovtol, rt))
        if (tret == 0):
            sys.stdout.write("\n  %s \033[92m pass (steps: %i > %i;  oversolve %g > %g)\033[94m [%.2g s]\033[0m" % (testlist[i], nst, nsttol, ov, ovtol, rt))
        iret += tret;
    if (iret == 0):
#        sys.stdout.write("  pass\n")
        sys.stdout.write("\n")
    else:
        sys.stdout.write("\n")
    return iret
    
#############


# set up a list of executable names to use in tests
testsI3 = ('ark_analytic.exe', 'ark_analytic_nonlin.exe', 'ark_analytic_nonlin_back.exe', 
           'ark_analytic_sys.exe', 'ark_brusselator.exe', 'ark_brusselator1D.exe',
           'ark_medakzo.exe' )
testsI4 = ('ark_analytic.exe', 'ark_analytic_nonlin.exe', 'ark_analytic_nonlin_back.exe', 
           'ark_analytic_sys.exe', 'ark_brusselator.exe', 'ark_bruss.exe', 
           'ark_brusselator1D.exe', 'ark_hires.exe', 'ark_medakzo.exe', 'ark_orego.exe',
           'ark_pollu.exe', 'ark_ringmod.exe', 'ark_rober.exe', 'ark_vdpol.exe',
           'ark_vdpolm.exe' )
testsI5 = ('ark_analytic.exe', 'ark_analytic_nonlin.exe', 'ark_analytic_nonlin_back.exe', 
           'ark_analytic_sys.exe', 'ark_brusselator.exe', 'ark_bruss.exe', 
           'ark_brusselator1D.exe', 'ark_medakzo.exe', 'ark_orego.exe',
           'ark_pollu.exe', 'ark_rober.exe', 'ark_vdpol.exe' )
testsI = (testsI3, testsI4, testsI5)
testsIFP = ('ark_analytic.exe', 'ark_analytic_nonlin.exe', 'ark_analytic_nonlin_back.exe', 
            'ark_analytic_sys.exe', 'ark_brusselator.exe' )
testsE = ('ark_analytic.exe', 'ark_analytic_nonlin.exe', 'ark_analytic_nonlin_back.exe', 
          'ark_analytic_sys.exe', 'ark_brusselator.exe' )
testsA3 = ('ark_analytic.exe', 'ark_analytic_sys.exe', 'ark_brusselator.exe', 'ark_bruss.exe',
           'ark_brusselator1D.exe', 'ark_medakzo.exe', 'ark_pollu.exe', 'ark_vdpol.exe', 
           'ark_vdpolm.exe' )
testsA4 = ('ark_analytic.exe', 'ark_analytic_sys.exe', 'ark_brusselator.exe', 'ark_bruss.exe',
           'ark_brusselator1D.exe', 'ark_hires.exe', 'ark_medakzo.exe', 'ark_pollu.exe',
           'ark_vdpol.exe' )
testsA5 = ('ark_analytic.exe', 'ark_analytic_sys.exe', 'ark_brusselator.exe', 'ark_bruss.exe',
           'ark_brusselator1D.exe', 'ark_hires.exe', 'ark_medakzo.exe', 'ark_pollu.exe',
           'ark_vdpol.exe', 'ark_vdpolm.exe' )
testsA = (testsA3, testsA4, testsA5)
nsttol = 10;
ovtol  = 0.01;
rtol = (1.e-3, 1.e-6);
atol = (1.e-11, 1.e-11);

itot = 0
ierr = 0

# run tests with base set of parameters to ensure everything runs
sys.stdout.write("Base tests (rtol = %g, atol = %g):" % (rtol[0], atol[0]))
p = ark.SolParams(-1.0, 0, -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 
                   0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[0], atol[0]);
ark.write_parameter_file(p);
iret = check_tests(testsI[1],nsttol,ovtol);
ierr += iret
itot += len(testsI[1])

# check ERK method orders {2,3,4,5,6}
ords = (2,3,4,5,6);
for j in range(len(rtol)):
  for i in range(len(ords)):
    sys.stdout.write("ERK order %i tests (rtol = %g, atol = %g):" % (ords[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, ords[i], -1, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(testsE,nsttol,ovtol);
    ierr += iret
    itot += len(testsE)


# check DIRK method orders {3,4,5}
ords = (3,4,5);
for j in range(len(rtol)):
  for i in range(len(ords)):
    sys.stdout.write("DIRK order %i tests (rtol = %g, atol = %g):" % (ords[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, ords[i], -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(testsI[i],nsttol,ovtol);
    ierr += iret
    itot += len(testsI[i])

# # check DIRK (fixed-point solver) method orders {3,4,5}
# ords = (3,4,5);
# for j in range(len(rtol)):
#   for i in range(len(ords)):
#     sys.stdout.write("DIRK order %i (fixed point solver) tests (rtol = %g, atol = %g):" 
#                      % (ords[i], rtol[j], atol[j]))
#     p = ark.SolParams(-1.0, ords[i], -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
#                        0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 1, 10, 50, 0.0, rtol[j], atol[j]);
#     ark.write_parameter_file(p);
#     iret = check_tests(testsIFP,nsttol,ovtol);
#     ierr += iret
#     itot += len(testsIFP)

# check ARK method orders {3,4,5}
ords = (3,4,5);
for j in range(len(rtol)):
  for i in range(len(ords)):
    sys.stdout.write("ARK order %i tests (rtol = %g, atol = %g):" % (ords[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, ords[i], -1, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(testsA[i],nsttol,ovtol);
    ierr += iret
    itot += len(testsA[i])

# check time step adaptivity methods {0,1,2,3,4,5} (DIRK only)
algs = (0,2,3);
for j in range(len(rtol)):
  for i in range(len(algs)):
    sys.stdout.write("H-adaptivity method %i tests (rtol = %g, atol = %g):" 
                     % (algs[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, 0, -1, 0, algs[i], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(testsI[1],nsttol,ovtol);
    ierr += iret
    itot += len(testsI[1])

# check predictor methods {0,1,2,3} (DIRK only)
algs = (0,2);
for j in range(len(rtol)):
  for i in range(len(algs)):
    sys.stdout.write("Predictor method %i tests (rtol = %g, atol = %g):" 
                     % (algs[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, 0, -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, algs[i], 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(testsI[1],nsttol,ovtol);
    ierr += iret
    itot += len(testsI[1])


if (ierr == 0):
    sys.stdout.write("\nPassed all %i tests\n" % (itot))
else:
    sys.stdout.write("\nFailed %i out of %i tests\n" % (ierr, itot))
    

##### end of script #####
