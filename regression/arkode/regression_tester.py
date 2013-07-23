#!/usr/bin/env python
# script to perform regression tests on ARKODE solvers.
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import sys
import time
import numpy as np
import arkode_tools as ark

#### Utility functions ####

def check_tests(testlist,ovtol):
    """ This routine takes in a string containing a set of      """
    """ executable names and a tolerance on the allowable       """
    """ oversolve (larger is better).  It then runs the desired """
    """ tests, and checks whether the tests pass.               """
    import shlex
    import os
    import subprocess
    iret = 0;
    for i in range(len(testlist)):
        tret = 0
        [fail,nst,ast,nfe,nfi,lset,nfi_lset,nJe,nnewt,ncf,nef,merr,rerr,ov,rt] = ark.run_test(testlist[i],0); 
        # check for integration failure
        if (fail == 1):
            tret = 1;
            sys.stdout.write("\n  %s \033[91m integration failure\033[0m" % (testlist[i]))
        # check for oversolve >= ovtol (fits within allowable error)
        if ((ov < ovtol) or (ov != ov)):
            tret = 1;
            sys.stdout.write("\n  %s \033[91m failure (too much error: %g < %g)\033[94m [%.2g s]\033[0m" % (testlist[i], ov, ovtol, rt))
        if (tret == 0):
            sys.stdout.write("\n  %s \033[92m pass (steps: %i;  oversolve %g)\033[94m [%.2g s]\033[0m" % (testlist[i], nst, ov, rt))
        iret += tret;
    if (iret == 0):
#        sys.stdout.write("  pass\n")
        sys.stdout.write("\n")
    else:
        sys.stdout.write("\n")
    return iret
    
#############


# set up a list of executable names to use in tests
testsI2 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
           './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_brusselator1D.exe', 
           './ark_heat1D.exe', './ark_pollu.exe' )
testsI3 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
           './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_brusselator1D.exe',
           './ark_bruss.exe', './ark_heat1D.exe', './ark_hires.exe', './ark_medakzo.exe', 
           './ark_pollu.exe', './ark_vdpol.exe' )
testsI4 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
           './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe', 
           './ark_brusselator1D.exe', './ark_heat1D.exe', './ark_hires.exe', './ark_medakzo.exe', 
           './ark_orego.exe', './ark_pollu.exe', './ark_ringmod.exe', './ark_rober.exe', 
           './ark_robertson.exe', './ark_vdpol.exe', './ark_vdpolm.exe' )
testsI5 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
           './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe', 
           './ark_brusselator1D.exe', './ark_heat1D.exe', './ark_hires.exe', './ark_medakzo.exe', 
           './ark_pollu.exe', './ark_rober.exe', './ark_robertson.exe', './ark_vdpol.exe', 
           './ark_vdpolm.exe' )
testsI = (testsI2, testsI3, testsI4, testsI5)

testsIF2 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
            './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_heat1D.exe' )
testsIF3 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
            './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe', 
            './ark_heat1D.exe' )
testsIF4 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
            './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_heat1D.exe' )
testsIF5 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
            './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe', 
            './ark_heat1D.exe' )
testsIF = (testsIF2, testsIF3, testsIF4, testsIF5)

testsA3 = ('./ark_analytic.exe', './ark_analytic_sys.exe', './ark_brusselator.exe', 
           './ark_bruss.exe', './ark_brusselator1D.exe', './ark_medakzo.exe', 
           './ark_pollu.exe', './ark_vdpol.exe' )
testsA4 = ('./ark_analytic.exe', './ark_analytic_sys.exe', './ark_brusselator.exe', 
           './ark_bruss.exe', './ark_brusselator1D.exe', './ark_hires.exe', './ark_medakzo.exe', 
           './ark_pollu.exe', './ark_vdpol.exe' )
testsA5 = ('./ark_analytic.exe', './ark_analytic_sys.exe', './ark_brusselator.exe', 
           './ark_bruss.exe', './ark_brusselator1D.exe', './ark_hires.exe', './ark_medakzo.exe', 
           './ark_pollu.exe', './ark_vdpol.exe' )
testsA = (testsA3, testsA4, testsA5)

testsAF3 = ('./ark_analytic.exe', './ark_analytic_sys.exe', './ark_brusselator.exe', 
            './ark_bruss.exe', './ark_brusselator1D.exe' )
testsAF4 = ('./ark_analytic.exe', './ark_analytic_sys.exe', './ark_brusselator.exe', 
            './ark_bruss.exe', './ark_brusselator1D.exe', './ark_hires.exe' )
testsAF5 = ('./ark_analytic.exe', './ark_analytic_sys.exe', './ark_brusselator.exe', 
            './ark_bruss.exe', './ark_brusselator1D.exe', './ark_hires.exe' )
testsAF = (testsAF3, testsAF4, testsAF5)

testsE = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
          './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_heat1D.exe' )

ovtol = 0.01;
rtol = (1.e-3, 1.e-6);
atol = (1.e-11, 1.e-11);

itot = 0
ierr = 0

# run tests with base set of parameters to ensure everything runs
sys.stdout.write("Base tests (rtol = %g, atol = %g):" % (rtol[0], atol[0]))
p = ark.SolParams(-1.0, -1, 0, -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 
                   0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[0], atol[0]);
ark.write_parameter_file(p);
iret = check_tests(testsI[2],ovtol);
ierr += iret
itot += len(testsI[2])

# check ERK method orders {2,3,4,5,6}
ords = (2,3,4,5,6);
for j in range(len(rtol)):
  for i in range(len(ords)):
    sys.stdout.write("ERK order %i tests (rtol = %g, atol = %g):" % (ords[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, -1, ords[i], -1, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(testsE,ovtol);
    ierr += iret
    itot += len(testsE)


# check DIRK method orders {2,3,4,5}
ords = (2,3,4,5);
for j in range(len(rtol)):
  for i in range(len(ords)):
    sys.stdout.write("DIRK order %i tests (rtol = %g, atol = %g):" % (ords[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, -1, ords[i], -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(testsI[i],ovtol);
    ierr += iret
    itot += len(testsI[i])

# check DIRK (fixed-point solver) method orders {2,3,4,5}
ords = (2,3,4,5);
for j in range(len(rtol)):
  for i in range(len(ords)):
    sys.stdout.write("DIRK order %i (fixed point solver) tests (rtol = %g, atol = %g):" 
                     % (ords[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, -1, ords[i], -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 1, 3, 50, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(testsIF[i],ovtol);
    ierr += iret
    itot += len(testsIF[i])

# check ARK method orders {3,4,5}
ords = (3,4,5);
for j in range(len(rtol)):
  for i in range(len(ords)):
    sys.stdout.write("ARK order %i tests (rtol = %g, atol = %g):" % (ords[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, -1, ords[i], -1, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(testsA[i],ovtol);
    ierr += iret
    itot += len(testsA[i])

# check ARK (fixed-point solver) method orders {3,4,5}
ords = (3,4,5);
for j in range(len(rtol)):
  for i in range(len(ords)):
    sys.stdout.write("ARK order %i (fixed point solver) tests (rtol = %g, atol = %g):" 
                     % (ords[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, -1, ords[i], -1, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 1, 3, 50, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(testsAF[i],ovtol);
    ierr += iret
    itot += len(testsAF[i])

# check time step adaptivity methods {0,1,2,3,4,5} (DIRK only)
algs = (0,1,2,3,4,5);
for j in range(len(rtol)):
  for i in range(len(algs)):
    sys.stdout.write("H-adaptivity method %i tests (rtol = %g, atol = %g):" 
                     % (algs[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, -1, 0, -1, 0, algs[i], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(testsI[2],ovtol);
    ierr += iret
    itot += len(testsI[2])

# check predictor methods {0,1,2,3} (DIRK only)
algs = (0,1,2,3);
for j in range(len(rtol)):
  for i in range(len(algs)):
    sys.stdout.write("Predictor method %i tests (rtol = %g, atol = %g):" 
                     % (algs[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, -1, 0, -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, algs[i], 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(testsI[2],ovtol);
    ierr += iret
    itot += len(testsI[2])


if (ierr == 0):
    sys.stdout.write("\nPassed all %i tests\n" % (itot))
else:
    sys.stdout.write("\nFailed %i out of %i tests\n" % (ierr, itot))
    

##### end of script #####
