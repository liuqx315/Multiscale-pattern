#!/usr/bin/env python
# script to perform regression tests on ARKODE solvers.
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import sys
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
            sys.stdout.write("\n  %s failure (too few steps: %i < %i)" % (testlist[i], nst, nsttol))
        # check for oversolve >= ovtol (fits within allowable error)
        if ((ov < ovtol) or (ov != ov)):
            tret = 1;
            sys.stdout.write("\n  %s failure (too much error: %g < %g)" % (testlist[i], ov, ovtol))
        if (tret == 0):
            sys.stdout.write("\n  %s pass (steps: %i > %i;  oversolve %g > %g)" % (testlist[i], nst, nsttol, ov, ovtol))
        iret += tret;
    if (iret == 0):
#        sys.stdout.write("  pass\n")
        sys.stdout.write("\n")
    else:
        sys.stdout.write("\n")
    return iret
    
#############


# set up a list of executable names to use in tests
tests = ('ark_analytic.exe', 'ark_analytic_nonlin.exe', 'ark_analytic_sys.exe', 'ark_brusselator.exe', 'ark_brusselator1D.exe');
tests2 = ('ark_analytic.exe', 'ark_analytic_sys.exe', 'ark_brusselator.exe', 'ark_brusselator1D.exe');
nsttol = 10;
ovtol  = 0.05;
rtol = 1.e-6;
atol = 1.e-10;

itot = 0
ierr = 0

# run tests with base set of parameters to ensure everything runs
sys.stdout.write("Base tests:")
p = ark.SolParams(-1.0, 0, -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                   0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0, rtol, atol);
ark.write_parameter_file(p);
iret = check_tests(tests,nsttol,ovtol);
ierr += iret
itot += len(tests)

# check ERK method orders {2,3,4,5,6}
ords = (2,3,4,5,6);
for i in range(len(ords)):
    sys.stdout.write("ERK order %i tests:" % (ords[i]))
    p = ark.SolParams(-1.0, ords[i], -1, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0, rtol, atol);
    ark.write_parameter_file(p);
    iret = check_tests(tests,nsttol,ovtol);
    ierr += iret
    itot += len(tests)


# check DIRK method orders {3,4,5}
ords = (3,4,5);
for i in range(len(ords)):
    sys.stdout.write("DIRK order %i tests:" % (ords[i]))
    p = ark.SolParams(-1.0, ords[i], -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0, rtol, atol);
    ark.write_parameter_file(p);
    iret = check_tests(tests,nsttol,ovtol);
    ierr += iret
    itot += len(tests)

# check ARK method orders {3,4,5}
ords = (3,4,5);
for i in range(len(ords)):
    sys.stdout.write("ARK order %i tests:" % (ords[i]))
    p = ark.SolParams(-1.0, ords[i], -1, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0, rtol, atol);
    ark.write_parameter_file(p);
    iret = check_tests(tests2,nsttol,ovtol);
    ierr += iret
    itot += len(tests2)

# check time step adaptivity methods {0,1,2,3,4,5} (DIRK only)
algs = (0,1,2,3,4,5);
for i in range(len(algs)):
    sys.stdout.write("H-adaptivity method %i tests:" % (algs[i]))
    p = ark.SolParams(-1.0, 0, -1, 0, algs[i], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0, rtol, atol);
    ark.write_parameter_file(p);
    iret = check_tests(tests,nsttol,ovtol);
    ierr += iret
    itot += len(tests)

# check predictor methods {0,1,2,3} (DIRK only)
algs = (0,1,2,3);
for i in range(len(algs)):
    sys.stdout.write("Predictor method %i tests:" % (algs[i]))
    p = ark.SolParams(-1.0, 0, -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, algs[i], 0, 0, 0.0, rtol, atol);
    ark.write_parameter_file(p);
    iret = check_tests(tests,nsttol,ovtol);
    ierr += iret
    itot += len(tests)


if (ierr == 0):
    sys.stdout.write("\nPassed all %i tests\n" % (itot))
else:
    sys.stdout.write("\nFailed %i out of %i tests\n" % (ierr, itot))
    

##### end of script #####
