#!/usr/bin/env python
#------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
#------------------------------------------------------------
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
#------------------------------------------------------------
# script to perform regression tests on ARKODE solvers.

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
        # check if test reported failure
        if (fail == 1):
            tret = 1;
            sys.stdout.write("\n  %30s integration failure" % (testlist[i]))
        # check for oversolve >= ovtol (fits within allowable error)
        if ((ov < ovtol) or (ov != ov)):
            tret = 1;
            sys.stdout.write("\n  %30s failure:  %g error,  %.2g sec" 
                             % (testlist[i], ov, rt))
        if (tret == 0):
            sys.stdout.write("\n  %30s pass:  %6i steps,  %.2e oversolve,  %.2g sec" 
                             % (testlist[i], nst, ov, rt))
        iret += tret;
    if (iret == 0):
#        sys.stdout.write("  pass\n")
        sys.stdout.write("\n")
    else:
        sys.stdout.write("\n")
    return iret
    
#############


# set up a list of all executable names to use in tests
tests = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
         './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_brusselator1D.exe', 
         './ark_bruss.exe', './ark_heat1D.exe', './ark_heat2D.exe', './ark_hires.exe', 
         './ark_medakzo.exe', './ark_orego.exe', './ark_pollu.exe', './ark_ringmod.exe', 
         './ark_rober.exe', './ark_robertson.exe', './ark_vdpol.exe', './ark_vdpolm.exe' )
testsA = ('./ark_analytic.exe', './ark_analytic_sys.exe', './ark_brusselator.exe', 
          './ark_brusselator1D.exe', './ark_bruss.exe', './ark_hires.exe', './ark_medakzo.exe', 
          './ark_orego.exe', './ark_pollu.exe', './ark_ringmod.exe', 
          './ark_rober.exe', './ark_robertson.exe', './ark_vdpol.exe', './ark_vdpolm.exe' )


ovtol = 0.01;
rtol = (1.e-1, 1.e-3, 1.e-5, 1.e-7);
atol = (1.e-11, 1.e-11, 1.e-11, 1.e-11);
nlscoef = 0.001
maxcor = 8

itot = 0
ierr = 0


# check DIRK tables
BTables = (11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)
for j in range(len(rtol)):
  for i in range(len(BTables)):
    sys.stdout.write("DIRK table %i tests (rtol = %g, atol = %g):" % (BTables[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, BTables[i], 0, -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, maxcor, nlscoef, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(tests,ovtol);



# check ARK tables
BTables = (2, 4, 9)
for j in range(len(rtol)):
  for i in range(len(BTables)):
    sys.stdout.write("ARK table %i tests (rtol = %g, atol = %g):" % (BTables[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, BTables[i], 0, -1, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, maxcor, nlscoef, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(testsA,ovtol);



# check ERK tables
BTables = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
for j in range(len(rtol)):
  for i in range(len(BTables)):
    sys.stdout.write("ERK table %i tests (rtol = %g, atol = %g):" % (BTables[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, BTables[i], 0, -1, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = check_tests(tests,ovtol);



##### end of script #####
