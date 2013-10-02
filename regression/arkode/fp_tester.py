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
        # check for test failure
        if (fail == 1):
            tret = 1;
            sys.stdout.write("\n  %30s \033[91m integration failure\033[0m" % (testlist[i]))
        # check for oversolve >= ovtol (fits within allowable error)
        if ((ov < ovtol) or (ov != ov)):
            tret = 1;
            sys.stdout.write("\n  %30s \033[91m failure (too much error: %g < %g)\033[0m" 
                             % (testlist[i], ov, ovtol))
        if (tret == 0):
            sys.stdout.write("\n  %30s: \033[92m %7i nst;  %8i nni;  %.2e ov; \033[94m %.2g sec\033[0m" 
                             % (testlist[i], nst, nnewt, ov, rt))
        iret += tret;
    if (iret == 0):
        sys.stdout.write("\n")
    else:
        sys.stdout.write("\n")
    return iret
    
#############


# set up a list of executable names to use in tests
testsI3 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
           './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe', 
           './ark_brusselator1D.exe', './ark_hires.exe', './ark_medakzo.exe' )
testsI4 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
           './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe', 
           './ark_brusselator1D.exe', './ark_hires.exe', './ark_medakzo.exe' )
testsI5 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
           './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe', 
           './ark_brusselator1D.exe', './ark_hires.exe', './ark_medakzo.exe' )
testsI = (testsI3, testsI4, testsI5)
testsA3 = ('./ark_analytic.exe', './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe',
           './ark_brusselator1D.exe', './ark_hires.exe', './ark_medakzo.exe', './ark_vdpolm.exe' )
testsA4 = ('./ark_analytic.exe', './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe',
           './ark_brusselator1D.exe', './ark_hires.exe', './ark_medakzo.exe', './ark_vdpolm.exe' )
testsA5 = ('./ark_analytic.exe', './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe',
           './ark_brusselator1D.exe', './ark_hires.exe', './ark_medakzo.exe', './ark_vdpolm.exe' )
testsA = (testsA3, testsA4, testsA5)
ovtol  = 0.01;



# check DIRK method orders {3,4,5}
rtol = (1.e-3, 1.e-6);
atol = (1.e-11, 1.e-11);
ords = (3, 4, 5);
maas = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
for k in range(len(maas)):
  for j in range(len(rtol)):
    for i in range(len(ords)):
      sys.stdout.write("DIRK tests, order = %i, rtol = %g, maa = %i:" % (ords[i], rtol[j], maas[k]))
      p = ark.SolParams(-1.0, -1, ords[i], -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                         0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 1, maas[k], 50, 0.0, rtol[j], atol[j]);
      ark.write_parameter_file(p);
      iret = check_tests(testsI[i],ovtol);

      sys.stdout.write("ARK tests, order = %i, rtol = %g, maa = %i:" % (ords[i], rtol[j], maas[k]))
      p = ark.SolParams(-1.0, -1, ords[i], -1, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                         0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 1, maas[k], 50, 0.0, rtol[j], atol[j]);
      ark.write_parameter_file(p);
      iret = check_tests(testsA[i],ovtol);



##### end of script #####
