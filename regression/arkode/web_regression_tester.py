#!/usr/bin/env python
# script to perform regression tests on ARKODE solvers, and generate web page with the results.
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import sys
import os
import numpy as np
from datetime import datetime
import arkode_tools as ark
import random
import string
import shutil as shutil


#### Utility functions ####

def run_tests(outdir,fptr,testlist,ovtol):
    """ This routine takes in a string containing the test      """
    """ category name, a file pointer to write output in, a     """
    """ string containing a set of executable names and a       """
    """ tolerance on the allowable oversolve (larger==better).  """
    """ It then runs the desired tests, and checks whether the  """
    """ tests pass, placing the results into an html table      """
    """ within the designated file.                             """
    import shlex
    import subprocess
    iret = 0;
    # start table
    fptr.write("<table border=\"1\" cellpadding=\"2\" cellspacing=\"0\" style=\"width: 100%; text-align: left;\">\n")
    fptr.write("<tbody>\n")
    fptr.write("<col width=\"100\">\n")
    fptr.write("<tr>\n")
    fptr.write("  <th style=\"vertical-align: top;\">\n")
    fptr.write("    Problem\n")
    fptr.write("  </th>\n")
    fptr.write("  <th style=\"vertical-align: top;\">\n")
    fptr.write("    Pass/Fail\n")
    fptr.write("  </th>\n")
    fptr.write("  <th style=\"vertical-align: top;\">\n")
    fptr.write("    Steps\n")
    fptr.write("  </th>\n")
    fptr.write("  <th style=\"vertical-align: top;\">\n")
    fptr.write("    Oversolve\n")
    fptr.write("  </th>\n")
    fptr.write("  <th style=\"vertical-align: top;\">\n")
    fptr.write("    Runtime\n")
    fptr.write("  </th>\n")
    fptr.write("</tr>\n")
    for i in range(len(testlist)):
        # start table row for this test
        fptr.write("<tr>\n")
        fptr.write("  <td style=\"vertical-align: top;\">\n")
        fptr.write("    %s \n" % (testlist[i]) )
        [fail,nst,ast,nfe,nfi,lset,nfi_lset,nJe,nnewt,ncf,nef,merr,rerr,ov,rt] = ark.run_test(testlist[i],1);
        ofile = outdir + "/output-" + random_hash() + ".txt"
        shutil.copy("output.txt", ofile)
        f.write("&nbsp; (<a href=\"" + ofile + "\">output</a>)\n")
        fptr.write("  </td>")
        fail_integration = 0;
        fail_error = 0
        # check for integration failure
        if (fail == 1):
            fail_integration = 1;
        # check for oversolve >= ovtol (fits within allowable error)
        if ((ov < ovtol) or (ov != ov)):
            fail_error = 1;
        # report on the pass/fail of the test
        fptr.write("  <td style=\"vertical-align: top;\">\n")
        if (fail_integration + fail_error == 0):
            fptr.write("    <span style=\"color:#298A08\">Pass</span>\n")
        else:
            fptr.write("    <span style=\"color:#B40404\">Fail</span>\n")
        fptr.write("  </td>")
        # report the number of steps
        fptr.write("  <td style=\"vertical-align: top;\">\n")
        fptr.write("    %i\n" % (nst))
        fptr.write("  </td>")
        # report the oversolve
        fptr.write("  <td style=\"vertical-align: top;\">\n")
        if (fail_error == 0):
            fptr.write("    <span style=\"color:#298A08\">%g</span>\n" % (ov))
        else:
            fptr.write("    <span style=\"color:#B40404\">%g</span>\n" % (ov))
        fptr.write("  </td>")
        # report the runtime
        fptr.write("  <td style=\"vertical-align: top;\">\n")
        fptr.write("    %g\n" % (rt))
        fptr.write("  </td>")
        fptr.write("</tr>\n")
    # finish off table
    fptr.write("</tbody></table><br>\n")
    return 

def random_hash():
    """ This routine generates a 32-character random hash. """
    pool = string.letters + string.digits
    return ''.join(random.choice(pool) for i in xrange(32))
    
#############


# set up a list of executable names to use in tests
testsI2 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
           './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_brusselator1D.exe', 
           './ark_heat1D.exe', './ark_pollu.exe' )
testsI3 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
           './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_brusselator1D.exe',
           './ark_bruss.exe', './ark_heat1D.exe', './ark_medakzo.exe', './ark_pollu.exe', 
           './ark_rober.exe', './ark_vdpol.exe' )
testsI4 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
           './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe', 
           './ark_brusselator1D.exe', './ark_heat1D.exe', './ark_hires.exe', './ark_medakzo.exe', 
           './ark_orego.exe', './ark_pollu.exe', './ark_ringmod.exe', './ark_rober.exe', 
           './ark_vdpol.exe', './ark_vdpolm.exe' )
testsI5 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
           './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe', 
           './ark_brusselator1D.exe', './ark_heat1D.exe', './ark_hires.exe', './ark_medakzo.exe', 
           './ark_pollu.exe', './ark_rober.exe', './ark_vdpol.exe', './ark_vdpolm.exe' )
testsI = (testsI2, testsI3, testsI4, testsI5)

testsIF2 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
            './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_heat1D.exe' )
testsIF3 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
            './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe', 
            './ark_heat1D.exe', './ark_medakzo.exe' )
testsIF4 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
            './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_heat1D.exe', 
            './ark_medakzo.exe', './ark_ringmod.exe' )
testsIF5 = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
            './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_bruss.exe', 
            './ark_brusselator1D.exe', './ark_heat1D.exe', './ark_hires.exe', './ark_medakzo.exe' )
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
            './ark_bruss.exe', './ark_brusselator1D.exe', './ark_medakzo.exe' )
testsAF4 = ('./ark_analytic.exe', './ark_analytic_sys.exe', './ark_brusselator.exe', 
            './ark_bruss.exe', './ark_brusselator1D.exe', './ark_hires.exe', './ark_medakzo.exe' )
testsAF5 = ('./ark_analytic.exe', './ark_analytic_sys.exe', './ark_brusselator.exe', 
            './ark_bruss.exe', './ark_brusselator1D.exe', './ark_hires.exe', './ark_medakzo.exe' )
testsAF = (testsAF3, testsAF4, testsAF5)

testsE = ('./ark_analytic.exe', './ark_analytic_nonlin.exe', './ark_analytic_nonlin_back.exe', 
          './ark_analytic_sys.exe', './ark_brusselator.exe', './ark_heat1D.exe' )

ovtol  = 0.01;
rtol = (1.e-3, 1.e-6);
atol = (1.e-11, 1.e-11);

# open test results file
now = datetime.today()
s = "./results/regression_results_%.4i_%.2i_%.2i_%.2d_%.2d_%.2d.html" % (now.year, now.month, now.day, now.hour, now.minute, now.second)
f = open(s,'w')

# create test results directory
outdir = "./results/regression_results_%.4i_%.2i_%.2i_%.2d_%.2d_%.2d" % (now.year, now.month, now.day, now.hour, now.minute, now.second)
os.mkdir(outdir)

# print header
f.write("<br>Tests run on %i/%i/%i at %.2i:%.2i:%.2i<br>\n" % (now.month, now.day, now.year, now.hour, now.minute, now.second))
f.write("<br><a href=\"regression_results_old.shtml\">Previous regression test results</a><br>")


# run tests with base set of parameters to ensure everything runs
test_string = "Base tests (defaults, with rtol = %g, atol = %g):" % (rtol[0], atol[0])
p = ark.SolParams(-1.0, -1, 0, -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 
                   0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[0], atol[0]);
ark.write_parameter_file(p);
pfile = outdir + "/solve_params-" + random_hash() + ".txt"
shutil.copy("./solve_params.txt", pfile)
f.write("<br><b>    " + test_string + "</b>  (<a href=\"" + pfile + "\">input parameters</a>)\n")
iret = run_tests(outdir,f,testsI[2],ovtol);

# check ERK method orders {2,3,4,5,6}
ords = (2,3,4,5,6);
for j in range(len(rtol)):
  for i in range(len(ords)):
    test_string = "ERK order %i tests (rtol = %g, atol = %g):" % (ords[i], rtol[j], atol[j])
    p = ark.SolParams(-1.0, -1, ords[i], -1, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    pfile = outdir + "/solve_params-" + random_hash() + ".txt"
    shutil.copy("./solve_params.txt", pfile)
    f.write("<br><b>    " + test_string + "</b>  (<a href=\"" + pfile + "\">input parameters</a>)\n")
    iret = run_tests(outdir,f,testsE,ovtol);

# check DIRK method orders {2,3,4,5}
ords = (2,3,4,5);
for j in range(len(rtol)):
  for i in range(len(ords)):
    test_string = "DIRK order %i tests (rtol = %g, atol = %g):" % (ords[i], rtol[j], atol[j])
    p = ark.SolParams(-1.0, -1, ords[i], -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    pfile = outdir + "/solve_params-" + random_hash() + ".txt"
    shutil.copy("./solve_params.txt", pfile)
    f.write("<br><b>    " + test_string + "</b>  (<a href=\"" + pfile + "\">input parameters</a>)\n")
    iret = run_tests(outdir,f,testsI[i],ovtol);

# check DIRK (fixed-point solver) method orders {2,3,4,5}
ords = (2,3,4,5);
for j in range(len(rtol)):
  for i in range(len(ords)):
    sys.stdout.write("DIRK order %i (fixed point solver) tests (rtol = %g, atol = %g):" 
                     % (ords[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, -1, ords[i], -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 1, 3, 50, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = run_tests(outdir,f,testsIF[i],ovtol);

# check ARK method orders {3,4,5}
ords = (3,4,5);
for j in range(len(rtol)):
  for i in range(len(ords)):
    test_string = "ARK order %i tests (rtol = %g, atol = %g):" % (ords[i], rtol[j], atol[j])
    p = ark.SolParams(-1.0, -1, ords[i], -1, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    pfile = outdir + "/solve_params-" + random_hash() + ".txt"
    shutil.copy("./solve_params.txt", pfile)
    f.write("<br><b>    " + test_string + "</b>  (<a href=\"" + pfile + "\">input parameters</a>)\n")
    iret = run_tests(outdir,f,testsA[i],ovtol);

# check ARK (fixed-point solver) method orders {3,4,5}
ords = (3,4,5);
for j in range(len(rtol)):
  for i in range(len(ords)):
    sys.stdout.write("ARK order %i (fixed point solver) tests (rtol = %g, atol = %g):" 
                     % (ords[i], rtol[j], atol[j]))
    p = ark.SolParams(-1.0, -1, ords[i], -1, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 1, 3, 50, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    iret = run_tests(outdir,f,testsAF[i],ovtol);
    ierr += iret
    itot += len(testsAF[i])

# check time step adaptivity methods {0,1,2,3,4,5} (DIRK only)
algs = (0,1,2,3,4,5);
for j in range(len(rtol)):
  for i in range(len(algs)):
    test_string = "H-adaptivity method %i tests (rtol = %g, atol = %g):" % (algs[i], rtol[j], atol[j])
    p = ark.SolParams(-1.0, -1, 0, -1, 0, algs[i], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    pfile = outdir + "/solve_params-" + random_hash() + ".txt"
    shutil.copy("./solve_params.txt", pfile)
    f.write("<br><b>    " + test_string + "</b>  (<a href=\"" + pfile + "\">input parameters</a>)\n")
    iret = run_tests(outdir,f,testsI[2],ovtol);

# check predictor methods {0,1,2,3} (DIRK only)
algs = (0,1,2,3);
for j in range(len(rtol)):
  for i in range(len(algs)):
    test_string = "Predictor method %i tests (rtol = %g, atol = %g):" % (algs[i], rtol[j], atol[j])
    p = ark.SolParams(-1.0, -1, 0, -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, algs[i], 0, 0, 0, 0, 0.0, rtol[j], atol[j]);
    ark.write_parameter_file(p);
    pfile = outdir + "/solve_params-" + random_hash() + ".txt"
    shutil.copy("./solve_params.txt", pfile)
    f.write("<br><b>    " + test_string + "</b>  (<a href=\"" + pfile + "\">input parameters</a>)\n")
    iret = run_tests(outdir,f,testsI[2],ovtol);

# write footer to page
f.write("Note: oversolve tolerance = %g\n" % (ovtol))
f.write("[\"oversolve\" is defined as tolerance/error (ideally greater than 1.0)]\n")


# close output file
f.close()
    

# copy output file to overwrite previous "new" results
s2 = "%.4i/%.2i/%.2i -- %.2d:%.2d:%.2d" % (now.year, now.month, now.day, now.hour, now.minute, now.second)
shutil.copy(s,"./regression_newresults.html")

# append new file name to "oldresults" list
f2 = open("./regression_oldresults.html","a")
f2.write("<a href=\"" + s + "\">" + s2 + "</a><br>\n")
f2.close()


##### end of script #####
