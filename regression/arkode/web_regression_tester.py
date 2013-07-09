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

def run_tests(outdir,fptr,testlist,nsttol,ovtol):
    """ This routine takes in a string containing the test      """
    """ category name, a file pointer to write output in, a     """
    """ string containing a set of executable names, a          """
    """ tolerance on the minimum number of steps that must be   """
    """ run to consider the test 'completed' and a tolerance on """
    """ the allowable oversolve (larger==better).  It then runs """
    """ the desired tests, and checks whether the tests pass,   """
    """ placing the results into an html table within the       """
    """ designated file.                                        """
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
        [nst,ast,nfe,nfi,lset,nfi_lset,nJe,nnewt,ncf,nef,merr,rerr,ov,rt] = ark.run_test(testlist[i],1);
        ofile = outdir + "/output-" + random_hash() + ".txt"
        shutil.copy("output.txt", ofile)
        f.write("&nbsp; (<a href=\"" + ofile + "\">output</a>)\n")
        fptr.write("  </td>")
        fail_steps = 0
        fail_error = 0
        # check for nst >= nsttol (in case something fails at initialization)
        if (nst < nsttol):
            fail_steps = 1;
        # check for oversolve >= ovtol (fits within allowable error)
        if ((ov < ovtol) or (ov != ov)):
            fail_error = 1;
        # report on the pass/fail of the test
        fptr.write("  <td style=\"vertical-align: top;\">\n")
        if (fail_steps + fail_error == 0):
            fptr.write("    <span style=\"color:#298A08\">Pass</span>\n")
        else:
            fptr.write("    <span style=\"color:#B40404\">Fail</span>\n")
        fptr.write("  </td>")
        # report the number of steps
        fptr.write("  <td style=\"vertical-align: top;\">\n")
        if (fail_steps == 0):
            fptr.write("    <span style=\"color:#298A08\">%i</span>\n" % (nst))
        else:
            fptr.write("    <span style=\"color:#B40404\">%i</span>\n" % (nst))
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
testsI3 = ('ark_analytic.exe', 'ark_analytic_nonlin.exe', 'ark_analytic_nonlin_back.exe', 
           'ark_analytic_sys.exe', 'ark_brusselator.exe', 'ark_brusselator1D.exe',
           'ark_medakzo.exe', 'ark_rober.exe' )
testsI4 = ('ark_analytic.exe', 'ark_analytic_nonlin.exe', 'ark_analytic_nonlin_back.exe', 
           'ark_analytic_sys.exe', 'ark_brusselator.exe', 'ark_bruss.exe', 
           'ark_brusselator1D.exe', 'ark_hires.exe', 'ark_medakzo.exe', 'ark_orego.exe',
           'ark_pollu.exe', 'ark_ringmod.exe', 'ark_rober.exe', 'ark_vdpol.exe',
           'ark_vdpolm.exe' )
testsI5 = ('ark_analytic.exe', 'ark_analytic_nonlin.exe', 'ark_analytic_nonlin_back.exe', 
           'ark_analytic_sys.exe', 'ark_brusselator.exe', 'ark_bruss.exe', 
           'ark_brusselator1D.exe', 'ark_hires.exe', 'ark_medakzo.exe', 'ark_orego.exe',
           'ark_pollu.exe', 'ark_rober.exe', 'ark_vdpol.exe' )
testsI = (testsI3, testsI4, testsI5)
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
p = ark.SolParams(-1.0, 0, -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 
                   0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol[0], atol[0]);
ark.write_parameter_file(p);
pfile = outdir + "/solve_params-" + random_hash() + ".txt"
shutil.copy("./solve_params.txt", pfile)
f.write("<br><b>    " + test_string + "</b>  (<a href=\"" + pfile + "\">input parameters</a>)\n")
iret = run_tests(outdir,f,testsI[1],nsttol,ovtol);

# check ERK method orders {2,3,4,5,6}
ords = (2,3,4,5,6);
for j in range(len(rtol)):
  for i in range(len(ords)):
    test_string = "ERK order %i tests (rtol = %g, atol = %g):" % (ords[i], rtol[j], atol[j])
    p = ark.SolParams(-1.0, ords[i], -1, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol, atol);
    ark.write_parameter_file(p);
    pfile = outdir + "/solve_params-" + random_hash() + ".txt"
    shutil.copy("./solve_params.txt", pfile)
    f.write("<br><b>    " + test_string + "</b>  (<a href=\"" + pfile + "\">input parameters</a>)\n")
    iret = run_tests(outdir,f,testsE,nsttol,ovtol);

# check DIRK method orders {3,4,5}
ords = (3,4,5);
for j in range(len(rtol)):
  for i in range(len(ords)):
    test_string = "DIRK order %i tests (rtol = %g, atol = %g):" % (ords[i], rtol[j], atol[j])
    p = ark.SolParams(-1.0, ords[i], -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol, atol);
    ark.write_parameter_file(p);
    pfile = outdir + "/solve_params-" + random_hash() + ".txt"
    shutil.copy("./solve_params.txt", pfile)
    f.write("<br><b>    " + test_string + "</b>  (<a href=\"" + pfile + "\">input parameters</a>)\n")
    iret = run_tests(outdir,f,testsI[i],nsttol,ovtol);

# check ARK method orders {3,4,5}
ords = (3,4,5);
for j in range(len(rtol)):
  for i in range(len(ords)):
    test_string = "ARK order %i tests (rtol = %g, atol = %g):" % (ords[i], rtol[j], atol[j])
    p = ark.SolParams(-1.0, ords[i], -1, 2, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol, atol);
    ark.write_parameter_file(p);
    pfile = outdir + "/solve_params-" + random_hash() + ".txt"
    shutil.copy("./solve_params.txt", pfile)
    f.write("<br><b>    " + test_string + "</b>  (<a href=\"" + pfile + "\">input parameters</a>)\n")
    iret = run_tests(outdir,f,testsA[i],nsttol,ovtol);

# check time step adaptivity methods {0,1,2,3,4,5} (DIRK only)
algs = (0,2,3);
for j in range(len(rtol)):
  for i in range(len(algs)):
    test_string = "H-adaptivity method %i tests (rtol = %g, atol = %g):" % (algs[i], rtol[j], atol[j])
    p = ark.SolParams(-1.0, 0, -1, 0, algs[i], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0.0, rtol, atol);
    ark.write_parameter_file(p);
    pfile = outdir + "/solve_params-" + random_hash() + ".txt"
    shutil.copy("./solve_params.txt", pfile)
    f.write("<br><b>    " + test_string + "</b>  (<a href=\"" + pfile + "\">input parameters</a>)\n")
    iret = run_tests(outdir,f,testsI[1],nsttol,ovtol);

# check predictor methods {0,1,2,3} (DIRK only)
algs = (0,2);
for j in range(len(rtol)):
  for i in range(len(algs)):
    test_string = "Predictor method %i tests (rtol = %g, atol = %g):" % (algs[i], rtol[j], atol[j])
    p = ark.SolParams(-1.0, 0, -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, algs[i], 0, 0, 0, 0, 0.0, rtol, atol);
    ark.write_parameter_file(p);
    pfile = outdir + "/solve_params-" + random_hash() + ".txt"
    shutil.copy("./solve_params.txt", pfile)
    f.write("<br><b>    " + test_string + "</b>  (<a href=\"" + pfile + "\">input parameters</a>)\n")
    iret = run_tests(outdir,f,testsI[1],nsttol,ovtol);

# write footer to page
f.write("Note: step tolerance = %i; oversolve tolerance = %g\n" % (nsttol, ovtol))
f.write("[\"oversolve\" is defined as tolerance/error (ideally greater than 1.0)\n")


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
