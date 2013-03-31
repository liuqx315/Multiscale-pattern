# analysis functions for running ARKODE tests with various parameters
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import sys
import numpy as np
import time

#### Data Structures ####

class SolParams:
    """ Holds a set of solver parameters and its associated cost. """
    def __init__(self, cost, order, dense_order, imex, adapt_method, cflfac, safety, bias, 
                 growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, etamxf, etacf, 
                 small_nef, crdown, rdiv, dgmax, predictor, msbp, maxcor, nlscoef):
        self.cost = cost;
        self.order = order;
        self.dense_order = dense_order;
        self.imex = imex;
        self.adapt_method = adapt_method;
        self.cflfac = cflfac;
        self.safety = safety;
        self.bias = bias;
        self.growth = growth;
        self.hfixed_lb = hfixed_lb;
        self.hfixed_ub = hfixed_ub;
        self.k1 = k1;
        self.k2 = k2;
        self.k3 = k3;
        self.etamx1 = etamx1;
        self.etamxf = etamxf;
        self.etacf = etacf;
        self.small_nef = small_nef;
        self.crdown = crdown;
        self.rdiv = rdiv;
        self.dgmax = dgmax;
        self.predictor = predictor;
        self.msbp = msbp;
        self.maxcor = maxcor;
        self.nlscoef = nlscoef;
    def WriteHeader(self):
        print '   cost    q  dq  imx  adp  cfl   safe  bias  grow  h0l  h0b    k1     k2     k3     emx1     emxf    ecf  smf  crdn   rdiv  dgmx  prd  msbp  mxcr   nlsc'
    def Write(self):
        sys.stdout.write(" %f %2i %3i  %2i  %3i   %3.1f  %5.3f  %4.2f  %4.1f  %3.1f  %3.1f  %5.3f  %5.3f  %5.3f  %6f  %5.2f  %5.2f %3i   %4.2f  %5.3f  %4.2f %3i  %3i   %3i   %.1e\n" % 
                         (self.cost, self.order, self.dense_order, self.imex, self.adapt_method, self.cflfac, 
                          self.safety, self.bias, self.growth, self.hfixed_lb, self.hfixed_ub, 
                          self.k1, self.k2, self.k3, self.etamx1, self.etamxf, self.etacf, 
                          self.small_nef, self.crdown, self.rdiv, self.dgmax, self.predictor,
                          self.msbp, self.maxcor, self.nlscoef))
        #print ' ',self.cost,' ',self.order,' ',self.dense_order,' ',self.imex,' ',self.adapt_method,' ',self.cflfac,' ',self.safety,' ',self.bias,' ',self.growth,' ',self.hfixed_lb,' ',self.hfixed_ub,' ',self.k1,' ',self.k2,' ',self.k3,' ',self.etamx1,' ',self.etamxf,' ',self.etacf,' ',self.small_nef,' ',self.crdown,' ',self.rdiv,' ',self.dgmax,' ',self.predictor,' ',self.msbp,' ',self.maxcor,' ',self.nlscoef


#### Utility functions ####

def run_test(testname,keep_output):
    """ This routine takes in a string containing an executable name, """
    """ runs that test sending all output to a temporary file, and    """
    """ processes that file to determine the overall run statistics,  """
    """ which are returned to the calling routine.                    """
    import shlex
    import os
    import subprocess
    from time import time
    #  run routine
    out = open('output.txt','w')
    null = open(os.devnull, 'w')
    tstart = time()
    subprocess.call(testname, stdout=out, stderr=null)
    runtime = (time() - tstart)
    null.close()
    out.close()
    #  initialize return values to "bad" values, to overwrite with actual data
    nsteps=100000
    asteps=100000
    nfe=100000
    nfi=100000 
    lsetups=100000
    nfi_lsetup=100000
    nJe=100000
    nnewt=100000
    ncf=100000
    nef=100000
    maxerr=1.0e5
    rmserr=1.0e5
    oversolve=1.0e-5
    #  parse output
    f = open('output.txt','r')
    for line in f:
        txt = shlex.split(line)
        if ("Internal" in txt):
            nsteps = int(txt[4]);
            asteps = int(txt[7].replace(')',''));
        elif ("system" in txt):
            nfi_lsetup = int(txt[10]);
        elif ("evals:" in txt):
            nfe = int(txt[5].replace(',',''));
            nfi = int(txt[8]);
        elif ("setups" in txt):
            lsetups = int(txt[5]);
        elif ("Jacobian" in txt):
            nJe = int(txt[6]);
        elif (("Newton" in txt) and ("iterations" in txt)):
            nnewt = int(txt[6]);
        elif (("convergence" in txt) and ("Total" in txt)):
            ncf = int(txt[8]);
        elif (("error" in txt) and ("test" in txt)):
            nef = int(txt[7]);
        elif ("rms" in txt):
            maxerr = float(txt[3].replace(',',''));
            rmserr = float(txt[6]);
        elif ("Oversolve" in txt):
            oversolve = float(txt[2]);
    f.close()
    if (keep_output==0):
        os.remove('output.txt')
    return [nsteps, asteps, nfe, nfi, lsetups, nfi_lsetup, nJe, nnewt, ncf, nef, maxerr, rmserr, oversolve, runtime]
    
##########
def write_parameter_file(params):
    """ This routine takes in a set of solver parameterrs and writes """
    """ an input file for the tests to supply these parameters to    """
    """ their solvers.                                               """
    f = open('solve_params.txt', 'w')
    f.write("order = %i\n" % (params.order)) 
    f.write("dense_order = %i\n" % (params.dense_order)) 
    f.write("imex = %i\n" % (params.imex)) 
    f.write("btable = -1\n") 
    f.write("adapt_method = %i\n" % (params.adapt_method)) 
    f.write("cflfac = %f\n" % (params.cflfac)) 
    f.write("safety = %f\n" % (params.safety)) 
    f.write("bias = %f\n" % (params.bias)) 
    f.write("growth = %f\n" % (params.growth)) 
    f.write("hfixed_lb = %f\n" % (params.hfixed_lb)) 
    f.write("hfixed_ub = %f\n" % (params.hfixed_ub)) 
    f.write("k1 = %f\n" % (params.k1)) 
    f.write("k2 = %f\n" % (params.k2)) 
    f.write("k3 = %f\n" % (params.k3))
    f.write("etamx1 = %f\n" % (params.etamx1)) 
    f.write("etamxf = %f\n" % (params.etamxf)) 
    f.write("etacf = %f\n" % (params.etacf)) 
    f.write("small_nef = %i\n" % (params.small_nef)) 
    f.write("crdown = %f\n" % (params.crdown)) 
    f.write("rdiv = %f\n" % (params.rdiv)) 
    f.write("dgmax = %f\n" % (params.dgmax)) 
    f.write("predictor = %i\n" % (params.predictor)) 
    f.write("msbp = %i\n" % (params.msbp)) 
    f.write("maxcor = %i\n" % (params.maxcor)) 
    f.write("nlscoef = %f\n" % (params.nlscoef)) 
    f.close()
    

##### end of script #####
