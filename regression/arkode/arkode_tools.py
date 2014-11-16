#------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
#------------------------------------------------------------
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
#------------------------------------------------------------
# analysis functions for running ARKODE tests with various parameters

# imports
import sys
import numpy as np
import time

#### Data Structures ####

class SolParams:
    """ Holds a set of solver parameters and its associated cost. """
    def __init__(self, cost=-1.0, btable=-1, order=0, dense_order=-1, imex=0, 
                 adapt_method=0, cflfac=0.0, safety=0.0, bias=0.0, growth=0.0, 
                 hfixed_lb=0.0, hfixed_ub=0.0, pq=0, k1=0.0, k2=0.0, k3=0.0, 
                 etamx1=0.0, etamxf=0.0, etacf=0.0, small_nef=0, crdown=0.0, 
                 rdiv=0.0, dgmax=0.0, predictor=0, msbp=0, fixedpt=0, m_aa=0, 
                 maxcor=0, nlscoef=0.0, rtol=1.e-3, atol=1.e-11):
        self.cost = cost;
        self.btable = btable;
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
        self.pq = pq;
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
        self.fixedpt = fixedpt;
        self.m_aa = m_aa;
        self.maxcor = maxcor;
        self.nlscoef = nlscoef;
        self.rtol = rtol;
        self.atol = atol;
    def WriteHeader(self,fstream=sys.stdout):
        fstream.write("   cost  tab q  dq  imx  adp  cfl   safe  bias  grow  h0l h0b  pq     k1     k2     k3     emx1     emxf    ecf  smf  crdn   rdiv  dgmx  prd  msbp  fxpt  m_aa  mxcr  nlsc    rtol  atol\n")
    def Write(self,fstream=sys.stdout):
        fstream.write(" %7.5f %2i %2i %3i  %2i  %3i   %3.1f  %5.3f  %4.2f  %4.1f  %3.1f  %3.1f  %2i  %5.3f  %5.3f  %5.3f  %6f  %5.2f  %5.2f %3i   %4.2f  %5.3f  %4.2f %3i  %3i   %3i   %3i   %3i   %.1e  %3.1e  %3.1e\n" % 
                         (self.cost, self.btable, self.order, self.dense_order, self.imex, 
                          self.adapt_method, self.cflfac, self.safety, self.bias, self.growth, 
                          self.hfixed_lb, self.hfixed_ub,  self.pq, self.k1, self.k2, self.k3, 
                          self.etamx1, self.etamxf, self.etacf, self.small_nef, self.crdown, 
                          self.rdiv, self.dgmax, self.predictor, self.msbp, self.fixedpt, 
                          self.m_aa, self.maxcor, self.nlscoef, self.rtol, self.atol))


class Test:
    """ Holds the name of an executable, as well as its associated """
    """ set of solver parameters. """
    def __init__(self, name='', executable='', params='', line=''):
        if (line=='' and (name=='' or executable=='' or params=='')):
            print 'Test constructor error: must supply either a line or contents'
        elif (line != ''):
            import shlex
            txt = shlex.split(line)
            if ("Test" in txt):
                self.name = txt[1]
                self.executable = txt[2]
                btable = int(txt[4])
                order = int(txt[5])
                dense_order = int(txt[6])
                imex = int(txt[7])
                adapt_method = int(txt[8])
                cflfac = float(txt[9])
                safety = float(txt[10])
                bias = float(txt[11])
                growth = float(txt[12])
                hfixed_lb = float(txt[13])
                hfixed_ub = float(txt[14])
                pq = int(txt[15])
                k1 = float(txt[16])
                k2 = float(txt[17])
                k3 = float(txt[18])
                etamx1 = float(txt[19])
                etamxf = float(txt[20])
                etacf = float(txt[21])
                small_nef = int(txt[22])
                crdown = float(txt[23])
                rdiv = float(txt[24])
                dgmax = float(txt[25])
                predictor = int(txt[26])
                msbp = int(txt[27])
                fixedpt = int(txt[28])
                m_aa = int(txt[29])
                maxcor = int(txt[30])
                nlscoef = float(txt[31])
                rtol = float(txt[32])
                atol = float(txt[33])
                self.params = SolParams(-1.0, btable, order, dense_order, imex, 
                                         adapt_method, cflfac, safety, bias, 
                                         growth, hfixed_lb, hfixed_ub, pq, k1, k2,
                                         k3, etamx1, etamxf, etacf, small_nef, 
                                         crdown, rdiv, dgmax, predictor, msbp, 
                                         fixedpt, m_aa, maxcor, nlscoef, rtol, atol)
            else:
                print 'Test constructor error: incompatible text line'
        else:
            self.name = name
            self.executable = executable
            self.params = params
    def Write(self, fstream=sys.stdout):
        fstream.write("Test %s %s" % (self.name, self.executable) )
        self.params.Write(fstream);
    def Run(self,keep_output):
        #  output parameter file
        write_parameter_file(self.params)
        #  call test-runner
        stats = run_test(self.executable, self.params)
        return stats


class RunStats:
    """ Holds the run statistics resulting from a given test problem. """
    def __init__(self):
        self.errfail=0
        self.nsteps=0
        self.asteps=0
        self.nfe=0 
        self.nfi=0
        self.lsetups=0
        self.nfi_lsetup=0
        self.nJe=0
        self.nnewt=0
        self.ncf=0
        self.nef=0
        self.maxerr=0
        self.rmserr=0
        self.oversolve=0
        self.runtime=0
    def Compare(self, other, fstream=sys.stdout, strict=False):
        same = 1
        if (strict):
            tol = 1.0
            rttol = 1.5
        else:
            tol = 1.05
            rttol = 2.0
        if (self.errfail != other.errfail):
            same = 0
            fstream.write("\n    Compare: errfail differs (%i vs %i)" % (self.errfail, other.errfail))
        if ((self.nsteps < other.nsteps/tol) or (self.nsteps > other.nsteps*tol)):
            same = 0
            fstream.write("\n    Compare: nsteps differs (%i vs %i)" % (self.nsteps, other.nsteps))
        if ((self.asteps < other.asteps/tol) or (self.asteps > other.asteps*tol)):
            same = 0
            fstream.write("\n    Compare: asteps differs (%i vs %i)" % (self.asteps, other.asteps))
        if ((self.nfe < other.nfe/tol) or (self.nfe > other.nfe*tol)):
            same = 0
            fstream.write("\n    Compare: nfe differs (%i vs %i)" % (self.nfe, other.nfe))
        if ((self.nfi < other.nfi/tol) or (self.nfi > other.nfi*tol)):
            same = 0
            fstream.write("\n    Compare: nfi differs (%i vs %i)" % (self.nfi, other.nfi))
        if ((self.lsetups < other.lsetups/tol) or (self.lsetups > other.lsetups*tol)):
            same = 0
            fstream.write("\n    Compare: lsetups differs (%i vs %i)" % (self.lsetups, other.lsetups))
        if ((self.nfi_lsetup < other.nfi_lsetup/tol) or (self.nfi_lsetup > other.nfi_lsetup*tol)):
            same = 0
            fstream.write("\n    Compare: nfi_lsetup differs (%i vs %i)" % (self.nfi_lsetup, other.nfi_lsetup))
        if ((self.nJe < other.nJe/tol) or (self.nJe > other.nJe*tol)):
            same = 0
            fstream.write("\n    Compare: nJe differs (%i vs %i)" % (self.nJe, other.nJe))
        if ((self.nnewt < other.nnewt/tol) or (self.nnewt > other.nnewt*tol)):
            same = 0
            fstream.write("\n    Compare: nnewt differs (%i vs %i)" % (self.nnewt, other.nnewt))
        if ((self.ncf < other.ncf/tol) or (self.ncf > other.ncf*tol)):
            same = 0
            fstream.write("\n    Compare: ncf differs (%i vs %i)" % (self.ncf, other.ncf))
        if ((self.nef < other.nef/tol) or (self.nef > other.nef*tol)):
            same = 0
            fstream.write("\n    Compare: nef differs (%i vs %i)" % (self.nef, other.nef))
        if ((self.maxerr < other.maxerr/tol) or (self.maxerr > other.maxerr*tol)):
            same = 0
            fstream.write("\n    Compare: maxerr differs (%g vs %g)" % (self.maxerr, other.maxerr))
        if ((self.rmserr < other.rmserr/tol) or (self.rmserr > other.rmserr*tol)):
            same = 0
            fstream.write("\n    Compare: rmserr differs (%g vs %g)" % (self.rmserr, other.rmserr))
        if ((self.runtime < other.runtime/rttol) or (self.runtime > other.runtime*rttol)):
            if (same == 1):
                same = 2
            fstream.write("\n    Compare: runtime differs (%g vs %g)" % (self.runtime, other.runtime))
        if (same != 1):
            fstream.write("\n")
        return same
    def Write(self, fstream=sys.stdout):
        fstream.write("Run stats:\n")
        if (errfail != 0):
            fstream.write("  errfail = %i\n" % (self.errfail) )
        if (nsteps != 0):
            fstream.write("  nsteps = %i\n" % (self.nsteps) )
        if (asteps != 0):
            fstream.write("  asteps = %i\n" % (self.asteps) )
        if (nfe != 0):
            fstream.write("  nfe = %i\n" % (self.nfe) )
        if (nfi != 0):
            fstream.write("  nfi = %i\n" % (self.nfi) )
        if (lsetups != 0):
            fstream.write("  lsetups = %i\n" % (self.lsetups) )
        if (nfi_lsetup != 0):
            fstream.write("  nfi_lsetup = %i\n" % (self.nfi_lsetup) )
        if (nJe != 0):
            fstream.write("  nJe = %i\n" % (self.nJe) )
        if (nnewt != 0):
            fstream.write("  nnewt = %i\n" % (self.nnewt) )
        if (ncf != 0):
            fstream.write("  ncf = %i\n" % (self.ncf) )
        if (nef != 0):
            fstream.write("  nef = %i\n" % (self.nef) )
        if (maxerr != 0):
            fstream.write("  maxerr = %g\n" % (self.maxerr) )
        if (rmserr != 0):
            fstream.write("  rmserr = %g\n" % (self.rmserr) )
        if (oversolve != 0):
            fstream.write("  oversolve = %g\n" % (self.oversolve) )
        if (runtime != 0):
            fstream.write("  runtime = %g\n" % (self.runtime) )



#### Utility functions ####

def ReadTests(fname):
    """ This routine takes in a string containing a filename """
    """ and reads a list of Test objects from that file.     """
    """ Comments are allowed, as long as those lines begin   """
    """ with the "#" character. """
    TestList = []
    for line in open(fname):
        li=line.strip()
        if not li.startswith("#"):
            TestList.append(Test(line=li))
    return TestList

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
    stats = RunStats()
    tstart = time()
    try:
        subprocess.call(testname, stdout=out, stderr=null)
    except OSError:
        print 'run_test error running executable ', testname
        print '   '
        print 'os.stat() on the requested executable:'
        os.stat(testname)
        print '   '
        print 'directory contents:'
        os.listdir('.')
        print '   '
    stats.runtime = (time() - tstart)
    null.close()
    out.close()
    #  parse output
    f = open('output.txt','r')
    for line in f:
        txt = shlex.split(line)
        if ("stopping" in txt): 
            stats.errfail = 1
            break
        elif ("Internal" in txt):
            stats.nsteps = int(txt[4]);
            stats.asteps = int(txt[7].replace(')',''));
        elif ("system" in txt):
            stats.nfi_lsetup = int(txt[10]);
        elif ("evals:" in txt):
            stats.nfe = int(txt[5].replace(',',''));
            stats.nfi = int(txt[8]);
        elif ("setups" in txt):
            stats.lsetups = int(txt[5]);
        elif ("Jacobian" in txt):
            stats.nJe = int(txt[6]);
        elif (("nonlinear" in txt) and ("iterations" in txt)):
            stats.nnewt = int(txt[6]);
        elif (("convergence" in txt) and ("nonlinear" in txt) and ("Total" in txt)):
            stats.ncf = int(txt[8]);
        elif (("error" in txt) and ("test" in txt)):
            stats.nef = int(txt[7]);
        elif ("rms" in txt):
            stats.maxerr = float(txt[3].replace(',',''));
            stats.rmserr = float(txt[6]);
        elif ("Oversolve" in txt):
            stats.oversolve = float(txt[2]);
    f.close()
    if (keep_output==0):
        os.remove('output.txt')
    return stats
    
    
##########
def write_parameter_file(params):
    """ This routine takes in a set of solver parameterrs and writes """
    """ an input file for the tests to supply these parameters to    """
    """ their solvers.                                               """
    f = open('solve_params.txt', 'w')
    f.write("order = %i\n" % (params.order)) 
    f.write("dense_order = %i\n" % (params.dense_order)) 
    f.write("imex = %i\n" % (params.imex)) 
    f.write("btable = %i\n" % (params.btable)) 
    f.write("adapt_method = %i\n" % (params.adapt_method)) 
    f.write("maxnef = 0\n") 
    f.write("maxncf = 0\n") 
    f.write("mxhnil = 0\n") 
    f.write("mxsteps = 0\n") 
    f.write("cflfac = %g\n" % (params.cflfac)) 
    f.write("safety = %g\n" % (params.safety)) 
    f.write("bias = %g\n" % (params.bias)) 
    f.write("growth = %g\n" % (params.growth)) 
    f.write("hfixed_lb = %g\n" % (params.hfixed_lb)) 
    f.write("hfixed_ub = %g\n" % (params.hfixed_ub)) 
    f.write("pq = %i\n" % (params.pq)) 
    f.write("k1 = %g\n" % (params.k1)) 
    f.write("k2 = %g\n" % (params.k2)) 
    f.write("k3 = %g\n" % (params.k3))
    f.write("etamx1 = %g\n" % (params.etamx1)) 
    f.write("etamxf = %g\n" % (params.etamxf)) 
    f.write("etacf = %g\n" % (params.etacf)) 
    f.write("small_nef = %i\n" % (params.small_nef)) 
    f.write("crdown = %g\n" % (params.crdown)) 
    f.write("rdiv = %g\n" % (params.rdiv)) 
    f.write("dgmax = %g\n" % (params.dgmax)) 
    f.write("predictor = %i\n" % (params.predictor)) 
    f.write("msbp = %i\n" % (params.msbp)) 
    f.write("fixedpt = %i\n" % (params.fixedpt)) 
    f.write("m_aa = %i\n" % (params.m_aa)) 
    f.write("maxcor = %i\n" % (params.maxcor)) 
    f.write("nlscoef = %g\n" % (params.nlscoef)) 
    f.write("h0 = 0.0\n") 
    f.write("hmin = 0.0\n") 
    f.write("hmax = 0.0\n") 
    f.write("rtol = %g\n" % (params.rtol)) 
    f.write("atol = %g\n" % (params.atol)) 
    f.close()

    f = open('fsolve_params.txt', 'w')
    f.write("&inputs\n")
    f.write("  order = %i,\n" % (params.order)) 
    f.write("  dense_order = %i,\n" % (params.dense_order)) 
    f.write("  imex = %i,\n" % (params.imex)) 
    f.write("  btable = %i,\n" % (params.btable))
    f.write("  adapt_method = %i,\n" % (params.adapt_method)) 
    f.write("  cflfac = %g,\n" % (params.cflfac)) 
    f.write("  safety = %g,\n" % (params.safety)) 
    f.write("  bias = %g,\n" % (params.bias)) 
    f.write("  growth = %g,\n" % (params.growth)) 
    f.write("  hfixed_lb = %g,\n" % (params.hfixed_lb)) 
    f.write("  hfixed_ub = %g,\n" % (params.hfixed_ub)) 
    f.write("  pq = %i,\n" % (params.pq)) 
    f.write("  k1 = %g,\n" % (params.k1)) 
    f.write("  k2 = %g,\n" % (params.k2)) 
    f.write("  k3 = %g,\n" % (params.k3))
    f.write("  etamx1 = %g,\n" % (params.etamx1)) 
    f.write("  etamxf = %g,\n" % (params.etamxf)) 
    f.write("  etacf = %g,\n" % (params.etacf)) 
    f.write("  small_nef = %i,\n" % (params.small_nef)) 
    f.write("  crdown = %g,\n" % (params.crdown)) 
    f.write("  rdiv = %g,\n" % (params.rdiv)) 
    f.write("  dgmax = %g,\n" % (params.dgmax)) 
    f.write("  predictor = %i,\n" % (params.predictor)) 
    f.write("  msbp = %i,\n" % (params.msbp)) 
    f.write("  fixedpt = %i,\n" % (params.fixedpt)) 
    f.write("  m_aa = %i,\n" % (params.m_aa)) 
    f.write("  maxcor = %i,\n" % (params.maxcor)) 
    f.write("  nlscoef = %g,\n" % (params.nlscoef)) 
    f.write("&end\n")
    f.close()
    

##### end of script #####
