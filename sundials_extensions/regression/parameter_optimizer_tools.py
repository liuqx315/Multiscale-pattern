# analysis script for ARKODE solver statistics
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import sys
import arkode_tools as ark
import numpy as np

#### Data Structures ####

class SolParams:
    """ Holds a set of solver parameters and its associated cost. """
    def __init__(self, cost, order, dense_order, adapt_method, cflfac, safety, bias, 
                 growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, etamxf, etacf, 
                 small_nef, crdown, rdiv, dgmax, predictor, msbp, maxcor, nlscoef):
        self.cost = cost;
        self.order = order;
        self.dense_order = dense_order;
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
        print '   cost   ord dord adpt cfl   safe  bias  grow  h0l  h0b    k1     k2     k3     emx1     emxf   ecf   smf  crdn  rdiv   dgmx  prd  msbp  mxcr  nlsc'
    def Write(self):
        sys.stdout.write(" %f %2i %3i  %3i   %3.1f  %5.3f  %4.2f  %4.1f  %3.1f  %3.1f  %5.3f  %5.3f  %5.3f  %6f  %5.2f  %5.2f %3i   %4.2f  %5.3f  %4.2f %3i  %3i   %3i   %.1e\n" % 
                         (self.cost, self.order, self.dense_order, self.adapt_method, self.cflfac, 
                          self.safety, self.bias, self.growth, self.hfixed_lb, self.hfixed_ub, 
                          self.k1, self.k2, self.k3, self.etamx1, self.etamxf, self.etacf, 
                          self.small_nef, self.crdown, self.rdiv, self.dgmax, self.predictor,
                          self.msbp, self.maxcor, self.nlscoef))
        #print ' ',self.cost,' ',self.order,' ',self.dense_order,' ',self.adapt_method,' ',self.cflfac,' ',self.safety,' ',self.bias,' ',self.growth,' ',self.hfixed_lb,' ',self.hfixed_ub,' ',self.k1,' ',self.k2,' ',self.k3,' ',self.etamx1,' ',self.etamxf,' ',self.etacf,' ',self.small_nef,' ',self.crdown,' ',self.rdiv,' ',self.dgmax,' ',self.predictor,' ',self.msbp,' ',self.maxcor,' ',self.nlscoef


#### Utility functions ####

def run_test(testname):
    """ This routine takes in a string containing an executable name, """
    """ runs that test sending all output to a temporary file, and    """
    """ processes that file to determine the overall run statistics,  """
    """ which are returned to the calling routine.                    """
    import shlex
    import os
    import subprocess
    #  run routine
    out = open('output.txt','w')
    null = open(os.devnull, 'w')
    subprocess.call(testname, stdout=out, stderr=null)
    null.close()
    out.close()
    #  initialize return values
    nsteps=0
    asteps=0
    csteps=0
    nfe=0
    nfi=0 
    lsetups=0
    nfi_lsetup=0
    nJe=0
    nnewt=0
    ncf=0
    nef=0
    maxerr=0.0
    rmserr=0.0
    oversolve=0.0
    #  parse output
    f = open('output.txt','r')
    for line in f:
        txt = shlex.split(line)
        if ("internal" in txt):
            nsteps = int(txt[5]);
            asteps = int(txt[8].replace(',',''));
            csteps = int(txt[11].replace(')',''));
        elif ("system" in txt):
            nfi_lsetup = int(txt[10]);
        elif ("evals:" in txt):
            nfe = int(txt[5].replace(',',''));
            nfi = int(txt[8]);
        elif ("setups" in txt):
            lsetups = int(txt[5]);
        elif ("Jacobian" in txt):
            nJe = int(txt[6]);
        elif ("Newton" in txt):
            nnewt = int(txt[6]);
        elif ("convergence" in txt):
            ncf = int(txt[8]);
        elif (("error" in txt) and ("test" in txt)):
            nef = int(txt[7]);
        elif ("rms" in txt):
            maxerr = float(txt[3].replace(',',''));
            rmserr = float(txt[6]);
        elif ("Oversolve" in txt):
            oversolve = float(txt[2]);
    f.close()
    os.remove('output.txt')
    return [nsteps, asteps, csteps, nfe, nfi, lsetups, nfi_lsetup, nJe, nnewt, ncf, nef, maxerr, rmserr, oversolve]
    
##########
def write_parameter_file(params):
    """ This routine takes in a set of solver parameterrs and writes """
    """ an input file for the tests to supply these parameters to    """
    """ their solvers.                                               """
    f = open('solve_params.txt', 'w')
    f.write("order = %i\n" % (params.order)) 
    f.write("dense_order = %i\n" % (params.dense_order)) 
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
    

##########
def run_cost(testname,CM):
    """ This routine takes in a string containing an executable name  """
    """ and a cost model array, runs that test sending all output to  """
    """ a temporary file, processes that file to determine the        """
    """ overall run statistics, and assesses a total run 'cost'.      """
    # run test and get statistics
    [nst,ast,cst,nfe,nfi,lset,nfi_lset,nJe,nnewt,ncf,nef,merr,rerr,ov] = run_test(testname)
    # generate total cost
    cost = CM[0]*nst + CM[1]*nfe + CM[2]*(nfi+nfi_lset) + CM[4]*nJe
    cost += CM[3]*lset + CM[5]*nnewt
    if (ov < 1.0):   # undersolve
        cost += CM[7]/ov
    else:            # oversolve
        cost += CM[6]*ov
    return cost

##########
def set_baseline(tests,CM):
    """ This routine takes in an array of strings containing executable """ 
    """ names and a cost model array.  It runs all of the tests and     """
    """ computes their costs according to the supplied cost model, and  """
    """ returns an array containing the baseline costs of each test.    """
    baseline = [];
    for i in range(len(tests)):
        baseline.append(run_cost(tests[i],CM));
    return baseline

##########
def calc_cost(tests,baseline,CM):
    """ This routine takes in an array of strings containing executable """ 
    """ names, an array of baseline costs for each test, and a cost     """
    """ model array.  It runs all of the tests and computes their costs """
    """ according to the supplied cost model, returning a single value  """
    """ corresponding to the effective cost of the runs.                """
    totalcost = 0.0
    ntests = len(tests);
    for i in range(len(tests)):
        cost = run_cost(tests[i],CM);
        totalcost += cost/baseline[i]/ntests;
    return totalcost

##########
def sort_params(param_list):
    """ This routine takes in an array of solver parameters, and sorts  """ 
    """ them by cost.                                                   """
    n = len(param_list);
    for i in range(n):
        # find remaining member of param_list with smallest cost
        mincost = param_list[i].cost;
        minloc  = i;
        for j in range(i,n):
            if (param_list[j].cost < mincost):
                mincost = param_list[j].cost;
                minloc = j;
        # swap lowest with current
        p = param_list[i];
        param_list[i] = param_list[minloc]
        param_list[minloc] = p
    return param_list

##########
def parameter_search(order, dense_order, adapt_method, cflfac, safety, bias, 
                     growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                     etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                     predictor, msbp, maxcor, nlscoef, nsaved, tests, CM):
    """ This routine iterates over all possible combinations of  """
    """ solver parameters, storing 'nsaved' of them to return to """
    """ the calling routine.                                     """
    # calculate total number of parameter options
    ntot = (len(dense_order) * len(adapt_method) * len(safety) * len(bias) * 
            len(growth) * len(k1) * len(k2) * len(k3) * len(etamx1) * 
            len(etamxf) * len(etacf) * len(small_nef) * len(crdown) * 
            len(rdiv) * len(dgmax) * len(predictor) * len(msbp) * 
            len(maxcor) * len(nlscoef));
    print 'Total number of tests for each cost model: ', ntot

    # create list of saved choices
    optimal = [];
    isaved = 0;
    cutoff = 100.0;

    # create parameter file of all defaults, set baseline cost
    p = SolParams(-1.0, 0, -1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                   0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0);
    write_parameter_file(p);
    base = set_baseline(tests, CM);

    # loop over all possible parameters combinations
    ntested = 0;
    sys.stdout.write("% complete:")
    for i1 in range(len(dense_order)):
     for i2 in range(len(adapt_method)):
      for i3 in range(len(safety)):
       for i4 in range(len(bias)):
        for i5 in range(len(growth)):
         for i6 in range(len(k1)):
          for i7 in range(len(k2)):
           for i8 in range(len(k3)):
            for i9 in range(len(etamx1)):
             for i10 in range(len(etamxf)):
              for i11 in range(len(etacf)):
               for i12 in range(len(small_nef)):
                for i13 in range(len(crdown)):
                 for i14 in range(len(rdiv)):
                  for i15 in range(len(dgmax)):
                   for i16 in range(len(predictor)):
                    for i17 in range(len(msbp)):
                     for i18 in range(len(maxcor)):
                      for i19 in range(len(nlscoef)):

                       # create parameter object
                       p = SolParams(-1.0, order, dense_order[i1], 
                                     adapt_method[i2], cflfac,
                                     safety[i3], bias[i4], 
                                     growth[i5], hfixed_lb, 
                                     hfixed_ub, k1[i6], k2[i7],
                                     k3[i8], etamx1[i9], 
                                     etamxf[i10], etacf[i11], 
                                     small_nef[i12], crdown[i13], 
                                     rdiv[i14], dgmax[i15], 
                                     predictor[i16], msbp[i17], 
                                     maxcor[i18], nlscoef[i19]);
                         
                       # create parameter file 
                       write_parameter_file(p);

                       # run tests to evaluate overall cost
                       p.cost = calc_cost(tests, base, CM)

                       # output progress
                       ntested += 1
                       if ((1.0*ntested/ntot >= 0.1) and (1.0*(ntested-1)/ntot < 0.1)):
                           sys.stdout.write(" 10")
                       if ((1.0*ntested/ntot >= 0.2) and (1.0*(ntested-1)/ntot < 0.2)):
                           sys.stdout.write(" 20")
                       if ((1.0*ntested/ntot >= 0.3) and (1.0*(ntested-1)/ntot < 0.3)):
                           sys.stdout.write(" 30")
                       if ((1.0*ntested/ntot >= 0.4) and (1.0*(ntested-1)/ntot < 0.4)):
                           sys.stdout.write(" 40")
                       if ((1.0*ntested/ntot >= 0.5) and (1.0*(ntested-1)/ntot < 0.5)):
                           sys.stdout.write(" 50")
                       if ((1.0*ntested/ntot >= 0.6) and (1.0*(ntested-1)/ntot < 0.6)):
                           sys.stdout.write(" 60")
                       if ((1.0*ntested/ntot >= 0.7) and (1.0*(ntested-1)/ntot < 0.7)):
                           sys.stdout.write(" 70")
                       if ((1.0*ntested/ntot >= 0.8) and (1.0*(ntested-1)/ntot < 0.8)):
                           sys.stdout.write(" 80")
                       if ((1.0*ntested/ntot >= 0.9) and (1.0*(ntested-1)/ntot < 0.9)):
                           sys.stdout.write(" 90")

                       # if structure is not full, just add entry and sort
                       if (isaved < nsaved):
                           optimal.append(p);
                           optimal = sort_params(optimal);
                           isaved = isaved + 1;

                       # otherwise, compare cost to see if it qualifies
                       elif (p.cost < optimal[-1].cost):
                           # store and sort
                           optimal[-1] = p;
                           optimal = sort_params(optimal);

    sys.stdout.write("\n")
    return optimal;
                           
##### end of script #####
