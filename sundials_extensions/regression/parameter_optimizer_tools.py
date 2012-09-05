# analysis functions for running ARKODE tests with various parameters
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import sys
import numpy as np
import arkode_tools as ark

##########
def run_cost(testname,CM):
    """ This routine takes in a string containing an executable name  """
    """ and a cost model array, runs that test sending all output to  """
    """ a temporary file, processes that file to determine the        """
    """ overall run statistics, and assesses a total run 'cost'.      """
    # run test and get statistics
    [nst,ast,cst,nfe,nfi,lset,nfi_lset,nJe,nnewt,ncf,nef,merr,rerr,ov] = ark.run_test(testname)
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
def parameter_search(order, dense_order, imex, adapt_method, cflfac, 
                     safety, bias, growth, hfixed_lb, hfixed_ub, k1, 
                     k2, k3, etamx1, etamxf, etacf, small_nef, 
                     crdown, rdiv, dgmax, predictor, msbp, maxcor, 
                     nlscoef, nsaved, tests, CM):
    """ This routine iterates over all possible combinations of  """
    """ solver parameters, storing 'nsaved' of them to return to """
    """ the calling routine.                                     """
    # calculate total number of parameter options
    ntot = (len(dense_order) * len(adapt_method) * len(safety) * len(bias) * 
            len(growth) * len(k1) * len(k2) * len(k3) * len(etamx1) * 
            len(etamxf) * len(etacf) * len(small_nef) * len(crdown) * 
            len(rdiv) * len(dgmax) * len(predictor) * len(msbp) * 
            len(maxcor) * len(nlscoef));
    print 'Total number of tests: ', ntot
    # print '  len(dense_order): ',len(dense_order)
    # print '  len(adapt_method): ',len(adapt_method)
    # print '  len(safety): ', len(safety)
    # print '  len(bias): ',len(bias)
    # print '  len(growth): ',len(growth)
    # print '  len(k1): ',len(k1)
    # print '  len(k2): ',len(k2)
    # print '  len(k3): ',len(k3)
    # print '  len(etamx1): ',len(etamx1)
    # print '  len(etamxf): ',len(etamxf)
    # print '  len(etacf): ',len(etacf)
    # print '  len(small_nef): ',len(small_nef)
    # print '  len(crdown): ',len(crdown)
    # print '  len(rdiv): ',len(rdiv)
    # print '  len(dgmax): ',len(dgmax)
    # print '  len(predictor): ',len(predictor)
    # print '  len(msbp): ',len(msbp)
    # print '  len(maxcor): ',len(maxcor)
    # print '  len(nlscoef): ',len(nlscoef)

    # create parameter file of all defaults, set baseline cost
    p = ark.SolParams(-1.0, 0, -1, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0);
    ark.write_parameter_file(p);
    base = set_baseline(tests, CM);

    # create list of saved choices
    optimal = [];  
    optimal.append(p);
    isaved = 0;
    cutoff = 100.0;

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
                       p = ark.SolParams(-1.0, order, dense_order[i1], 
                                          imex, adapt_method[i2], cflfac,
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
                       ark.write_parameter_file(p);

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
