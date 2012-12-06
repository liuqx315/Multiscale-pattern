# analysis functions for running ARKODE tests with various parameters
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import sys
import numpy as np
import arkode_tools as ark
import time as time

##########
def run_cost(testname,CM):
    """ This routine takes in a string containing an executable name  """
    """ and a cost model array, runs that test sending all output to  """
    """ a temporary file, processes that file to determine the        """
    """ overall run statistics, and assesses a total run 'cost'.      """
    # run test and get statistics
    [nst,ast,cst,nfe,nfi,lset,nfi_lset,nJe,nnewt,ncf,nef,merr,rerr,ov] = ark.run_test(testname)
    # generate total cost
    cost = 1.0 + CM[0]*nst + CM[1]*nfe + CM[2]*(nfi+nfi_lset) + CM[4]*nJe
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
    print 'Total number of tests:\n', ntot
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
    p = ark.SolParams(-1.0, order, -1, imex, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0);
    ark.write_parameter_file(p);
    base = set_baseline(tests, CM);
    p.cost = 1.0;

    # create list of saved choices
    optimal = [];  
    optimal.append(p);
    isaved = 0;

    # loop over all possible parameters combinations
    ntested = 0;
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
                       if (p.cost != p.cost):
                           p.cost = 10000.0;

                       # output progress
                       ntested += 1
                       if ((1.0*ntested/ntot >= 0.1) and (1.0*(ntested-1)/ntot < 0.1)):
                           sys.stdout.write(" 10% complete\n")
                       if ((1.0*ntested/ntot >= 0.2) and (1.0*(ntested-1)/ntot < 0.2)):
                           sys.stdout.write(" 20% complete\n")
                       if ((1.0*ntested/ntot >= 0.3) and (1.0*(ntested-1)/ntot < 0.3)):
                           sys.stdout.write(" 30% complete\n")
                       if ((1.0*ntested/ntot >= 0.4) and (1.0*(ntested-1)/ntot < 0.4)):
                           sys.stdout.write(" 40% complete\n")
                       if ((1.0*ntested/ntot >= 0.5) and (1.0*(ntested-1)/ntot < 0.5)):
                           sys.stdout.write(" 50% complete\n")
                       if ((1.0*ntested/ntot >= 0.6) and (1.0*(ntested-1)/ntot < 0.6)):
                           sys.stdout.write(" 60% complete\n")
                       if ((1.0*ntested/ntot >= 0.7) and (1.0*(ntested-1)/ntot < 0.7)):
                           sys.stdout.write(" 70% complete\n")
                       if ((1.0*ntested/ntot >= 0.8) and (1.0*(ntested-1)/ntot < 0.8)):
                           sys.stdout.write(" 80% complete\n")
                       if ((1.0*ntested/ntot >= 0.9) and (1.0*(ntested-1)/ntot < 0.9)):
                           sys.stdout.write(" 90% complete\n")

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
                           
##########
def parameter_rand_search(order, dense_order, imex, adapt_method, cflfac, 
                          safety, biases, growth, hfixed_lb, hfixed_ub, k1vals, 
                          k2vals, k3vals, etamx1, etamxf, etacf, small_nef, 
                          crdown, rdiv, dgmax, predictor, msbpvals, maxcor, 
                          nlscoef, nsaved, tests, CM, ntries):
    """ This routine performs a stochastic optimization over the """
    """ set of possible solver parameters.  Each input should be """ 
    """ an array of length 2, storing the bounds of allowed      """
    """ values for the given parameter.  The only exceptions to  """
    """ this rule are order (a single integer), imex (a single   """
    """ integer) and CM (a single cost model).                   """

    # create parameter file of all defaults, set baseline cost
    p = ark.SolParams(-1.0, order, -1, imex, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                       0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0);
    ark.write_parameter_file(p);
    base = set_baseline(tests, CM);
    p.cost = 1.0

    # create list of saved choices
    optimal = [];  
    optimal.append(p);
    isaved = 0;

    # bound update frequency, counter
    nupdate = max(nsaved, ntries/10);
    iupdate = 0;

    # set scalar variables for the current bounds on each value
    dord_l = dense_order[0];    dord_r = dense_order[1];
    hmeth_l = adapt_method[0];  hmeth_r = adapt_method[1];
    snef_l = small_nef[0];      snef_r = small_nef[1];
    pred_l = predictor[0];      pred_r = predictor[1];
    msbp_l = msbpvals[0];       msbp_r = msbpvals[1];
    mxcr_l = maxcor[0];         mxcr_r = maxcor[1];
    cfl_l = cflfac[0];          cfl_r = cflfac[1];
    safe_l = safety[0];         safe_r = safety[1];
    bias_l = biases[0];         bias_r = biases[1];
    grow_l = growth[0];         grow_r = growth[1];
    hf_lb_l = hfixed_lb[0];     hf_lb_r = hfixed_lb[1];
    hf_ub_l = hfixed_ub[0];     hf_ub_r = hfixed_ub[1];
    k1_l = k1vals[0];           k1_r = k1vals[1];
    k2_l = k2vals[0];           k2_r = k2vals[1];
    k3_l = k3vals[0];           k3_r = k3vals[1];
    emx1_l = etamx1[0];         emx1_r = etamx1[1];
    emxf_l = etamxf[0];         emxf_r = etamxf[1];
    ecf_l = etacf[0];           ecf_r = etacf[1];
    crd_l = crdown[0];          crd_r = crdown[1];
    rdv_l = rdiv[0];            rdv_r = rdiv[1];
    dgmx_l = dgmax[0];          dgmx_r = dgmax[1];
    nlsc_l = nlscoef[0];        nlsc_r = nlscoef[1];

    # loop over all possible parameters combinations
    ntested = 0;
    sys.stdout.write("Total number of tests to run: %i\n" % (ntries));
    tstart = time.time()
    for ntested in range(ntries):

        # generate real and integer random parameter values in input interval
        dord  = np.random.random_integers(dord_l,  dord_r);
        hmeth = np.random.random_integers(hmeth_l, hmeth_r);
        snef  = np.random.random_integers(snef_l,  snef_r);
        pred  = np.random.random_integers(pred_l,  pred_r);
        msbp  = np.random.random_integers(msbp_l,  msbp_r);
        mxcr  = np.random.random_integers(mxcr_l,  mxcr_r);
        cfl   = cfl_l   + np.random.random_sample()*(cfl_r   - cfl_l);
        safe  = safe_l  + np.random.random_sample()*(safe_r  - safe_l);
        bias  = bias_l  + np.random.random_sample()*(bias_r  - bias_l);
        grow  = grow_l  + np.random.random_sample()*(grow_r  - grow_l);
        hf_lb = hf_lb_l + np.random.random_sample()*(hf_lb_r - hf_lb_l);
        hf_ub = hf_ub_l + np.random.random_sample()*(hf_ub_r - hf_ub_l);
        k1    = k1_l    + np.random.random_sample()*(k1_r    - k1_l);
        k2    = k2_l    + np.random.random_sample()*(k2_r    - k2_l);
        k3    = k3_l    + np.random.random_sample()*(k3_r    - k3_l);
        emx1  = emx1_l  + np.random.random_sample()*(emx1_r  - emx1_l);
        emxf  = emxf_l  + np.random.random_sample()*(emxf_r  - emxf_l);
        ecf   = ecf_l   + np.random.random_sample()*(ecf_r   - ecf_l);
        crd   = crd_l   + np.random.random_sample()*(crd_r   - crd_l);
        rdv   = rdv_l   + np.random.random_sample()*(rdv_r   - rdv_l);
        dgmx  = dgmx_l  + np.random.random_sample()*(dgmx_r  - dgmx_l);
        nlsc  = nlsc_l  + np.random.random_sample()*(nlsc_r  - nlsc_l);

        # create parameter object
        p = ark.SolParams(-1.0, order, dord, imex, hmeth, cfl, safe, bias, grow, 
                           hf_lb, hf_ub, k1, k2, k3, emx1, emxf, ecf, snef, crd, 
                           rdv, dgmx, pred, msbp, mxcr, nlsc);
                         
        # create parameter file 
        ark.write_parameter_file(p);

        # run tests to evaluate overall cost
        p.cost = calc_cost(tests, base, CM)
        if (p.cost != p.cost):
            p.cost = 10000.0;

        # output progress
        if ((1.0*ntested/ntries >= 0.1) and (1.0*(ntested-1)/ntries < 0.1)):
            elapsed = (time.time() - tstart)
            eta = 9.0/1.0*elapsed;
            sys.stdout.write(" 1/10 complete (%g seconds, %g remaining)\n" % (elapsed,eta))
        if ((1.0*ntested/ntries >= 0.2) and (1.0*(ntested-1)/ntries < 0.2)):
            elapsed = (time.time() - tstart)
            eta = 8.0/2.0*elapsed;
            sys.stdout.write(" 2/10 complete (%g seconds, %g remaining)\n" % (elapsed,eta))
        if ((1.0*ntested/ntries >= 0.3) and (1.0*(ntested-1)/ntries < 0.3)):
            elapsed = (time.time() - tstart)
            eta = 7.0/3.0*elapsed;
            sys.stdout.write(" 3/10 complete (%g seconds, %g remaining)\n" % (elapsed,eta))
        if ((1.0*ntested/ntries >= 0.4) and (1.0*(ntested-1)/ntries < 0.4)):
            elapsed = (time.time() - tstart)
            eta = 6.0/4.0*elapsed;
            sys.stdout.write(" 4/10 complete (%g seconds, %g remaining)\n" % (elapsed,eta))
        if ((1.0*ntested/ntries >= 0.5) and (1.0*(ntested-1)/ntries < 0.5)):
            elapsed = (time.time() - tstart)
            eta = 5.0/5.0*elapsed;
            sys.stdout.write(" 5/10 complete (%g seconds, %g remaining)\n" % (elapsed,eta))
        if ((1.0*ntested/ntries >= 0.6) and (1.0*(ntested-1)/ntries < 0.6)):
            elapsed = (time.time() - tstart)
            eta = 4.0/6.0*elapsed;
            sys.stdout.write(" 6/10 complete (%g seconds, %g remaining)\n" % (elapsed,eta))
        if ((1.0*ntested/ntries >= 0.7) and (1.0*(ntested-1)/ntries < 0.7)):
            elapsed = (time.time() - tstart)
            eta = 3.0/7.0*elapsed;
            sys.stdout.write(" 7/10 complete (%g seconds, %g remaining)\n" % (elapsed,eta))
        if ((1.0*ntested/ntries >= 0.8) and (1.0*(ntested-1)/ntries < 0.8)):
            elapsed = (time.time() - tstart)
            eta = 2.0/8.0*elapsed;
            sys.stdout.write(" 8/10 complete (%g seconds, %g remaining)\n" % (elapsed,eta))
        if ((1.0*ntested/ntries >= 0.9) and (1.0*(ntested-1)/ntries < 0.9)):
            elapsed = (time.time() - tstart)
            eta = 1.0/9.0*elapsed;
            sys.stdout.write(" 9/10 complete (%g seconds, %g remaining)\n" % (elapsed,eta))

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

        # increment the interval update counter
        iupdate += 1
        
        # update intervals if necessary
        if (iupdate >= nupdate):
            iupdate = 0;
            oldbounds = ((dord_l, dord_r), (hmeth_l, hmeth_r), (snef_l, snef_r), 
                         (pred_l, pred_r), (msbp_l, msbp_r), (mxcr_l, mxcr_r), 
                         (cfl_l, cfl_r), (safe_l, safe_r), (bias_l, bias_r), 
                         (grow_l, grow_r), (hf_lb_l, hf_lb_r), (hf_ub_l, hf_ub_r), 
                         (k1_l, k1_r), (k2_l, k2_r), (k3_l, k3_r), (emx1_l, emx1_r), 
                         (emxf_l, emxf_r), (ecf_l, ecf_r), (crd_l, crd_r), 
                         (rdv_l, rdv_r), (dgmx_l, dgmx_r), (nlsc_l, nlsc_r));
            newbounds = update_param_bounds(optimal, oldbounds);
            dord_l = newbounds[0][0];  dord_r = newbounds[0][1];
            hmeth_l = newbounds[1][0];  hmeth_r = newbounds[1][1];
            snef_l = newbounds[2][0];   snef_r = newbounds[2][1];
            pred_l = newbounds[3][0];   pred_r = newbounds[3][1];
            msbp_l = newbounds[4][0];   msbp_r = newbounds[4][1];
            mxcr_l = newbounds[5][0];   mxcr_r = newbounds[5][1];
            cfl_l = newbounds[6][0];    cfl_r = newbounds[6][1];
            safe_l = newbounds[7][0];   safe_r = newbounds[7][1];
            bias_l = newbounds[8][0];   bias_r = newbounds[8][1];
            grow_l = newbounds[9][0];   grow_r = newbounds[9][1];
            hf_lb_l = newbounds[10][0];  hf_lb_r = newbounds[10][1];
            hf_ub_l = newbounds[11][0];  hf_ub_r = newbounds[11][1];
            k1_l = newbounds[12][0];     k1_r = newbounds[12][1];
            k2_l = newbounds[13][0];     k2_r = newbounds[13][1];
            k3_l = newbounds[14][0];     k3_r = newbounds[14][1];
            emx1_l = newbounds[15][0];   emx1_r = newbounds[15][1];
            emxf_l = newbounds[16][0];   emxf_r = newbounds[16][1];
            ecf_l = newbounds[17][0];    ecf_r = newbounds[17][1];
            crd_l = newbounds[18][0];    crd_r = newbounds[18][1];
            rdv_l = newbounds[19][0];    rdv_r = newbounds[19][1];
            dgmx_l = newbounds[20][0];   dgmx_r = newbounds[20][1];
            nlsc_l = newbounds[21][0];   nlsc_r = newbounds[21][1];
            
    sys.stdout.write("\n")
    return optimal;
                           
##########
def update_param_bounds(parameters, oldbounds):
    """ This routine searches through a set of parameters and existing  """                
    """ bounds to determine a reduced set of bounds on those parameters """
    """ that encompass everything in the parameter set.                 """

    # set initial lower/upper bounds on each parameter from first allowable set
    dord_l = 100;   dord_r = -100;
    hmeth_l = 100;  hmeth_r = -100;
    cfl_l = 100;    cfl_r = -100;
    safe_l = 100;   safe_r = -100;
    bias_l = 100;   bias_r = -100;
    grow_l = 100;   grow_r = -100;
    hflb_l = 100;   hflb_r = -100;
    hfub_l = 100;   hfub_r = -100;
    k1_l = 100;     k1_r = -100;
    k2_l = 100;     k2_r = -100;
    k3_l = 100;     k3_r = -100;
    em1_l = 100;    em1_r = -100;
    emf_l = 100;    emf_r = -100;
    ecf_l = 100;    ecf_r = -100;
    snf_l = 100;    snf_r = -100;
    crd_l = 100;    crd_r = -100;
    rdv_l = 100;    rdv_r = -100;
    dgm_l = 100;    dgm_r = -100;
    prd_l = 100;    prd_r = -100;
    msb_l = 10000;  msb_r = -100;
    mxc_l = 1000;   mxc_r = -100;
    nls_l = 1000;   nls_r = -100;

    # iterate through parameter list, setting interval based on min/max allowed vals
    for i in range(len(parameters)):
        if ((parameters[i].dense_order > oldbounds[0][0]) and 
            (parameters[i].dense_order < oldbounds[0][1])):
            dord_l = min(dord_l, parameters[i].dense_order)
            dord_r = max(dord_r, parameters[i].dense_order)
        if ((parameters[i].adapt_method > oldbounds[1][0]) and 
            (parameters[i].adapt_method < oldbounds[1][1])):
            hmeth_l = min(hmeth_l, parameters[i].adapt_method)
            hmeth_r = max(hmeth_r, parameters[i].adapt_method)
        if ((parameters[i].cflfac > oldbounds[6][0]) and 
            (parameters[i].cflfac < oldbounds[6][1])):
            cfl_l = min(cfl_l, parameters[i].cflfac)
            cfl_r = max(cfl_r, parameters[i].cflfac)
        if ((parameters[i].safety > oldbounds[7][0]) and 
            (parameters[i].safety < oldbounds[7][1])):
            safe_l = min(safe_l, parameters[i].safety)
            safe_r = max(safe_r, parameters[i].safety)
        if ((parameters[i].bias > oldbounds[8][0]) and 
            (parameters[i].bias < oldbounds[8][1])):
            bias_l = min(bias_l, parameters[i].bias)
            bias_r = max(bias_r, parameters[i].bias)
        if ((parameters[i].growth > oldbounds[9][0]) and 
            (parameters[i].growth < oldbounds[9][1])):
            grow_l = min(grow_l, parameters[i].growth)
            grow_r = max(grow_r, parameters[i].growth)
        if ((parameters[i].hfixed_lb > oldbounds[10][0]) and 
            (parameters[i].hfixed_lb < oldbounds[10][1])):
            hflb_l = min(hflb_l, parameters[i].hfixed_lb)
            hflb_r = max(hflb_r, parameters[i].hfixed_lb)
        if ((parameters[i].hfixed_ub > oldbounds[11][0]) and 
            (parameters[i].hfixed_ub < oldbounds[11][1])):
            hfub_l = min(hfub_l, parameters[i].hfixed_ub)
            hfub_r = max(hfub_r, parameters[i].hfixed_ub)
        if ((parameters[i].k1 > oldbounds[12][0]) and 
            (parameters[i].k1 < oldbounds[12][1])):
            k1_l = min(k1_l, parameters[i].k1)
            k1_r = max(k1_r, parameters[i].k1)
        if ((parameters[i].k2 > oldbounds[13][0]) and 
            (parameters[i].k2 < oldbounds[13][1])):
            k2_l = min(k2_l, parameters[i].k2)
            k2_r = max(k2_r, parameters[i].k2)
        if ((parameters[i].k3 > oldbounds[14][0]) and 
            (parameters[i].k3 < oldbounds[14][1])):
            k3_l = min(k3_l, parameters[i].k3)
            k3_r = max(k3_r, parameters[i].k3)
        if ((parameters[i].etamx1 > oldbounds[15][0]) and 
            (parameters[i].etamx1 < oldbounds[15][1])):
            em1_l = min(em1_l, parameters[i].etamx1)
            em1_r = max(em1_r, parameters[i].etamx1)
        if ((parameters[i].etamxf > oldbounds[16][0]) and 
            (parameters[i].etamxf < oldbounds[16][1])):
            emf_l = min(emf_l, parameters[i].etamxf)
            emf_r = max(emf_r, parameters[i].etamxf)
        if ((parameters[i].etacf > oldbounds[17][0]) and 
            (parameters[i].etacf < oldbounds[17][1])):
            ecf_l = min(ecf_l, parameters[i].etacf)
            ecf_r = max(ecf_r, parameters[i].etacf)
        if ((parameters[i].small_nef > oldbounds[2][0]) and 
            (parameters[i].small_nef < oldbounds[2][1])):
            snf_l = min(snf_l, parameters[i].small_nef)
            snf_r = max(snf_r, parameters[i].small_nef)
        if ((parameters[i].crdown > oldbounds[18][0]) and 
            (parameters[i].crdown < oldbounds[18][1])):
            crd_l = min(crd_l, parameters[i].crdown)
            crd_r = max(crd_r, parameters[i].crdown)
        if ((parameters[i].rdiv > oldbounds[19][0]) and 
            (parameters[i].rdiv < oldbounds[19][1])):
            rdv_l = min(rdv_l, parameters[i].rdiv)
            rdv_r = max(rdv_r, parameters[i].rdiv)
        if ((parameters[i].dgmax > oldbounds[20][0]) and 
            (parameters[i].dgmax < oldbounds[20][1])):
            dgm_l = min(dgm_l, parameters[i].dgmax)
            dgm_r = max(dgm_r, parameters[i].dgmax)
        if ((parameters[i].predictor > oldbounds[3][0]) and 
            (parameters[i].predictor < oldbounds[3][1])):
            prd_l = min(prd_l, parameters[i].predictor)
            prd_r = max(prd_r, parameters[i].predictor)
        if ((parameters[i].msbp > oldbounds[4][0]) and 
            (parameters[i].msbp < oldbounds[4][1])):
            msb_l = min(msb_l, parameters[i].msbp)
            msb_r = max(msb_r, parameters[i].msbp)
        if ((parameters[i].maxcor > oldbounds[5][0]) and 
            (parameters[i].maxcor < oldbounds[5][1])):
            mxc_l = min(mxc_l, parameters[i].maxcor)
            mxc_r = max(mxc_r, parameters[i].maxcor)
        if ((parameters[i].nlscoef > oldbounds[21][0]) and 
            (parameters[i].nlscoef < oldbounds[21][1])):
            nls_l = min(nls_l, parameters[i].nlscoef)
            nls_r = max(nls_r, parameters[i].nlscoef)
        
    # update intervals 
    if (dord_l > dord_r):
        dord_l = oldbounds[0][0];  dord_r = oldbounds[0][1];
    if (hmeth_l > hmeth_r):
        hmeth_l = oldbounds[1][0];  hmeth_r = oldbounds[1][1];
    if (cfl_l > cfl_r):
        cfl_l = oldbounds[6][0];  cfl_r = oldbounds[6][1];
    if (safe_l > safe_r):
        safe_l = oldbounds[7][0]; safe_r = oldbounds[7][1];
    if (bias_l > bias_r):
        bias_l = oldbounds[8][0]; bias_r = oldbounds[8][1];
    if (grow_l > grow_r):
        grow_l = oldbounds[9][0]; grow_r = oldbounds[9][1];
    if (hflb_l > hflb_r):
        hflb_l = oldbounds[10][0]; hflb_r = oldbounds[10][1];
    if (hfub_l > hfub_r):
        hfub_l = oldbounds[11][0]; hfub_r = oldbounds[11][1];
    if (k1_l > k1_r):
        k1_l = oldbounds[12][0]; k1_r = oldbounds[12][1];
    if (k2_l > k2_r):
        k2_l = oldbounds[13][0]; k2_r = oldbounds[13][1];
    if (k3_l > k3_r):
        k3_l = oldbounds[14][0]; k3_r = oldbounds[14][1];
    if (em1_l > em1_r):
        em1_l = oldbounds[15][0]; em1_r = oldbounds[15][1];
    if (emf_l > emf_r):
        emf_l = oldbounds[16][0]; emf_r = oldbounds[16][1];
    if (ecf_l > ecf_r):
        ecf_l = oldbounds[17][0]; ecf_r = oldbounds[17][1];
    if (snf_l > snf_r):
        snf_l = oldbounds[2][0]; snf_r = oldbounds[2][1];
    if (crd_l > crd_r):
        crd_l = oldbounds[18][0]; crd_r = oldbounds[18][1];
    if (rdv_l > rdv_r):
        rdv_l = oldbounds[19][0]; rdv_r = oldbounds[19][1];
    if (dgm_l > dgm_r):
        dgm_l = oldbounds[20][0]; dgm_r = oldbounds[20][1];
    if (prd_l > prd_r):
        prd_l = oldbounds[3][0]; prd_r = oldbounds[3][1];
    if (msb_l > msb_r):
        msb_l = oldbounds[4][0]; msb_r = oldbounds[4][1];
    if (mxc_l > mxc_r):
        mxc_l = oldbounds[5][0]; mxc_r = oldbounds[5][1];
    if (nls_l > nls_r):
        nls_l = oldbounds[21][0]; nls_r = oldbounds[21][1];
    newbounds = ((dord_l, dord_r), (hmeth_l, hmeth_r), (snf_l, snf_r), 
                 (prd_l, prd_r), (msb_l, msb_r), (mxc_l, mxc_r), 
                 (cfl_l, cfl_r), (safe_l, safe_r), (bias_l, bias_r), 
                 (grow_l, grow_r), (hflb_l, hflb_r), (hfub_l, hfub_r), 
                 (k1_l, k1_r), (k2_l, k2_r), (k3_l, k3_r), (em1_l, em1_r), 
                 (emf_l, emf_r), (ecf_l, ecf_r), (crd_l, crd_r), 
                 (rdv_l, rdv_r), (dgm_l, dgm_r), (nls_l, nls_r));
    return newbounds;
                           
##### end of script #####
