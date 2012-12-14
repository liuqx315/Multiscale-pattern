# analysis script for ARKODE solver statistics
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import sys
import numpy as np
import parameter_optimizer_tools as po


# set up a list of executable names to use in tests
tests = ('ark_analytic.exe', 'ark_analytic_nonlin.exe', 'ark_analytic_sys.exe', 'ark_brusselator.exe', 'ark_brusselator1D.exe');

# set up a few cost models: order of entries are the weight factors per:
#     nstep,  nfe, nfi,  lset, nJe, nnewt, oversolve, undersolve
CM1 = (0.0,   1.0, 1.0,   1.0, 2.0,  10.0,    0.0,    100000.0);  # matrix-free
CM2 = (0.0,   1.0, 1.0, 100.0, 5.0,  10.0,    0.0,    100000.0);  # dense Jacobian
CM3 = (0.0, 100.0, 1.0,   1.0, 1.0,  10.0,    0.0,    100000.0);  # imex, costly fe, matrix-free

# set type of integrator to optimize, and number of trials to take, number of params to store
imex = 0;           # scalar: 0=>implicit, 1=>explicit, 2=>imex
ntries = 10000;
nsaved = 20;


# order 3 

# set intervals of available parameters to search
safety = (0.951796, 0.956735);
bias = (1.28973, 1.35246);
growth = (21.6526, 24.3975);
k1 = (0.43088, 0.433799);
k2 = (0.32923, 0.331317);
etamxf = (0.363892, 0.38368);
dense_order  = (-1, -1);
adapt_method = (3, 3);
cflfac = (0.0, 0.0);
hfixed_lb = (0.0, 0.0);
hfixed_ub = (0.0, 0.0);
k3 = (0.0, 0.0);
etamx1 = (0.0, 0.0);
etacf = (0.0, 0.0);
small_nef = (0, 0);
crdown = (0.0, 0.0);
rdiv = (0.0, 0.0);
dgmax = (0.0, 0.0);
predictor = (3, 3);
msbp = (0, 0);
maxcor = (0, 0);
nlscoef = (0.3, 0.3);


# run parameter search for cost model 1, order 3
opt_params13 = po.parameter_rand_search(3, dense_order, imex, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM1, ntries);
# output saved parameters and costs
print '\nOrder 3, cost model 1, the ',nsaved,' best sets of solver parameters are:'
opt_params13[0].WriteHeader();
for i in range(min(nsaved,len(opt_params13))):
    opt_params13[i].Write();


# run parameter search for cost model 2, order 3
opt_params23 = po.parameter_rand_search(3, dense_order, imex, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM2, ntries);
# output saved parameters and costs
print '\nOrder 3, cost model 2, the ',nsaved,' best sets of solver parameters are:'
opt_params23[0].WriteHeader();
for i in range(min(nsaved,len(opt_params23))):
    opt_params23[i].Write();


# run parameter search for cost model 3, order 3
opt_params33 = po.parameter_rand_search(3, dense_order, imex, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM3, ntries);
# output saved parameters and costs
print '\nOrder 3, cost model 3, the ',nsaved,' best sets of solver parameters are:'
opt_params33[0].WriteHeader();
for i in range(min(nsaved,len(opt_params33))):
    opt_params33[i].Write();





# order 4

# set intervals of available parameters to search
safety = (0.989858, 0.991592);
bias = (1.14474, 1.16081);
growth = (25.3979, 26.8551);
k1 = (0.464956, 0.467155);
k2 = (0.316591, 0.317827);
etamxf = (0.261513, 0.291913);
dense_order  = (-1, -1);
adapt_method = (3, 3);
cflfac = (0.0, 0.0);
hfixed_lb = (0.0, 0.0);
hfixed_ub = (0.0, 0.0);
k3 = (0.0, 0.0);
etamx1 = (0.0, 0.0);
etacf = (0.0, 0.0);
small_nef = (0, 0);
crdown = (0.0, 0.0);
rdiv = (0.0, 0.0);
dgmax = (0.0, 0.0);
predictor = (3, 3);
msbp = (0, 0);
maxcor = (0, 0);
nlscoef = (0.3, 0.3);


# run parameter search for cost model 1, order 4
opt_params14 = po.parameter_rand_search(4, dense_order, imex, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM1, ntries);
# output saved parameters and costs
print '\nOrder 4, cost model 1, the ',nsaved,' best sets of solver parameters are:'
opt_params14[0].WriteHeader();
for i in range(min(nsaved,len(opt_params14))):
    opt_params14[i].Write();


# run parameter search for cost model 2, order 4
opt_params24 = po.parameter_rand_search(4, dense_order, imex, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM2, ntries);
# output saved parameters and costs
print '\nOrder 4, cost model 2, the ',nsaved,' best sets of solver parameters are:'
opt_params24[0].WriteHeader();
for i in range(min(nsaved,len(opt_params24))):
    opt_params24[i].Write();


# run parameter search for cost model 3, order 4
opt_params34 = po.parameter_rand_search(4, dense_order, imex, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM3, ntries);
# output saved parameters and costs
print '\nOrder 4, cost model 3, the ',nsaved,' best sets of solver parameters are:'
opt_params34[0].WriteHeader();
for i in range(min(nsaved,len(opt_params34))):
    opt_params34[i].Write();





# order 5

# set intervals of available parameters to search
safety = (0.960199, 0.963675);
bias = (2.85007, 3.39367);
growth = (26.2181, 29.007);
k1 = (0.377189, 0.399051);
k2 = (0.3495, 0.361343);
etamxf = (0.350226, 0.383884);
dense_order  = (-1, -1);
adapt_method = (3, 3);
cflfac = (0.0, 0.0);
hfixed_lb = (0.0, 0.0);
hfixed_ub = (0.0, 0.0);
k3 = (0.0, 0.0);
etamx1 = (0.0, 0.0);
etacf = (0.0, 0.0);
small_nef = (0, 0);
crdown = (0.0, 0.0);
rdiv = (0.0, 0.0);
dgmax = (0.0, 0.0);
predictor = (3, 3);
msbp = (0, 0);
maxcor = (0, 0);
nlscoef = (0.3, 0.3);


# run parameter search for cost model 1, order 5
opt_params15 = po.parameter_rand_search(5, dense_order, imex, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM1, ntries);
# output saved parameters and costs
print '\nOrder 5, cost model 1, the ',nsaved,' best sets of solver parameters are:'
opt_params15[0].WriteHeader();
for i in range(min(nsaved,len(opt_params15))):
    opt_params15[i].Write();


# run parameter search for cost model 2, order 5
opt_params25 = po.parameter_rand_search(5, dense_order, imex, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM2, ntries);
# output saved parameters and costs
print '\nOrder 5, cost model 2, the ',nsaved,' best sets of solver parameters are:'
opt_params25[0].WriteHeader();
for i in range(min(nsaved,len(opt_params25))):
    opt_params25[i].Write();


# run parameter search for cost model 1, order 5
opt_params35 = po.parameter_rand_search(5, dense_order, imex, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM3, ntries);
# output saved parameters and costs
print '\nOrder 5, cost model 3, the ',nsaved,' best sets of solver parameters are:'
opt_params35[0].WriteHeader();
for i in range(min(nsaved,len(opt_params35))):
    opt_params35[i].Write();



##### end of script #####
