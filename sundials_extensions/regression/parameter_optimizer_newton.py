#!/usr/bin/env python
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

# set tolerances
rtol = 1.e-6;
atol = 1.e-10;


# order 3

# set intervals of available parameters to search
maxcor = (3, 3);
crdown = (0.164942, 0.165608);
rdiv = (2.27675, 2.37037);
nlscoef = (0.21559, 0.219672);
dense_order = (-1, -1);
adapt_method = (0, 0);
cflfac = (0.0, 0.0);
safety = (0.0, 0.0);
bias = (0.0, 0.0);
growth = (0.0, 0.0);
hfixed_lb = (0.0, 0.0);
hfixed_ub = (0.0, 0.0);
k1 = (0.0, 0.0);
k2 = (0.0, 0.0);
k3 = (0.0, 0.0);
etamx1 = (0.0, 0.0);
etamxf = (0.0, 0.0);
etacf = (0.0, 0.0);
small_nef = (0, 0);
dgmax = (0.0, 0.0);
predictor = (1, 1);
msbp = (0, 0);


# run parameter search for cost model 1, order 3
opt_params13 = po.parameter_rand_search(3, dense_order, imex, adapt_method, cflfac, 
                                        safety, bias, growth, hfixed_lb, hfixed_ub, 
                                        k1, k2, k3, etamx1, etamxf, etacf, small_nef, 
                                        crdown, rdiv, dgmax, predictor, msbp, maxcor, 
                                        nlscoef, nsaved, tests, CM1, rtol, atol, ntries);
# output saved parameters and costs
print '\nOrder 3, cost model 1, the ',nsaved,' best sets of solver parameters are:'
opt_params13[0].WriteHeader();
for i in range(min(nsaved,len(opt_params13))):
    opt_params13[i].Write();


# run parameter search for cost model 2, order 3
opt_params23 = po.parameter_rand_search(3, dense_order, imex, adapt_method, cflfac, 
                                        safety, bias, growth, hfixed_lb, hfixed_ub, 
                                        k1, k2, k3, etamx1, etamxf, etacf, small_nef, 
                                        crdown, rdiv, dgmax, predictor, msbp, maxcor, 
                                        nlscoef, nsaved, tests, CM2, rtol, atol, ntries);
# output saved parameters and costs
print '\nOrder 3, cost model 2, the ',nsaved,' best sets of solver parameters are:'
opt_params23[0].WriteHeader();
for i in range(min(nsaved,len(opt_params23))):
    opt_params23[i].Write();


# run parameter search for cost model 3, order 3
opt_params33 = po.parameter_rand_search(3, dense_order, imex, adapt_method, cflfac, 
                                        safety, bias, growth, hfixed_lb, hfixed_ub, 
                                        k1, k2, k3, etamx1, etamxf, etacf, small_nef, 
                                        crdown, rdiv, dgmax, predictor, msbp, maxcor, 
                                        nlscoef, nsaved, tests, CM3, rtol, atol, ntries);
# output saved parameters and costs
print '\nOrder 3, cost model 3, the ',nsaved,' best sets of solver parameters are:'
opt_params33[0].WriteHeader();
for i in range(min(nsaved,len(opt_params33))):
    opt_params33[i].Write();




# order 4

# set intervals of available parameters to search
maxcor = (3, 4);
crdown = (0.255285, 0.25564);
rdiv = (2.25, 2.32726);
nlscoef = (0.24319, 0.243755);
dense_order = (-1, -1);
adapt_method = (0, 0);
cflfac = (0.0, 0.0);
safety = (0.0, 0.0);
bias = (0.0, 0.0);
growth = (0.0, 0.0);
hfixed_lb = (0.0, 0.0);
hfixed_ub = (0.0, 0.0);
k1 = (0.0, 0.0);
k2 = (0.0, 0.0);
k3 = (0.0, 0.0);
etamx1 = (0.0, 0.0);
etamxf = (0.0, 0.0);
etacf = (0.0, 0.0);
small_nef = (0, 0);
dgmax = (0.0, 0.0);
predictor = (1, 1);
msbp = (0, 0);


# run parameter search for cost model 1, order 4
opt_params14 = po.parameter_rand_search(4, dense_order, imex, adapt_method, cflfac, 
                                        safety, bias, growth, hfixed_lb, hfixed_ub, 
                                        k1, k2, k3, etamx1, etamxf, etacf, small_nef, 
                                        crdown, rdiv, dgmax, predictor, msbp, maxcor, 
                                        nlscoef, nsaved, tests, CM1, rtol, atol, ntries);
# output saved parameters and costs
print '\nOrder 4, cost model 1, the ',nsaved,' best sets of solver parameters are:'
opt_params14[0].WriteHeader();
for i in range(min(nsaved,len(opt_params14))):
    opt_params14[i].Write();


# run parameter search for cost model 2, order 4
opt_params24 = po.parameter_rand_search(4, dense_order, imex, adapt_method, cflfac, 
                                        safety, bias, growth, hfixed_lb, hfixed_ub, 
                                        k1, k2, k3, etamx1, etamxf, etacf, small_nef, 
                                        crdown, rdiv, dgmax, predictor, msbp, maxcor, 
                                        nlscoef, nsaved, tests, CM2, rtol, atol, ntries);
# output saved parameters and costs
print '\nOrder 4, cost model 2, the ',nsaved,' best sets of solver parameters are:'
opt_params24[0].WriteHeader();
for i in range(min(nsaved,len(opt_params24))):
    opt_params24[i].Write();


# run parameter search for cost model 3, order 4
opt_params34 = po.parameter_rand_search(4, dense_order, imex, adapt_method, cflfac, 
                                        safety, bias, growth, hfixed_lb, hfixed_ub, 
                                        k1, k2, k3, etamx1, etamxf, etacf, small_nef, 
                                        crdown, rdiv, dgmax, predictor, msbp, maxcor, 
                                        nlscoef, nsaved, tests, CM3, rtol, atol, ntries);
# output saved parameters and costs
print '\nOrder 4, cost model 3, the ',nsaved,' best sets of solver parameters are:'
opt_params34[0].WriteHeader();
for i in range(min(nsaved,len(opt_params34))):
    opt_params34[i].Write();




# order 5

# set intervals of available parameters to search
maxcor = (3, 3);
crdown = (0.397253, 0.397406);
rdiv = (2.19831, 2.35371);
nlscoef = (0.254376, 0.254698);
dense_order = (-1, -1);
adapt_method = (0, 0);
cflfac = (0.0, 0.0);
safety = (0.0, 0.0);
bias = (0.0, 0.0);
growth = (0.0, 0.0);
hfixed_lb = (0.0, 0.0);
hfixed_ub = (0.0, 0.0);
k1 = (0.0, 0.0);
k2 = (0.0, 0.0);
k3 = (0.0, 0.0);
etamx1 = (0.0, 0.0);
etamxf = (0.0, 0.0);
etacf = (0.0, 0.0);
small_nef = (0, 0);
dgmax = (0.0, 0.0);
predictor = (1, 1);
msbp = (0, 0);


# run parameter search for cost model 1, order 5
opt_params15 = po.parameter_rand_search(5, dense_order, imex, adapt_method, cflfac, 
                                        safety, bias, growth, hfixed_lb, hfixed_ub, 
                                        k1, k2, k3, etamx1, etamxf, etacf, small_nef, 
                                        crdown, rdiv, dgmax, predictor, msbp, maxcor, 
                                        nlscoef, nsaved, tests, CM1, rtol, atol, ntries);
# output saved parameters and costs
print '\nOrder 5, cost model 1, the ',nsaved,' best sets of solver parameters are:'
opt_params15[0].WriteHeader();
for i in range(min(nsaved,len(opt_params15))):
    opt_params15[i].Write();


# run parameter search for cost model 2, order 5
opt_params25 = po.parameter_rand_search(5, dense_order, imex, adapt_method, cflfac, 
                                        safety, bias, growth, hfixed_lb, hfixed_ub, 
                                        k1, k2, k3, etamx1, etamxf, etacf, small_nef, 
                                        crdown, rdiv, dgmax, predictor, msbp, maxcor, 
                                        nlscoef, nsaved, tests, CM2, rtol, atol, ntries);
# output saved parameters and costs
print '\nOrder 5, cost model 2, the ',nsaved,' best sets of solver parameters are:'
opt_params25[0].WriteHeader();
for i in range(min(nsaved,len(opt_params25))):
    opt_params25[i].Write();


# run parameter search for cost model 1, order 5
opt_params35 = po.parameter_rand_search(5, dense_order, imex, adapt_method, cflfac, 
                                        safety, bias, growth, hfixed_lb, hfixed_ub, 
                                        k1, k2, k3, etamx1, etamxf, etacf, small_nef,
                                        crdown, rdiv, dgmax, predictor, msbp, maxcor,
                                        nlscoef, nsaved, tests, CM3, rtol, atol, ntries);
# output saved parameters and costs
print '\nOrder 5, cost model 3, the ',nsaved,' best sets of solver parameters are:'
opt_params35[0].WriteHeader();
for i in range(min(nsaved,len(opt_params35))):
    opt_params35[i].Write();



##### end of script #####
