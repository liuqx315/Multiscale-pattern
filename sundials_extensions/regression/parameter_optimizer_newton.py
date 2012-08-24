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

# set lists of available parameters for search
# dense_order = (0, 1, 2, 3, 4, 5);
# adapt_method = (0, 1, 2, 3, 4, 5);
# cflfac = 0.0;     # must be a scalar
# safety = (0.9, 0.925, 0.95, 0.975, 1.0);
# bias = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
# growth = (2.0, 5.0, 10.0, 20.0, 50.0);
# hfixed_lb = 0.0;  # must be a scalar
# hfixed_ub = 0.0;  # must be a scalar
# k1 = []; k1.append(0.0);
# k2 = []; k2.append(0.0);
# k3 = []; k3.append(0.0);
# etamx1 = []; etamx1.append(0.0);
# etamxf = (0.1, 0.2, 0.5, 0.75);
# etacf = (0.1, 0.25, 0.5);
# small_nef = (1, 2, 3);
# crdown = (0.1, 0.2, 0.3, 0.5, 0.75, 1.0);
# rdiv = (1.5, 2.0, 2.5, 3.0);
# dgmax = (0.1, 0.2, 0.3, 0.5, 0.75, 1.0);
# predictor = (0, 1, 2, 3);
# msbp = (2, 5, 10, 20, 50, 100);
# maxcor = (2, 3, 4, 5);
# nlscoef = (0.3, 0.1, 0.03, 0.01, 0.003, 0.001);
dense_order  = [];  dense_order.append(3);
adapt_method = [];  adapt_method.append(0);
cflfac = 0.0;       # must be a scalar
safety = [];        safety.append(0.0);
bias = [];          bias.append(0.0);
growth = [];        growth.append(0.0);
hfixed_lb = 0.0;    # must be a scalar
hfixed_ub = 0.0;    # must be a scalar
k1 = [];            k1.append(0.0);
k2 = [];            k2.append(0.0);
k3 = [];            k3.append(0.0);
etamx1 = [];        etamx1.append(0.0);
etamxf = [];        etamxf.append(0.0);
etacf = [];         etacf.append(0.0);
small_nef = [];     small_nef.append(0);
crdown = (0.1, 0.2, 0.3, 0.5, 0.75, 1.0);
rdiv = (1.5, 2.0, 2.5, 3.0);
dgmax = [];         dgmax.append(0.0);
predictor = (0, 1, 2, 3);
msbp = [];          msbp.append(0);
maxcor = (2, 3, 4, 5);
nlscoef = (0.3, 0.1, 0.03, 0.01, 0.003, 0.001);


# set number of parameter sets to store
nsaved = 25;

# run parameter search for cost model 1, order 3
opt_params13 = po.parameter_search(3, dense_order, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM1);
# output saved parameters and costs
print '\nOrder 3, cost model 1, the ',nsaved,' best sets of solver parameters are:'
opt_params13[0].WriteHeader();
for i in range(min(nsaved,len(opt_params13))):
    opt_params13[i].Write();


# run parameter search for cost model 1, order 4
opt_params14 = po.parameter_search(4, dense_order, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM1);
# output saved parameters and costs
print '\nOrder 4, cost model 1, the ',nsaved,' best sets of solver parameters are:'
opt_params14[0].WriteHeader();
for i in range(min(nsaved,len(opt_params14))):
    opt_params14[i].Write();


# run parameter search for cost model 1, order 5
opt_params15 = po.parameter_search(5, dense_order, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM1);
# output saved parameters and costs
print '\nOrder 5, cost model 1, the ',nsaved,' best sets of solver parameters are:'
opt_params15[0].WriteHeader();
for i in range(min(nsaved,len(opt_params15))):
    opt_params15[i].Write();




# run parameter search for cost model 2, order 3
opt_params23 = po.parameter_search(3, dense_order, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM2);
# output saved parameters and costs
print '\nOrder 3, cost model 2, the ',nsaved,' best sets of solver parameters are:'
opt_params23[0].WriteHeader();
for i in range(min(nsaved,len(opt_params23))):
    opt_params23[i].Write();


# run parameter search for cost model 2, order 4
opt_params24 = po.parameter_search(4, dense_order, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM2);
# output saved parameters and costs
print '\nOrder 4, cost model 2, the ',nsaved,' best sets of solver parameters are:'
opt_params24[0].WriteHeader();
for i in range(min(nsaved,len(opt_params24))):
    opt_params24[i].Write();


# run parameter search for cost model 2, order 5
opt_params25 = po.parameter_search(5, dense_order, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM2);
# output saved parameters and costs
print '\nOrder 5, cost model 2, the ',nsaved,' best sets of solver parameters are:'
opt_params25[0].WriteHeader();
for i in range(min(nsaved,len(opt_params25))):
    opt_params25[i].Write();




# run parameter search for cost model 3, order 3
opt_params33 = po.parameter_search(3, dense_order, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM3);
# output saved parameters and costs
print '\nOrder 3, cost model 3, the ',nsaved,' best sets of solver parameters are:'
opt_params33[0].WriteHeader();
for i in range(min(nsaved,len(opt_params33))):
    opt_params33[i].Write();


# run parameter search for cost model 3, order 4
opt_params34 = po.parameter_search(4, dense_order, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM3);
# output saved parameters and costs
print '\nOrder 4, cost model 3, the ',nsaved,' best sets of solver parameters are:'
opt_params34[0].WriteHeader();
for i in range(min(nsaved,len(opt_params34))):
    opt_params34[i].Write();


# run parameter search for cost model 1, order 5
opt_params35 = po.parameter_search(5, dense_order, adapt_method, cflfac, safety, bias, 
                                growth, hfixed_lb, hfixed_ub, k1, k2, k3, etamx1, 
                                etamxf, etacf, small_nef, crdown, rdiv, dgmax, 
                                predictor, msbp, maxcor, nlscoef, nsaved, tests, CM3);
# output saved parameters and costs
print '\nOrder 5, cost model 3, the ',nsaved,' best sets of solver parameters are:'
opt_params35[0].WriteHeader();
for i in range(min(nsaved,len(opt_params35))):
    opt_params35[i].Write();



##### end of script #####
