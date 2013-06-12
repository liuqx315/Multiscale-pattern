#!/usr/bin/env python
# analysis script for ARKODE solver statistics
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import sys
import diag_tools as diags

# load the diagnostics output file (passed on the command line)
fname = sys.argv[len(sys.argv)-1];

# load the time step data 
Tdiags = diags.load_diags(fname);

# generate plots 
diags.plot_h_vs_t(Tdiags,'h_vs_t.png');
diags.plot_h_vs_iter(Tdiags,'h_vs_iter.png');
diags.plot_work_vs_t(Tdiags,'work_vs_t.png');
diags.plot_work_vs_h(Tdiags,'work_vs_h.png');
diags.plot_oversolve_vs_t(Tdiags,'oversolve_vs_t.png');

# print solve statistics to screen
diags.etest_stats(Tdiags, sys.stdout);
diags.solver_stats(Tdiags, sys.stdout);

##### end of script #####
