#!/usr/bin/env python
# matplotlib-based plotting script for RMHD3d tests
# Daniel R. Reynolds, reynolds@smu.edu

# imports
import sys
import arkode_tools as ark

# load the diagnostics output file (passed on the command line)
fname = sys.argv[len(sys.argv)-1];

# load the time step data 
Tdiags = ark.load_diags(fname);

# generate plots 
ark.plot_h_vs_t(Tdiags,'h_vs_t.png');
ark.plot_h_vs_iter(Tdiags,'h_vs_iter.png');
ark.plot_work_vs_t(Tdiags,'work_vs_t.png');
ark.plot_work_vs_h(Tdiags,'work_vs_h.png');
ark.plot_oversolve_vs_t(Tdiags,'oversolve_vs_t.png');

# print solve statistics to screen
ark.etest_stats(Tdiags, sys.stdout);
ark.solver_stats(Tdiags, sys.stdout);

##### end of script #####
