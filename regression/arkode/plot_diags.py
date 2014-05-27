#!/usr/bin/env python
#------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
#------------------------------------------------------------
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
#------------------------------------------------------------
# analysis script for ARKODE solver statistics

# imports
import sys
import diag_tools as diags

# load the diagnostics output file (passed on the command line)
fname = sys.argv[len(sys.argv)-1];

# load the time step data 
try:
    Tdiags = diags.load_diags(fname);
except ValueError:
    # if no command-line argument provided, return with an error
    sys.exit("\nCalling syntax error: filename containing diagnostics must be provided, e.g.\n    $ plot_ark_diags.py fname.txt\n")

# generate plots 
diags.plot_h_vs_t(Tdiags,'h_vs_t.png');
diags.plot_h_vs_iter(Tdiags,'h_vs_iter.png');
diags.plot_nonlin_vs_t(Tdiags,'nonlin_vs_t.png');
diags.plot_nonlin_vs_h(Tdiags,'nonlin_vs_h.png');
diags.plot_krylov_vs_t(Tdiags,'krylov_vs_t.png');
diags.plot_krylov_vs_h(Tdiags,'krylov_vs_h.png');
diags.plot_oversolve_vs_t(Tdiags,'oversolve_vs_t.png');
diags.plot_mass_vs_t(Tdiags,'mass_vs_t.png');
diags.plot_mass_vs_h(Tdiags,'mass_vs_h.png');

# print solve statistics to screen
diags.etest_stats(Tdiags, sys.stdout);
diags.solver_stats(Tdiags, sys.stdout);

# write the Tdiags structure to screen
#print 'Internal representation of solver diagnostics structure:'
#diags.write_diags(Tdiags);

##### end of script #####
