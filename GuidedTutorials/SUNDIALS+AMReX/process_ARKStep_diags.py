#!/usr/bin/env python3
#------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
#------------------------------------------------------------
# Copyright (c) 2019, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
#------------------------------------------------------------
# analysis script for ARKODE solver statistics

# imports
import sys
import diag_tools as diags

# output error if no diagnostics output file is provided
if (len(sys.argv) == 1):
    sys.exit("\nCalling syntax error: filename containing diagnostics must be provided, e.g.\n    $ process_ARKStep_diags.py fname.txt\n")

# load the diagnostics output file (passed on the command line)
fname = sys.argv[len(sys.argv)-1]

# load the time step data
Tdiags = diags.load_diags(fname)

# generate plots
print('\nGenerating plots h_vs_t.png, h_vs_iter.png, and oversolve_vs_t.png')
diags.plot_h_vs_t(Tdiags,'h_vs_t.png')
diags.plot_h_vs_tstep(Tdiags,'h_vs_iter.png')
diags.plot_oversolve_vs_t(Tdiags,'oversolve_vs_t.png')

# print solve statistics to screen
diags.etest_stats(Tdiags, sys.stdout)

##### end of script #####
