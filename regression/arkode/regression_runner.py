#!/usr/bin/env python
#------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
#------------------------------------------------------------
# Copyright (c) 2014, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
#------------------------------------------------------------
# Function to perform regression tests on ARKODE solvers.  If
# a file local_standards.dat exists, it will load from those; 
# otherwise it will load from the file 
# regression_standards.dat.

# main routine
def main():

    # imports
    import sys
    import os
    import os.path
    from time import time
    import arkode_tools as ark
    import pickle
    from optparse import OptionParser
    #from argparse import ArgumentParser

    # parse arguments from command-line
    parser = OptionParser(usage="usage: %prog [options]")
    parser.add_option("-s", "--strict", action="store_true", 
                      dest="strict", default=True,
                      help="Run in strict mode.  This checks for identical"
                      + " solver statistics, and very similar run times as"
                      + " reference values. [default]")
    parser.add_option("-l", "--lenient", action="store_false", 
                      dest="strict", default=True,
                      help="Run in lenient mode (opposite of strict).  This checks" 
                      + " for similar results to reference values; useful when"
                      + " testing on new machines where architecture differences"
                      + " may slightly affect results.")
    parser.add_option("-f", "--full", action="store_false", 
                      dest="quick", default=False,
                      help="Runs the full suite of regression tests stored in the"
                      + " reference set.  These tests may take some time to finish,"
                      + " though an estimate from the reference results is displayed"
                      + " at start. [default]")
    parser.add_option("-q", "--quick", action="store_true", 
                      dest="quick", default=False,
                      help="Runs an abbreviated set of tests from the reference set"
                      + " (opposite of full).  These should encompass most of the"
                      + " tests from the full suite, but are downselected to choose"
                      + " only those faster than a specified threshold.")
    (options, args) = parser.parse_args()

    # access command-line variables
    strict = options.strict
    quick = options.quick

    ## parse arguments from command-line
    #parser = ArgumentParser(usage="usage: %prog [options]")
    #parser.add_argument("-s", "--strict", action="store_true", default=True,
    #                    help="Run in strict mode")
    #parser.add_argument("-l", "--lenient", action="store_false", default=False,
    #                    help="Run in lenient mode")
    #parser.add_argument("-f", "--full", action="store_false", default=True, 
    #                    help="Runs the full suite of regression tests")
    #parser.add_argument("-q", "--quick", action="store_true", default=False,
    #                    help="Runs an abbreviated set of regression tests")
    #args = parser.parse_args()

    ## access command-line variables
    #strict = True
    #if args.lenient:
    #    strict = False
    #quick = False
    #if args.quick:
    #    quick = True

    # set threshold time for 'quick' tests
    tquick = 1.0

    # load standards dictionary from disk
    if os.path.isfile("local_standards.dat") and os.access("local_standards.dat", os.R_OK):
        print '\nReading reference results from: local_standards.dat'
        gold_standard = pickle.load( open( "local_standards.dat", "rb" ) )
    else:
        print '\nReading reference results from: regression_standards.dat'
        gold_standard = pickle.load( open( "regression_standards.dat", "rb" ) )


    # read in list of tests
    AllTests = ark.ReadTests('regression_tests.txt')

    print '\nRunning regression tests:'
    if strict:
        print '  strict error-checking'
    else:
        print '  lenient error-checking'
    if quick:
        print '  quick suite'
    else:
        print '  full suite'

    # determine total number of tests to run
    numtests = 0
    time_estimate = 0.0
    for test in AllTests:
        if ((not quick) or (gold_standard[test.name].runtime < tquick)):
            numtests += 1
            time_estimate += gold_standard[test.name].runtime
    sys.stdout.write("\nRunning %i tests, estimated runtime = %g seconds\n\n" 
                     % (numtests, time_estimate))
            

    # run requested tests in list
    tstart = time()
    numsuccess = 0
    numfailed = 0
    itest = 0
    failedtests = []
    for test in AllTests:

        # with quick suite, check test length
        if ((not quick) or (gold_standard[test.name].runtime < tquick)):

            # begin test line output
            itest += 1
            sys.stdout.write("%4i/%i:  %30s" % (itest, numtests, test.name))

            # run test
            stats = test.Run(0)

            # check results against gold standard
            result = stats.Compare(gold_standard[test.name], strict=strict)

            # output results to stdout (use color if available)
            if os.environ["TERM"] == 'xterm':
                if (result == 1):
                    numsuccess += 1
                    sys.stdout.write(" \033[92m pass  \033[94m(steps: %6i,  runtime: %.2g s)\033[0m\n" 
                                     % (stats.nsteps, stats.runtime))
                elif (result == 2):
                    numsuccess += 1
                    sys.stdout.write(" \033[93m warn  \033[94m(steps: %6i,  runtime: %.2g s)\033[0m\n" 
                                     % (stats.nsteps, stats.runtime))
                else:
                    numfailed += 1
                    failedtests.append(test.name)
                    sys.stdout.write(" \033[91m fail  \033[94m(runtime: %.2g s)\033[0m\n" % (stats.runtime))
            else:
                if (result == 1):
                    numsuccess += 1
                    sys.stdout.write("  pass  (steps: %6i,  runtime: %.2g s)\n" 
                                     % (stats.nsteps, stats.runtime))
                elif (result == 2):
                    numsuccess += 1
                    sys.stdout.write("  warn  (steps: %6i,  runtime: %.2g s)\n" 
                                     % (stats.nsteps, stats.runtime))
                else:
                    numfailed += 1
                    failedtests.append(test.name)
                    sys.stdout.write("  fail  (runtime: %.2g s)\n" % (stats.runtime))

            
    tend = time()

    # print final statistics
    sys.stdout.write("\nRan %i total regression tests:\n" % (numtests) )
    sys.stdout.write("  Passed: %i\n" % (numsuccess))
    sys.stdout.write("  Failed: %i\n" % (numfailed))
    sys.stdout.write("Total test time: %g seconds\n" % (tend-tstart))
    if (numfailed > 0):
        sys.stdout.write("Failed tests were:\n")
        for test in failedtests:
            sys.stdout.write("  %s\n" % (test))
        sys.exit(1)



# just run the main routine
if __name__ == '__main__':
    main()


#### end of script ####
