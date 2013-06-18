# matplotlib-based plotting and analysis utilities for ARKode solver diagnostics
# Daniel R. Reynolds
# SMU Mathematics


#### Data Structures ####

class ErrorTest:
    """ For each error test performed, we keep track of the time step """
    """ index, the time step size, the estimate of the local error    """
    """ error, and a flag denoting whether the error test failed.     """
    def __init__(self, step, h, dsm):
        self.index    = step;
        self.h        = h;
        if (dsm <= 0.0):
            self.estimate = 1.e-10
        else:
            self.estimate = dsm;
        self.errfail  = 0;
        if (dsm > 1.0):
            self.errfail = 1;
    def Write(self):
        print '  ErrorTest: index =',self.index,', h =',self.h,', estimate =',self.estimate,', errfail =',self.errfail

##########
class AdaptH:
    """ For each time step adaptivity calculation, we keep track of the """
    """ biased error history array (eh0, eh1 and eh2), the raw accuracy """
    """ and stability time step estimates prior to placing bounds and   """
    """ safety factors on them (h_accuracy0 and h_stability0), the same """
    """ estimates after bounds and safety factors have been enforced    """
    """ (h_accuracy1 and h_stability1), a flag denoting whether the     """
    """ step was limited by the stability estimate (stab_restrict),     """
    """ and the resulting time step growth factor (eta).                """
    def __init__(self, eh0, eh1, eh2, hh0, hh1, hh2, ha0, hs0, ha1, hs1, eta):
        self.eh0 = eh0;
        self.eh1 = eh1;
        self.eh2 = eh2;
        self.hh0 = hh0;
        self.hh1 = hh1;
        self.hh2 = hh2;
        self.h_accuracy0 = ha0;
        self.h_accuracy1 = ha1;
        self.h_stability0 = hs0;
        self.h_stability1 = hs1;
        self.eta = eta;
        if (hs0 < ha0):
            self.stab_restrict = 1;
        else:
            self.stab_restrict = 0;
    def Write(self):
        print '  AdaptH: errhist =',self.eh0,self.eh1,self.eh2,', stephist =',self.hh0,self.hh1,self.hh2,', ha0 =',self.h_accuracy0,', hs0 =',self.h_stability0,', ha1 =',self.h_accuracy1,', hs1 =',self.h_stability1,', stabrestrict =',self.stabrestrict,', eta =',self.eta

##########
class NewtonStep:
    """ Each Newton step holds the step iteration, the WRMS norm of """
    """ the Newton update (delta), and the convergence test (dcon). """
    def __init__(self, iter, delta, dcon):
        self.iter = iter;
        self.delta = delta;
        self.dcon  = dcon;
    def Write(self):
        print '      Nstep: iter =',self.iter,', delta =',self.delta,', dcon =',self.dcon

##########
class NewtonSolve:
    """ Each Newton solve holds the total number of Newton steps taken """
    """ (iters), an array of the Newton steps themselves (steps, of    """
    """ type NewtonStep), and a flag denoting whether the solve failed """
    """ (nonconv).                                                     """
    def __init__(self):
        self.steps   = [];    # container for NewtonStep objects
        self.iters   = -1;    # post-processing statistics
        self.nonconv = -1;
        self.dcon    = 0.0;
    def AddStep(self, NewtStep):
        self.steps.append(NewtStep);
    def Cleanup(self):
        self.iters = 0;
        for i in range(len(self.steps)):
            self.iters += 1;
        self.dcon = self.steps[-1].dcon;
        if (self.dcon < 1.0):
            self.nonconv = 0;
        else:
            self.nonconv = 1;
    def Write(self):
        print '    Nsolve: iters =',self.iters,', dcon =',self.dcon,', nonconv =',self.nonconv
        for i in range(len(self.steps)):
            self.steps[i].Write();

##########
class StageStep:
    """ For every RK stage of every time step, we keep track of the time  """
    """ step index (step), the time step size (h), the stage index        """
    """ (stage), the stage time (tn), each Newton solve used for          """
    """ calculation of the stage solution (NewtSolves -- we can have more """
    """ than one if the first fails and lsetup is called), the total      """
    """ number of Newton iterations required for convergence (NewtIters), """
    """ and the number of lsetups that occured within the step (lsetups). """
    def __init__(self, step, h, stage, tn):
        self.step  = step;
        self.h     = h;
        self.stage = stage;
        self.tn    = tn;
        self.NewtSolves = [];
        self.NewtIters = 0;
        self.NewtDcon = 0.0;
        self.lsetups = 0;
        self.conv_fails = 0;
    def AddLSetup(self):
        self.lsetups += 1;
    def AddNewton(self, NStep):
        if (NStep.iter == 0):    # create empty solve object
            self.NewtSolves.append(NewtonSolve());
        self.NewtSolves[-1].AddStep(NStep);
    def Cleanup(self):
        self.NewtIters  = 0;
        self.conv_fails = 0;
        for i in range(len(self.NewtSolves)):
            self.NewtSolves[i].Cleanup();
            self.NewtIters  += self.NewtSolves[i].iters;
            self.NewtDcon    = self.NewtSolves[i].dcon;
            self.conv_fails += self.NewtSolves[i].nonconv;
    def Write(self):
        print '  StageStep: step =',self.step,', h =',self.h,', stage =',self.stage,', tn =',self.tn,', conv_fails =',self.conv_fails,', lsetups =',self.lsetups,', NewtDcon =',self.NewtDcon,', NewtIters =',self.NewtIters
        for i in range(len(self.NewtSolves)):
            self.NewtSolves[i].Write();

##########
class TimeStep:
    """ For every time step, we keep track of the time step index (step),  """
    """ the different time step sizes that were attempted (h_attempts --   """
    """ typically only one, unless convergence or error test failures      """
    """ occur), the final successful time step size (h_final), and we      """
    """ store all StageSteps that comprised the time step.                 """
    def __init__(self):
        self.StageSteps = [];
        self.h_attempts = [];
        self.step       = -1;
        self.tn         = -1.0;
        self.h_final    = -1.0;
        self.NewtIters  = -1;
        self.lsetups    = -1;
        self.err_fails  =  0;
        self.conv_fails = -1;
    def AddStage(self, Stage):
        self.step = Stage.step
        if (self.h_final != Stage.h):
            self.h_final = Stage.h;
            self.h_attempts.append(Stage.h);
        self.StageSteps.append(Stage);
        self.tn = max(self.tn,Stage.tn);   # store maximum stage time in tn
    def AddErrorTest(self, ETest):
        self.ErrTest = ETest;
        self.err_fails += ETest.errfail;
    def AddHAdapt(self, HAdapt):
        self.HAdapt = HAdapt;
    def AddLSetup(self, stage):
        self.StageSteps[stage].AddLSetup();
    def AddNewton(self, stage, NStep):
        self.StageSteps[stage].AddNewton(NStep);
    def Cleanup(self):
        self.NewtIters  = 0;
        self.lsetups    = 0;
        self.conv_fails = 0;
        for i in range(len(self.StageSteps)):
            self.StageSteps[i].Cleanup();
            self.NewtIters += self.StageSteps[i].NewtIters;
            self.lsetups += self.StageSteps[i].lsetups;
            self.conv_fails += self.StageSteps[i].conv_fails;
    def Write(self):
        print 'TimeStep: step =',self.step,', tn =',self.tn,', h_attempts =',self.h_attempts,', h_final =',self.h_final,', NewtIters =',self.NewtIters,', lsetups =',self.lsetups,', err_fails =',self.err_fails,', conv_fails =',self.conv_fails
        for i in range(len(self.StageSteps)):
            self.StageSteps[i].Write();

#### Utility functions ####

def load_line(line):
    """ This routine parses a line of the diagnostics output file to  """
    """ determine what type of data it contains (an RK stage step, an """
    """ lsetup, a Newton iteration, an error test, or a time step     """
    """ adaptivity calculation), and creates the relevant object for  """
    """ that data line.  Each of these output types are indexed by a  """
    """ specific linetype for use by the calling routine.             """
    """                                                               """
    """ The function returns [linetype, entry].                       """
    import shlex
    txt = shlex.split(line)
    if ("step" in txt):
        linetype = 0;
        step  = int(txt[1]);
        h     = float(txt[2]);
        stage = int(txt[3]);
        tn    = float(txt[4]);
        entry = StageStep(step, h, stage, tn);
    elif ("lsetup" in txt):
        linetype = 1;
        entry = 0;
    elif ("newt" in txt):
        linetype = 2;
        iter  = int(txt[1]);
        delta = float(txt[2]);
        dcon  = float(txt[3]);
        entry = NewtonStep(iter, delta, dcon);
    elif ("etest" in txt):
        linetype = 3;
        step  = int(txt[1]);
        h     = float(txt[2]);
        dsm   = float(txt[3]);
        entry = ErrorTest(step, h, dsm);
    elif ("adapt" in txt):
        linetype = 4;
        eh0 = float(txt[1]);
        eh1 = float(txt[2]);
        eh2 = float(txt[3]);
        hh0 = float(txt[4]);
        hh1 = float(txt[5]);
        hh2 = float(txt[6]);
        ha0 = float(txt[7]);
        hs0 = float(txt[8]);
        ha1 = float(txt[9]);
        hs1 = float(txt[10]);
        eta = float(txt[11]);
        entry = AdaptH(eh0, eh1, eh2, hh0, hh1, hh2, ha0, hs0, ha1, hs1, eta);
    else:
        linetype = -1;
        entry = 0;
    return [linetype, entry]

##########
def load_diags(fname):
    """ This routine opens the diagnostics output file, loads all lines """
    """ to create an array of TimeSteps with all of the relevant data   """
    """ contained therein.                                              """
    f = open(fname)
    step  = -1;
    stage = -1;
    TimeSteps = [];
    for line in f:
        linetype, entry = load_line(line);
        if (linetype == 0):   # stage step
            if (entry.step != step):   # new step
                step = entry.step;
                TimeSteps.append(TimeStep());
            TimeSteps[step].AddStage(entry);
            stage = entry.stage;
        elif (linetype == 1):   # lsetup
            TimeSteps[step].AddLSetup(stage);
        elif (linetype == 2):   # newton step
            TimeSteps[step].AddNewton(stage,entry)
        elif (linetype == 3):   # error test
            TimeSteps[step].AddErrorTest(entry);
        elif (linetype == 4):   # h adaptivity
            TimeSteps[step].AddHAdapt(entry);
    f.close()
    for i in range(len(TimeSteps)):
        TimeSteps[i].Cleanup();
    return TimeSteps

##########
def write_diags(TimeSteps):
    """ This routine takes in the array of TimeSteps (returned from  """
    """ load_diags), and writes out the internal representation of   """
    """ the time step history to stdout.                             """
    for i in range(len(TimeSteps)):
        print '  '
        TimeSteps[i].Write();

##########
def plot_h_vs_t(TimeSteps,fname):
    """ This routine takes in the array of TimeSteps (returned from  """
    """ load_diags), and plots the time step sizes h as a function   """
    """ of the simulation time t.  Failed time steps are marked on   """
    """ the plot with either a red X or green O, where X corresponds """
    """ to an error test failure, and an O corresponds to a Newton   """
    """ solver convergence failure.                                  """
    """                                                              """
    """ The resulting plot is stored in the file <fname>, that       """
    """ should include an extension appropriate for the matplotlib   """
    """ 'savefig' command.                                           """
    import pylab as plt
    import numpy as np
    hvals  = [];
    tvals  = [];
    EfailH = [];
    EfailT = [];
    CfailH = [];
    CfailT = [];
    for istep in range(len(TimeSteps)):

        # store successful step size and time
        hvals.append(TimeSteps[istep].h_final);
        tvals.append(TimeSteps[istep].tn);

        # account for convergence failures and error test failures
        if (TimeSteps[istep].err_fails > 0):
            EfailH.append(TimeSteps[istep].h_attempts[0]);
            EfailT.append(TimeSteps[istep].tn);
        if (TimeSteps[istep].conv_fails > 0):
            CfailH.append(TimeSteps[istep].h_attempts[0]);
            CfailT.append(TimeSteps[istep].tn);

    # convert data to numpy arrays
    h = np.array(hvals);
    t = np.array(tvals);
    eh = np.array(EfailH);
    et = np.array(EfailT);
    ch = np.array(CfailH);
    ct = np.array(CfailT);

    # generate plot
    plt.figure()
    plt.plot(t,h,'b-')
    plt.plot(et,eh,'rx')
    plt.plot(ct,ch,'go')
    plt.xlabel('time')
    plt.ylabel('step size')
    plt.title('Step size versus time')
    plt.legend(('successful','error fails','conv. fails'), loc='upper left', shadow=True)
    plt.grid()
    plt.savefig(fname)


##########
def plot_h_vs_iter(TimeSteps,fname):
    """ This routine takes in the array of TimeSteps (returned from  """
    """ load_diags), and plots the time step sizes h as a function   """
    """ of the time step iteration index.  Failed time steps are     """
    """ marked on the plot with either a red X or green O, where X   """
    """ corresponds to an error test failure, and an O corresponds   """
    """ to a Newton solver convergence failure.                      """
    """                                                              """
    """ The resulting plot is stored in the file <fname>, that       """
    """ should include an extension appropriate for the matplotlib   """
    """ 'savefig' command.                                           """
    import pylab as plt
    import numpy as np
    hvals  = [];
    ivals  = [];
    EfailH = [];
    EfailI = [];
    CfailH = [];
    CfailI = [];
    for istep in range(len(TimeSteps)):

        # store successful step size and time
        hvals.append(TimeSteps[istep].h_final);
        ivals.append(istep);

        # account for convergence failures and error test failures
        if (TimeSteps[istep].err_fails > 0):
            EfailH.append(TimeSteps[istep].h_attempts[0]);
            EfailI.append(istep);
        if (TimeSteps[istep].conv_fails > 0):
            CfailH.append(TimeSteps[istep].h_attempts[0]);
            CfailI.append(istep);

    # convert data to numpy arrays
    h = np.array(hvals);
    I = np.array(ivals);
    eh = np.array(EfailH);
    eI = np.array(EfailI);
    ch = np.array(CfailH);
    cI = np.array(CfailI);

    # generate plot
    plt.figure()
    plt.plot(I,h,'b-')
    plt.plot(eI,eh,'rx')
    plt.plot(cI,ch,'go')
    plt.xlabel('time step')
    plt.ylabel('step size')
    plt.title('Step size versus iteration')
    plt.legend(('successful','error fails','conv. fails'), loc='upper left', shadow=True)
    plt.grid()
    plt.savefig(fname)


##########
def plot_work_vs_t(TimeSteps,fname):
    """ This routine takes in the array of TimeSteps (returned from  """
    """ load_diags), and plots the total number of Newton iterations """
    """ as a function of the simulation time t.  Solves that failed  """
    """ the nonlinear convergence test are marked on the plot with   """
    """ a red X, and lsetup calls are marked on the plot with a      """
    """ green O.                                                     """
    """                                                              """
    """ The resulting plot is stored in the file <fname>, that       """
    """ should include an extension appropriate for the matplotlib   """
    """ 'savefig' command.                                           """
    import pylab as plt
    import numpy as np
    Nvals   = [];
    tvals   = [];
    LsetupN = [];
    LsetupT = [];
    CfailN  = [];
    CfailT  = [];
    for istep in range(len(TimeSteps)):
        
        # store Newton iterations and time
        Nvals.append(TimeSteps[istep].NewtIters);
        tvals.append(TimeSteps[istep].tn);

        # account for convergence failures and lsetups
        if (TimeSteps[istep].conv_fails > 0):
            CfailN.append(TimeSteps[istep].NewtIters);
            CfailT.append(TimeSteps[istep].tn);
        if (TimeSteps[istep].lsetups > 0):
            LsetupN.append(TimeSteps[istep].NewtIters);
            LsetupT.append(TimeSteps[istep].tn);

    # convert data to numpy arrays
    N = np.array(Nvals);
    t = np.array(tvals);
    lN = np.array(LsetupN);
    lt = np.array(LsetupT);
    cN = np.array(CfailN);
    ct = np.array(CfailT);

    # generate plot
    plt.figure()
    plt.plot(t,N,'b-')
    plt.plot(ct,cN,'rx')
    plt.plot(lt,lN,'go')
    plt.xlabel('time')
    plt.ylabel('Newton iters')
    plt.title('Newton iterations versus time')
    plt.legend(('successful','conv. fails','lsetups'), loc='upper left', shadow=True)
    plt.grid()
    plt.savefig(fname)


##########
def plot_work_vs_h(TimeSteps,fname):
    """ This routine takes in the array of TimeSteps (returned from  """
    """ load_diags), and plots the total number of Newton iterations """
    """ as a function of the time step size h.  Solves that failed   """
    """ the nonlinear convergence test are marked on the plot with   """
    """ a red X, whereas all other data is plotted with a blue dot.  """
    """                                                              """
    """ The resulting plot is stored in the file <fname>, that       """
    """ should include an extension appropriate for the matplotlib   """
    """ 'savefig' command.                                           """
    import pylab as plt
    import numpy as np
    Nvals   = [];
    hvals   = [];
    LsetupN = [];
    LsetupH = [];
    CfailN  = [];
    CfailH  = [];
    for istep in range(len(TimeSteps)):
        
        # store Newton iterations and time
        Nvals.append(TimeSteps[istep].NewtIters);
        hvals.append(TimeSteps[istep].h_final);

        # account for convergence failures
        if (TimeSteps[istep].conv_fails > 0):
            CfailH.append(TimeSteps[istep].h_attempts[0]);
            CfailN.append(TimeSteps[istep].NewtIters);

    # convert data to numpy arrays
    N = np.array(Nvals);
    h = np.array(hvals);
    cN = np.array(CfailN);
    ch = np.array(CfailH);

    # generate plot
    plt.figure()
    plt.plot(h,N,'b.')
    plt.plot(ch,cN,'rx')
    plt.xlabel('step size')
    plt.ylabel('Newton iters')
    plt.title('Newton iterations versus step size')
    plt.legend(('successful','conv. fails'), loc='upper left', shadow=True)
    plt.grid()
    plt.savefig(fname)


##########
def plot_oversolve_vs_t(TimeSteps,fname):
    """ This routine takes in the array of TimeSteps (returned from  """
    """ load_diags), and plots the oversolve as a function of the    """
    """ simulation time t. We cap the computed oversolve value at    """
    """ 1000 to get more intuitive plots since first few time steps  """
    """ are purposefully too small.                                  """
    """                                                              """
    """ The resulting plot is stored in the file <fname>, that       """
    """ should include an extension appropriate for the matplotlib   """
    """ 'savefig' command.                                           """
    import pylab as plt
    import numpy as np
    Ovals = [];
    tvals = [];
    for istep in range(len(TimeSteps)):

        # store oversolve and time
        Ovals.append(min(1e3,1.0 / TimeSteps[istep].ErrTest.estimate));
        tvals.append(TimeSteps[istep].tn);

    # convert data to numpy arrays
    O = np.array(Ovals);
    t = np.array(tvals);

    # generate plot
    plt.figure()
    plt.semilogy(t,O,'b-')
    plt.xlabel('time')
    plt.ylabel('oversolve')
    plt.title('Oversolve versus time')
    plt.grid()
    plt.savefig(fname)


##########
def etest_stats(TimeSteps,fptr):
    """ This routine takes in the array of TimeSteps (returned from  """
    """ load_diags), and computes statistics on how well the time    """
    """ step adaptivity estimations predicted step values that met   """
    """ the accuracy requirements.                                   """
    """                                                              """
    """ Oversolves: we define 'oversolveX' as a step in which the    """
    """ resulting time step accuracy was X times more accurate than  """
    """ required.                                                    """
    """                                                              """
    """ Note: we ignore steps immediately following a convergence or """
    """ error test failure, since etamax is bounded above by 1.      """
    """                                                              """
    """ The resulting data is appended to the stream corresponding   """
    """ to fptr (either the result of 'open' or sys.stdout).         """
    import numpy as np
    oversolve2   = 0;
    oversolve4   = 0;
    oversolve10  = 0;
    oversolve30  = 0;
    oversolve100 = 0;
    oversolve_max   = 0.0;
    oversolve_min   = 1.0e200;
    oversolves = [];
    for istep in range(len(TimeSteps)):
        # if this or previous step had an error test failure, skip oversolve results
        ignore = 0;
        if (istep > 0):
            if( (TimeSteps[istep].err_fails > 0) or (TimeSteps[istep-1].err_fails > 0)):
                ignore = 1;
        else:
            if(TimeSteps[istep].err_fails > 0):
                ignore = 1;
            
        if (ignore == 0):
            over = 1.0 / TimeSteps[istep].ErrTest.estimate;
            oversolves.append(over);
            if (over > 100.0):
                oversolve100 += 1;
            elif (over > 30.0):
                oversolve30 += 1;
            elif (over > 10.0):
                oversolve10 += 1;
            elif (over > 4.0):
                oversolve4 += 1;
            elif (over > 2.0):
                oversolve2 += 1;
            oversolve_max = max(oversolve_max, over)
            oversolve_min = min(oversolve_min, over)

    # generate means and percentages
    ov = np.array(oversolves);
    oversolve_amean = np.mean(ov);
    oversolve_gmean = np.prod(np.power(ov,1.0/len(ov)));
    oversolve2   = 1.0*oversolve2/len(TimeSteps);
    oversolve4   = 1.0*oversolve4/len(TimeSteps);
    oversolve10  = 1.0*oversolve10/len(TimeSteps);
    oversolve30  = 1.0*oversolve30/len(TimeSteps);
    oversolve100 = 1.0*oversolve100/len(TimeSteps);

    fptr.write("\n")
    fptr.write("Simulation took %i steps\n" % (len(TimeSteps)))
    fptr.write("  Oversolve fractions:\n")
    fptr.write("        2x fraction = %g\n" % (oversolve2))
    fptr.write("        4x fraction = %g\n" % (oversolve4))
    fptr.write("       10x fraction = %g\n" % (oversolve10))
    fptr.write("       30x fraction = %g\n" % (oversolve30))
    fptr.write("      100x fraction = %g\n" % (oversolve100))
    fptr.write("  Oversolve extrema:\n")
    fptr.write("       min = %g\n" % (oversolve_min))
    fptr.write("       max = %g\n" % (oversolve_max))
    fptr.write("  Oversolve means:\n")
    fptr.write("       arithmetic = %g\n" % (oversolve_amean))
    fptr.write("        geometric = %g\n" % (oversolve_gmean))
    fptr.write("\n")
    
##########
def solver_stats(TimeSteps,fptr):
    """ This routine takes in the array of TimeSteps (returned from  """
    """ load_diags), and computes statistics on how the solver       """
    """ performed.                                                   """
    """                                                              """
    """ The resulting data is appended to the stream corresponding   """
    """ to fptr (either the result of 'open' or sys.stdout).         """
    import numpy as np
    newtits = 0;
    stages  = 0;
    newtfails = 0;
    lsetups = 0;
    dcon_final = [];
    for istep in range(len(TimeSteps)):
        newtits += TimeSteps[istep].NewtIters;
        stages += len(TimeSteps[istep].StageSteps);
        newtfails += TimeSteps[istep].conv_fails;
        lsetups += TimeSteps[istep].lsetups;
        for istage in range(len(TimeSteps[istep].StageSteps)):
            dcon_final.append(TimeSteps[istep].StageSteps[istage].NewtDcon);

    # generate means and percentages
    dcon = np.array(dcon_final);
    if (lsetups == 0):
        newt_per_lsetup = 0;
    else:
        newt_per_lsetup = (1.0*newtits/lsetups)

    fptr.write("\n")
    fptr.write("Simulation took %i steps with %i stages\n" % (len(TimeSteps),stages))
    fptr.write("    Avg Newton / step  = %g\n" % (1.0*newtits/len(TimeSteps)))
    fptr.write("    Avg Newton / stage = %g\n" % (1.0*newtits/stages))
    fptr.write("    Total Newton failures = %i\n" % (newtfails))
    fptr.write("    Average Newton / lsetup = %g\n" % (newt_per_lsetup))
    fptr.write("  Newton convergence tests (dcon):\n")
    fptr.write("       min = %g\n" % (np.min(dcon)))
    fptr.write("       max = %g\n" % (np.max(dcon)))
    fptr.write("       avg = %g\n" % (np.mean(dcon)))
    fptr.write("     gmean = %g\n" % (np.prod(np.power(dcon,1.0/len(dcon)))))
    fptr.write("\n")
    
   


##########

