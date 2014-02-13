#------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
#------------------------------------------------------------
# Copyright (c) 2013, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
#------------------------------------------------------------
# matplotlib-based plotting and analysis utilities for ARKode solver diagnostics


#### Data Structures ####

class ErrorTest:
    """ An ErrorTest object stores, for each error test performed: """
    """    index    -- the time step index """
    """    h        -- the time step size """
    """    estimate -- the estimate of the local error """
    """    errfail  -- a flag denoting whether the error test failed """
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
    """ An AdaptH object stores, for each time step adaptivity calculation: """
    """    eh[0,1,2]        -- the biased error history array """
    """    hh[0,1,2]        -- the time step history array """
    """    h_accuracy[0,1]  -- the accuracy step estimates, before  """
    """                        and after applying step size bounds """
    """    h_stability[0,1] -- the stability step estimates, before """
    """                        and after applying cfl & stability bounds """
    """    stab_restrict    -- flag whether step was stability-limited """
    """    eta              -- the final time step growth factor """
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
class KrylovSolve:
    """ A KrylovSolve object stores, for each Krylov solve: """
    """    bnorm   -- the norm of the right-hand side"""
    """    resnorm -- the final linear residual norm"""
    """    iters   -- the total number of steps taken """
    """    psolves -- the total number of preconditioner solves taken """
    def __init__(self, bnorm, resnorm, iters, psolves):
        self.bnorm   = bnorm;    # container for KrylovSolve objects
        self.resnorm = resnorm;  # post-processing statistics
        self.iters   = iters;
        self.psolves = psolves;
    def Write(self):
        print '      KSolve: iters =',self.iters,', psolves =',self.psolves,', bnorm =',self.bnorm,', resnorm =',self.resnorm

##########
class MassSolve:
    """ A MassSolve object stores, for each Krylov mass-matrix solve: """
    """    resnorm -- the final linear residual norm"""
    """    iters   -- the total number of steps taken """
    """    psolves -- the total number of preconditioner solves taken """
    def __init__(self, resnorm, iters, psolves):
        self.resnorm = resnorm;  # container for MassSolve objects
        self.iters   = iters;    # post-processing statistics
        self.psolves = psolves;
    def Write(self):
        print '      Msolve: iters =',self.iters,', psolves =',self.psolves,', resnorm =',self.resnorm

##########
class NonlinStep:
    """ A NonlinStep object stores, for each nonlinear iteration: """
    """    iter  -- the step iteration """
    """    delta -- the WRMS norm of the nonlinear update """
    """    dcon  -- the convergence test value"""
    def __init__(self, iter, delta, dcon):
        self.iter = iter;
        self.delta = delta;
        self.dcon  = dcon;
    def Write(self):
        print '      Nstep: iter =',self.iter,', delta =',self.delta,', dcon =',self.dcon

##########
class NonlinSolve:
    """ A NonlinSolve object stores, for each nonlinear solve: """
    """    iters    -- the total number of steps taken """
    """    steps    -- array of nonlinear steps taken (of type NonlinStep) """
    """    KSolves  -- array of linear solves taken (of type KrylovSolve) """
    """    liniters -- the total number of linear steps taken """ 
    """    nonconv  -- a flag denoting whether the solve failed """
    """    dcon     -- the final convergence test value """
    def __init__(self):
        self.steps    = [];
        self.KSolves  = [];
        self.iters    = -1;
        self.liniters = -1;
        self.nonconv  = -1;
        self.dcon     = 0.0;
    def AddStep(self, NonlinStep):
        self.steps.append(NonlinStep);
    def AddKrylov(self, KrylovSolve):
        self.KSolves.append(KrylovSolve);
    def Cleanup(self):
        self.iters = 0;
        self.liniters = 0;
        for i in range(len(self.steps)):
            self.iters += 1;
        for i in range(len(self.KSolves)):
            self.liniters += self.KSolves[i].iters;
        self.dcon = self.steps[-1].dcon;
        if (self.dcon < 1.0):
            self.nonconv = 0;
        else:
            self.nonconv = 1;
    def Write(self):
        print '    Nsolve: iters =',self.iters,', liniters =',self.liniters,', dcon =',self.dcon,', nonconv =',self.nonconv
        for i in range(len(self.steps)):
            self.steps[i].Write();
        for i in range(len(self.KSolves)):
            self.KSolves[i].Write();

##########
class StageStep:
    """ A StageStep object stores, for each RK stage of every time step: """
    """    step         -- the time step index """
    """    h            -- the time step size """
    """    stage        -- the stage index """
    """    tn           -- the stage time """
    """    NonlinSolves -- array of nonlinear solves used for stage """
    """                    solution (of type NonlinSolves); we can have  """
    """                    multiple if the first fails and lsetup is called) """
    """    NonlinIters  -- the total number of nonlinear iterations """
    """    NonlinDcon   -- the final nonlinear convergence test value """
    """    lsetups      -- the number of lsetups that occured within the step """
    """    KrylovIters  -- the total number of Krylov iterations """
    """    conv_fails   -- the number of nonlinear convergence failures """
    def __init__(self, step, h, stage, tn):
        self.step  = step;
        self.h     = h;
        self.stage = stage;
        self.tn    = tn;
        self.NonlinSolves = [];
        self.NonlinSolves.append(NonlinSolve());  # create first nonlinear solver
        self.NSolveData = False;            # set flag that solve data not set
        self.NonlinIters = -1;
        self.KrylovIters = -1;
        self.NonlinDcon = 0.0;
        self.lsetups = 0;
        self.conv_fails = 0;
    def AddLSetup(self):
        self.lsetups += 1;
    def AddNonlin(self, NStep):
        if (NStep.iter == 0 and self.NSolveData): # add new solve if starting over
            self.NonlinSolves.append(NonlinSolve());
        self.NSolveData = True                # update flag
        self.NonlinSolves[-1].AddStep(NStep);
    def AddKrylov(self, KSolve):
        self.NonlinSolves[-1].AddKrylov(KSolve);  # add to newest nonlinear solver
    def Cleanup(self):
        self.NonlinIters  = 0;
        self.KrylovIters  = 0;
        self.conv_fails = 0;
        for i in range(len(self.NonlinSolves)):
            self.NonlinSolves[i].Cleanup();
            self.NonlinIters  += self.NonlinSolves[i].iters;
            self.KrylovIters  += self.NonlinSolves[i].liniters;
            self.NonlinDcon    = self.NonlinSolves[i].dcon;
            self.conv_fails += self.NonlinSolves[i].nonconv;
    def Write(self):
        print '  StageStep: step =',self.step,', h =',self.h,', stage =',self.stage,', tn =',self.tn,', conv_fails =',self.conv_fails,', lsetups =',self.lsetups,', NonlinDcon =',self.NonlinDcon,', NonlinIters =',self.NonlinIters,', KrylovIters =',self.KrylovIters
        for i in range(len(self.NonlinSolves)):
            self.NonlinSolves[i].Write();

##########
class TimeStep:
    """ A TimeStep object stores, for every time step: """
    """    StageSteps  -- array of StageStep objects comprising the step"""
    """    MassSolve   -- MassSolve object used to finish the step """
    """    h_attempts  -- array of step sizes attempted (typically only one, """
    """                   unless convergence or error failures occur) """
    """    step        -- the time step index """
    """    tn          -- maximum stage time in step """
    """    h_final     -- the final successful time step size """
    """    NonlinIters -- total nonlinear iters for step """
    """    lsetups     -- total lsetup calls in step """
    """    ErrTest     -- an ErrorTest object for the step """
    """    err_fails   -- total error test failures in step """
    """    conv_fails  -- total number of convergence failures in step """
    def __init__(self):
        self.StageSteps  = [];
        self.MassSolve   = 0;
        self.h_attempts  = [];
        self.step        = -1;
        self.tn          = -1.0;
        self.h_final     = -1.0;
        self.NonlinIters = -1;
        self.KrylovIters = -1;
        self.lsetups     = -1;
        self.err_fails   =  0;
        self.conv_fails  = -1;
    def AddStage(self, Stage):
        self.step = Stage.step
        if (self.h_final != Stage.h):
            self.h_final = Stage.h;
            self.h_attempts.append(Stage.h);
        self.StageSteps.append(Stage);
        self.tn = max(self.tn,Stage.tn);
    def AddErrorTest(self, ETest):
        self.ErrTest = ETest;
        self.err_fails += ETest.errfail;
    def AddHAdapt(self, HAdapt):
        self.HAdapt = HAdapt;
    def AddLSetup(self, stage):
        self.StageSteps[stage].AddLSetup();
    def AddNonlin(self, stage, NStep):
        self.StageSteps[stage].AddNonlin(NStep);
    def AddKrylov(self, stage, KSolve):
        self.StageSteps[stage].AddKrylov(KSolve);
    def AddMass(self, MSolve):
        self.MassSolve = MSolve;
    def Cleanup(self):
        self.NonlinIters = 0;
        self.KrylovIters = 0;
        self.lsetups     = 0;
        self.conv_fails  = 0;
        for i in range(len(self.StageSteps)):
            self.StageSteps[i].Cleanup();
            self.NonlinIters += self.StageSteps[i].NonlinIters;
            self.KrylovIters += self.StageSteps[i].KrylovIters;
            self.lsetups += self.StageSteps[i].lsetups;
            self.conv_fails += self.StageSteps[i].conv_fails;
    def Write(self):
        print 'TimeStep: step =',self.step,', tn =',self.tn,', h_attempts =',self.h_attempts,', h_final =',self.h_final,', NonlinIters =',self.NonlinIters,', KrylovIters =',self.KrylovIters,', lsetups =',self.lsetups,', err_fails =',self.err_fails,', conv_fails =',self.conv_fails
        for i in range(len(self.StageSteps)):
            self.StageSteps[i].Write();
        if (self.MassSolve != 0):
            self.MassSolve.Write();


#### Utility functions ####

def load_line(line):
    """ This routine parses a line of the diagnostics output file to  """
    """ determine what type of data it contains (an RK stage step, an """
    """ lsetup, a nonlinear iteration, an error test, or a time step  """
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
    elif (("newt" in txt) or ("fp" in txt)):
        linetype = 2;
        iter  = int(txt[1]);
        delta = float(txt[2]);
        dcon  = float(txt[3]);
        entry = NonlinStep(iter, delta, dcon);
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
    elif ("kry" in txt):
        linetype = 5;
        bnorm = float(txt[1]);
        resnorm = float(txt[2]);
        iters = int(txt[3]);
        psolves = int(txt[4]);
        entry = KrylovSolve(bnorm, resnorm, iters, psolves);
    elif ("mass" in txt):
        linetype = 6;
        resnorm = float(txt[1]);
        iters = int(txt[2]);
        psolves = int(txt[3]);
        entry = MassSolve(resnorm, iters, psolves);
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
        elif (linetype == 2):   # nonlinear step
            TimeSteps[step].AddNonlin(stage,entry)
        elif (linetype == 3):   # error test
            TimeSteps[step].AddErrorTest(entry);
        elif (linetype == 4):   # h adaptivity
            TimeSteps[step].AddHAdapt(entry);
        elif (linetype == 5):   # Krylov solve
            TimeSteps[step].AddKrylov(stage,entry);
        elif (linetype == 6):   # mass Krylov solve
            TimeSteps[step].AddMass(entry);
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
    """ to an error test failure, and an O corresponds to a          """
    """ nonlinear solver convergence failure.                        """
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
    plt.semilogy(t,h,'b-')
    plt.semilogy(et,eh,'rx')
    plt.semilogy(ct,ch,'go')
    plt.xlabel('time')
    plt.ylabel('step size')
    plt.title('Step size versus time')
    plt.legend(('successful','error fails','conv. fails'), loc='lower right', shadow=True)
    plt.grid()
    plt.savefig(fname)


##########
def plot_h_vs_iter(TimeSteps,fname):
    """ This routine takes in the array of TimeSteps (returned from  """
    """ load_diags), and plots the time step sizes h as a function   """
    """ of the time step iteration index.  Failed time steps are     """
    """ marked on the plot with either a red X or green O, where X   """
    """ corresponds to an error test failure, and an O corresponds   """
    """ to a nonlinear solver convergence failure.                   """
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
    plt.semilogy(I,h,'b-')
    plt.semilogy(eI,eh,'rx')
    plt.semilogy(cI,ch,'go')
    plt.xlabel('time step')
    plt.ylabel('step size')
    plt.title('Step size versus iteration')
    plt.legend(('successful','error fails','conv. fails'), loc='lower right', shadow=True)
    plt.grid()
    plt.savefig(fname)


##########
def plot_work_vs_t(TimeSteps,fname):
    """ This routine takes in the array of TimeSteps (returned from  """
    """ load_diags), and plots the total number of nonlinear         """
    """ iterations as a function of the simulation time t.  Solves   """
    """ that failed the nonlinear convergence test are marked on the """
    """ plot with a red X, and lsetup calls are marked on the plot   """
    """ with a green O.                                              """
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
        
        # store nonlinear iterations and time
        Nvals.append(TimeSteps[istep].NonlinIters);
        tvals.append(TimeSteps[istep].tn);

        # account for convergence failures and lsetups
        if (TimeSteps[istep].conv_fails > 0):
            CfailN.append(TimeSteps[istep].NonlinIters);
            CfailT.append(TimeSteps[istep].tn);
        if (TimeSteps[istep].lsetups > 0):
            LsetupN.append(TimeSteps[istep].NonlinIters);
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
    if (len(lN) > 0):
        plt.plot(lt,lN,'go')
    plt.xlabel('time')
    plt.ylabel('Nonlinear iters')
    plt.title('Nonlinear iterations per step versus time')
    if (len(lN) > 0):
        plt.legend(('successful','conv. fails','lsetups'), loc='lower right', shadow=True)
    else:
        plt.legend(('successful','conv. fails'), loc='lower right', shadow=True)
    plt.grid()
    plt.savefig(fname)


##########
def plot_work_vs_h(TimeSteps,fname):
    """ This routine takes in the array of TimeSteps (returned from  """
    """ load_diags), and plots the total number of nonlinear         """
    """ iterations as a function of the time step size h.  Solves    """
    """ that failed the nonlinear convergence test are marked on the """
    """ plot with a red X, whereas all other data is plotted with a  """
    """ blue dot.                                                    """
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
        
        # store nonlinear iterations and time
        Nvals.append(TimeSteps[istep].NonlinIters);
        hvals.append(TimeSteps[istep].h_final);

        # account for convergence failures
        if (TimeSteps[istep].conv_fails > 0):
            CfailH.append(TimeSteps[istep].h_attempts[0]);
            CfailN.append(TimeSteps[istep].NonlinIters);

    # convert data to numpy arrays
    N = np.array(Nvals);
    h = np.array(hvals);
    cN = np.array(CfailN);
    ch = np.array(CfailH);

    # generate plot
    plt.figure()
    plt.semilogx(h,N,'b.')
    plt.semilogx(ch,cN,'rx')
    plt.xlabel('step size')
    plt.ylabel('Nonlinear iters')
    plt.title('Nonlinear iterations per step versus step size')
    plt.legend(('successful','conv. fails'), loc='lower right', shadow=True)
    plt.grid()
    plt.savefig(fname)


##########
def plot_krylov_vs_t(TimeSteps,fname):
    """ This routine takes in the array of TimeSteps (returned from  """
    """ load_diags), and plots the total number of Krylov            """
    """ iterations as a function of the simulation time t.           """
    """                                                              """
    """ The resulting plot is stored in the file <fname>, that       """
    """ should include an extension appropriate for the matplotlib   """
    """ 'savefig' command.                                           """
    import pylab as plt
    import numpy as np
    Kvals   = [];
    tvals   = [];
    for istep in range(len(TimeSteps)):
        
        # store nonlinear iterations and time
        Kvals.append(TimeSteps[istep].KrylovIters);
        tvals.append(TimeSteps[istep].tn);

    # convert data to numpy arrays
    K = np.array(Kvals);
    t = np.array(tvals);

    # generate plot
    plt.figure()
    plt.plot(t,K,'b-')
    plt.xlabel('time')
    plt.ylabel('Krylov iters')
    plt.title('Krylov iterations per step versus time')
    plt.grid()
    plt.savefig(fname)


##########
def plot_krylov_vs_h(TimeSteps,fname):
    """ This routine takes in the array of TimeSteps (returned from  """
    """ load_diags), and plots the total number of Krylov            """
    """ iterations as a function of the time step size h.            """
    """                                                              """
    """ The resulting plot is stored in the file <fname>, that       """
    """ should include an extension appropriate for the matplotlib   """
    """ 'savefig' command.                                           """
    import pylab as plt
    import numpy as np
    Kvals   = [];
    hvals   = [];
    for istep in range(len(TimeSteps)):
        
        # store nonlinear iterations and time
        Kvals.append(TimeSteps[istep].KrylovIters);
        hvals.append(TimeSteps[istep].h_final);

    # convert data to numpy arrays
    K = np.array(Kvals);
    h = np.array(hvals);

    # generate plot
    plt.figure()
    plt.semilogx(h,K,'b.')
    plt.xlabel('step size')
    plt.ylabel('Krylov iters')
    plt.title('Krylov iterations per step versus step size')
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
    nonlinits = 0;
    krylovits = 0;
    stages  = 0;
    nonlinfails = 0;
    lsetups = 0;
    dcon_final = [];
    nsteps = len(TimeSteps)
    for istep in range(nsteps):
        nonlinits += TimeSteps[istep].NonlinIters;
        krylovits += TimeSteps[istep].KrylovIters;
        stages += len(TimeSteps[istep].StageSteps);
        nonlinfails += TimeSteps[istep].conv_fails;
        lsetups += TimeSteps[istep].lsetups;
        for istage in range(len(TimeSteps[istep].StageSteps)):
            dcon_final.append(TimeSteps[istep].StageSteps[istage].NonlinDcon);

    # generate means and percentages
    dcon = np.array(dcon_final);
    if (lsetups == 0):
        nonlin_per_lsetup = 0;
        krylov_per_lsetup = 0;
    else:
        nonlin_per_lsetup = (1.0*nonlinits/lsetups)
        krylov_per_lsetup = (1.0*krylovits/lsetups)

    fptr.write("\nSimulation stats:\n")
    fptr.write("   steps = %i\n" % (nsteps))
    fptr.write("   stages = %i\n" % (stages))
    fptr.write("   nonlinear iterations = %i\n" % (nonlinits))
    fptr.write("   Krylov iterations = %i\n" % (krylovits))
    fptr.write("   LSetup calls = %i\n" % (lsetups))

    fptr.write("\nNonlinear stats:\n")
    fptr.write("   Avg Nonlinear / step  = %g\n" % (1.0*nonlinits/nsteps))
    fptr.write("   Avg Nonlinear / stage = %g\n" % (1.0*nonlinits/stages))
    fptr.write("   Total Nonlinear failures = %i\n" % (nonlinfails))
    fptr.write("   Average Nonlinear / lsetup = %g\n" % (nonlin_per_lsetup))

    if (krylovits > 0):
        fptr.write("\nKrylov stats:\n")
        fptr.write("   Avg Krylov / step  = %g\n" % (1.0*krylovits/nsteps))
        fptr.write("   Avg Krylov / stage = %g\n" % (1.0*krylovits/stages))
        fptr.write("   Average Krylov / lsetup = %g\n" % (krylov_per_lsetup))

    fptr.write("\nNonlinear convergence tests (dcon):\n")
    fptr.write("     min = %g\n" % (np.min(dcon)))
    fptr.write("     max = %g\n" % (np.max(dcon)))
    fptr.write("     avg = %g\n" % (np.mean(dcon)))
    fptr.write("   gmean = %g\n" % (np.prod(np.power(dcon,1.0/len(dcon)))))
    fptr.write("\n")
    
   


##########

