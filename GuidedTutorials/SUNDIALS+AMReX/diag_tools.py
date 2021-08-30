#------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ SMU
#------------------------------------------------------------
# Copyright (c) 2019, Southern Methodist University.
# All rights reserved.
# For details, see the LICENSE file.
#------------------------------------------------------------
# matplotlib-based plotting and analysis utilities for ARKode
# solver diagnostics


#### Data Structures ####

class ErrorTest:
    """
    An ErrorTest object stores, for each error test performed:
        index    -- the time step index
        h        -- the time step size
        estimate -- the estimate of the local error
        errfail  -- a flag denoting whether the error test failed
    """
    def __init__(self, step, h, dsm):
        self.index    = step
        self.h        = h
        if (dsm <= 0.0):
            self.estimate = 1.e-10
        else:
            self.estimate = dsm
        self.errfail  = 0
        if (dsm > 1.0):
            self.errfail = 1
    def Write(self):
        print('  ErrorTest: index =',self.index,', h =',self.h,', estimate =',self.estimate,', errfail =',self.errfail)

##########
class AdaptH:
    """
    An AdaptH object stores, for each time step adaptivity calculation:
        eh[0,1,2]        -- the biased error history array
        hh[0,1,2]        -- the time step history array
        h_accuracy[0,1]  -- the accuracy step estimates, before
                            and after applying step size bounds
        h_stability[0,1] -- the stability step estimates, before
                            and after applying cfl & stability bounds
        stab_restrict    -- flag whether step was stability-limited
        eta              -- the final time step growth factor
    """
    def __init__(self, eh0, eh1, eh2, hh0, hh1, hh2, ha0, hs0, ha1, hs1, eta):
        self.eh0 = eh0
        self.eh1 = eh1
        self.eh2 = eh2
        self.hh0 = hh0
        self.hh1 = hh1
        self.hh2 = hh2
        self.h_accuracy0 = ha0
        self.h_accuracy1 = ha1
        self.h_stability0 = hs0
        self.h_stability1 = hs1
        self.eta = eta
        if (hs0 < ha0):
            self.stab_restrict = 1
        else:
            self.stab_restrict = 0
    def Write(self):
        print('  AdaptH: errhist =',self.eh0,self.eh1,self.eh2,', stephist =',self.hh0,self.hh1,self.hh2,', ha0 =',self.h_accuracy0,', hs0 =',self.h_stability0,', ha1 =',self.h_accuracy1,', hs1 =',self.h_stability1,', stabrestrict =',self.stabrestrict,', eta =',self.eta)

##########
class StageStep:
    """
    A StageStep object stores, for each RK stage of every time step:
        step         -- the time step index
        h            -- the time step size
        stage        -- the stage index
        tn           -- the stage time
    """
    def __init__(self, step, h, stage, tn):
        self.step  = step
        self.h     = h
        self.stage = stage
        self.tn    = tn
    def Write(self):
        print('  StageStep: step =',self.step,', stage =',self.stage)
        print('    h =',self.h,', tn =',self.tn)

##########
class TimeStep:
    """
    A TimeStep object stores, for every time step:
        StageSteps  -- array of StageStep objects comprising the step
        h_attempts  -- array of step sizes attempted (typically only one,
                       unless convergence or error failures occur)
        step        -- the time step index
        tn          -- maximum stage time in step
        h_final     -- the final successful time step size
        ErrTest     -- an ErrorTest object for the step
        err_fails   -- total error test failures in step
    """
    def __init__(self):
        self.StageSteps  = []
        self.h_attempts  = []
        self.step        = -1
        self.tn          = -1.0
        self.h_final     = -1.0
        self.err_fails   =  0
    def AddStage(self, Stage):
        self.step = Stage.step
        if (self.h_final != Stage.h):
            self.h_final = Stage.h
            self.h_attempts.append(Stage.h)
        self.StageSteps.append(Stage)
        self.tn = max(self.tn,Stage.tn)
    def AddErrorTest(self, ETest):
        self.ErrTest = ETest
        self.err_fails += ETest.errfail
    def AddHAdapt(self, HAdapt):
        self.HAdapt = HAdapt
    def Write(self):
        print('TimeStep: step =',self.step,', tn =',self.tn)
        print('  h_attempts =',self.h_attempts)
        print('  h_final =',self.h_final)
        print('  err_fails =',self.err_fails)
        for i in range(len(self.StageSteps)):
            self.StageSteps[i].Write()


#### Utility functions ####

def load_line(line):
    """
    This routine parses a line of the diagnostics output file to
    determine what type of data it contains (an RK stage step,
    an error test, or a time step adaptivity calculation), and
    creates the relevant object for that data line.  Each of
    these output types are indexed by a specific linetype for
    use by the calling routine.

    The function returns [linetype, entry].
    """
    import shlex
    txt = shlex.split(line)
    if ("step" in txt):
        linetype = 0
        step  = int(txt[2])
        h     = float(txt[3])
        stage = int(txt[4])
        tn    = float(txt[5])
        entry = StageStep(step, h, stage, tn)
    elif ("etest" in txt):
        linetype = 3
        step  = int(txt[2])
        h     = float(txt[3])
        dsm   = float(txt[4])
        entry = ErrorTest(step, h, dsm)
    elif ("adapt" in txt):
        linetype = 4
        eh0 = float(txt[2])
        eh1 = float(txt[3])
        eh2 = float(txt[4])
        hh0 = float(txt[5])
        hh1 = float(txt[6])
        hh2 = float(txt[7])
        ha0 = float(txt[8])
        hs0 = float(txt[9])
        ha1 = float(txt[10])
        hs1 = float(txt[11])
        eta = float(txt[12])
        entry = AdaptH(eh0, eh1, eh2, hh0, hh1, hh2, ha0, hs0, ha1, hs1, eta)
    else:
        linetype = -1
        entry = 0
    return [linetype, entry]

##########
def load_diags(fname):
    """
    This routine opens the diagnostics output file, loads all lines
    to create an array of TimeSteps with all of the relevant data
    contained therein.
    """
    f = open(fname)
    step  = -1
    stage = -1
    TimeSteps = []
    for line in f:
        linetype, entry = load_line(line)
        if (linetype == 0):   # stage step
            if (entry.step != step):   # new step
                step = entry.step
                TimeSteps.append(TimeStep())
            TimeSteps[step].AddStage(entry)
            stage = entry.stage
        elif (linetype == 3):   # error test
            TimeSteps[step].AddErrorTest(entry)
        elif (linetype == 4):   # h adaptivity
            TimeSteps[step].AddHAdapt(entry)
    f.close()
    return TimeSteps

##########
def write_diags(TimeSteps):
    """
    This routine takes in the array of TimeSteps (returned from
    load_diags), and writes out the internal representation of
    the time step history to stdout.
    """
    for i in range(len(TimeSteps)):
        print('  ')
        TimeSteps[i].Write()

##########
def plot_h_vs_t(TimeSteps,fname):
    """
    This routine takes in the array of TimeSteps (returned from
    load_diags), and plots the time step sizes h as a function
    of the simulation time t.  Failed time steps are marked on
    the plot with either a red X or green O, where X corresponds
    to an error test failure.

    The resulting plot is stored in the file <fname>, that
    should include an extension appropriate for the matplotlib
    'savefig' command.
    """
    import pylab as plt
    import numpy as np
    hvals  = []
    tvals  = []
    EfailH = []
    EfailT = []
    CfailH = []
    CfailT = []
    for istep in range(len(TimeSteps)):

        # store successful step size and time
        hvals.append(TimeSteps[istep].h_final)
        tvals.append(TimeSteps[istep].tn)

        # account for convergence failures and error test failures
        if (TimeSteps[istep].err_fails > 0):
            EfailH.append(TimeSteps[istep].h_attempts[0])
            EfailT.append(TimeSteps[istep].tn)

    # convert data to numpy arrays
    h = np.array(hvals)
    t = np.array(tvals)
    eh = np.array(EfailH)
    et = np.array(EfailT)

    # generate plot
    plt.figure()
    plt.semilogy(t,h,'b-')
    if (len(eh) > 0):
        plt.semilogy(et,eh,'rx')
    plt.xlabel('time')
    plt.ylabel('step size')
    plt.title('Step size versus time')
    if (len(eh) > 0):
        plt.legend(('successful','error fails'), loc='lower right', shadow=True)
    plt.grid()
    plt.savefig(fname)


##########
def plot_h_vs_tstep(TimeSteps,fname):
    """
    This routine takes in the array of TimeSteps (returned from
    load_diags), and plots the time step sizes h as a function
    of the time step iteration index.  Failed time steps are
    marked on the plot a red X if an error test failure occurred.

    The resulting plot is stored in the file <fname>, that
    should include an extension appropriate for the matplotlib
    'savefig' command.
    """
    import pylab as plt
    import numpy as np
    hvals  = []
    ivals  = []
    EfailH = []
    EfailI = []
    for istep in range(len(TimeSteps)):

        # store successful step size and time
        hvals.append(TimeSteps[istep].h_final)
        ivals.append(istep)

        # account for convergence failures and error test failures
        if (TimeSteps[istep].err_fails > 0):
            EfailH.append(TimeSteps[istep].h_attempts[0])
            EfailI.append(istep)

    # convert data to numpy arrays
    h = np.array(hvals)
    I = np.array(ivals)
    eh = np.array(EfailH)
    eI = np.array(EfailI)

    # generate plot
    plt.figure()
    plt.semilogy(I,h,'b-')
    if (len(eI) > 0):
        plt.semilogy(eI,eh,'rx')
    plt.xlabel('time step')
    plt.ylabel('step size')
    plt.title('Step size versus time step')
    if (len(eI) > 0):
        plt.legend(('successful','error fails'), loc='lower right', shadow=True)
    plt.grid()
    plt.savefig(fname)


##########
def plot_oversolve_vs_t(TimeSteps,fname):
    """
    This routine takes in the array of TimeSteps (returned from
    load_diags), and plots the oversolve as a function of the
    simulation time t. We cap the computed oversolve value at
    1000 to get more intuitive plots since first few time steps
    are purposefully too small.

    The resulting plot is stored in the file <fname>, that
    should include an extension appropriate for the matplotlib
    'savefig' command.
    """
    import pylab as plt
    import numpy as np
    Ovals = []
    tvals = []
    for istep in range(len(TimeSteps)):

        # store oversolve and time
        Ovals.append(min(1e3,1.0 / TimeSteps[istep].ErrTest.estimate))
        tvals.append(TimeSteps[istep].tn)

    # convert data to numpy arrays
    O = np.array(Ovals)
    t = np.array(tvals)

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
    """
    This routine takes in the array of TimeSteps (returned from
    load_diags), and computes statistics on how well the time
    step adaptivity estimations predicted step values that met
    the accuracy requirements.

    Note: we ignore steps immediately following an
    error test failure, since etamax is bounded above by 1.

    The resulting data is appended to the stream corresponding
    to fptr (either the result of 'open' or sys.stdout).
    """
    import numpy as np
    errfails     = 0
    oversolve10  = 0
    oversolve100 = 0
    oversolve_max = 0.0
    oversolve_min = 1.0e200
    oversolves = []
    nsteps = len(TimeSteps)
    hvals = []
    for istep in range(nsteps):
        if (TimeSteps[istep].err_fails > 0):
            errfails += 1
        hvals.append(TimeSteps[istep].h_final)

        # if this or previous step had an error test failure, skip oversolve results
        ignore = 0
        if (istep > 0):
            if( (TimeSteps[istep].err_fails > 0) or (TimeSteps[istep-1].err_fails > 0)):
                ignore = 1
        else:
            if(TimeSteps[istep].err_fails > 0):
                ignore = 1

        if (ignore == 0):
            over = 1.0 / TimeSteps[istep].ErrTest.estimate
            oversolves.append(over)
            if (over > 100.0):
                oversolve100 += 1
            elif (over > 10.0):
                oversolve10 += 1
            oversolve_max = max(oversolve_max, over)
            oversolve_min = min(oversolve_min, over)

    # generate means and percentages
    ov = np.array(oversolves)
    hv = np.array(hvals)
    hval_mean = np.mean(hv)
    hval_max = np.max(hv)
    hval_min = np.min(hv)
    oversolve_mean = np.mean(ov)
    oversolve10  = 100.0*oversolve10/nsteps
    oversolve100 = 100.0*oversolve100/nsteps

    fptr.write("\n")
    fptr.write("Simulation took %i steps\n" % (nsteps))
    fptr.write("  Stepsize statistics:\n")
    fptr.write("           min = %.3g\n" % (hval_min))
    fptr.write("           max = %.3g\n" % (hval_max))
    fptr.write("          mean = %.3g\n" % (hval_mean))
    fptr.write("  Failed steps = %i (%.2g%%)\n" % (errfails, 100.0*errfails/nsteps))
    fptr.write("  Oversolve fractions:\n")
    fptr.write("           10x = %g%%\n" % (oversolve10))
    fptr.write("          100x = %g%%\n" % (oversolve100))
    fptr.write("  Oversolve statistics:\n")
    fptr.write("           min = %.3g\n" % (oversolve_min))
    fptr.write("           max = %.3g\n" % (oversolve_max))
    fptr.write("          mean = %.3g\n" % (oversolve_mean))
    fptr.write("\n")


##########
