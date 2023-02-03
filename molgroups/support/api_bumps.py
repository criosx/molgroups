from __future__ import print_function
from os import path
from random import seed, random
from re import VERBOSE, IGNORECASE, compile
from sys import stdout

import matplotlib
from matplotlib import pyplot as plt

from bumps.cli import load_model, save_best
from bumps.mapper import MPMapper
from bumps.fitters import fit, FitDriver, DreamFit, LevenbergMarquardtFit, MPFit

import numpy
import shutil
import glob
import os

from molgroups.support import api_base


class CBumpsAPI(api_base.CBaseAPI):
    def __init__(self, spath=".", mcmcpath=".", runfile="", state=None, problem=None, load_state=True):
        super().__init__(spath, mcmcpath, runfile)
        if load_state:
            self.state = self.fnRestoreState() if state is None else state
        self.problem = self.fnRestoreFitProblem() if problem is None else problem

    def fnBackup(self, origin=None, target=None):
        if origin is None:
            origin = self.spath
        if target is None:
            target = self.spath + '/rsbackup'
        if not path.isdir(target):
            os.mkdir(target)
        for file in glob.glob(origin + r'/*.dat'):
            shutil.copy(file, target)
        for file in glob.glob(origin + r'/*.py'):
            shutil.copy(file, target)
        for file in glob.glob(origin + r'/*.pyc'):
            shutil.copy(file, target)

    def fnLoadMCMCResults(self):
        # load Parameter
        if self.diParameters == {}:
            self.diParameters, _ = self.fnLoadParameters()
        lParName = self.diParameters.keys()

        # parvars is a list of variables(parameters) to return for each point
        parvars = [i for i in range(len(lParName))]
        draw = self.state.draw(1, parvars, None)
        points = draw.points
        logp = draw.logp

        return points, lParName, logp

    # LoadStatResults returns a list of variable names, a logP array, and a numpy.ndarray
    # [values,var_numbers].
    def fnLoadParameters(self):
        if self.state is None:
            self.state = self.fnRestoreState()

        p = self.state.best()[0]
        self.problem.setp(p)

        # from bumps.cli import load_best
        # load_best(problem, os.path.join(self.mcmcpath, self.runfile) + '.par')

        # distinguish between fitproblem and multifitproblem
        if "models" in dir(self.problem):
            i = 0
            overall = 0.0
            for M in self.problem.models:
                overall += M.chisq()
                i += 1
            overall /= float(i)
            pnamekeys = []
            pars = self.problem._parameters
            for par in pars:
                pnamekeys.append(par.name)
        else:
            overall = self.problem.chisq()
            pnamekeys = self.problem.labels()

        # Do not accept parameter names with spaces, replace with underscore
        for i in range(len(pnamekeys)):
            pnamekeys[i] = pnamekeys[i].replace(" ", "_")

        for element in pnamekeys:
            self.diParameters[element] = {}
        bounds = self.problem.bounds()

        for key in pnamekeys:
            parindex = pnamekeys.index(key)
            self.diParameters[key]["number"] = parindex
            self.diParameters[key]["lowerlimit"] = float(bounds[0][parindex])
            self.diParameters[key]["upperlimit"] = float(bounds[1][parindex])
            self.diParameters[key]["value"] = float(p[parindex])
            self.diParameters[key]["relval"] = (
                self.diParameters[key]["value"] - self.diParameters[key]["lowerlimit"]
            ) / (
                self.diParameters[key]["upperlimit"]
                - self.diParameters[key]["lowerlimit"]
            )
            # TODO: Do we still need this? Set to dummy value
            self.diParameters[key]["variable"] = key
            # TODO: Do we still need this? Would have to figure out how to get the confidence limits from state
            self.diParameters[key]["error"] = 0.01

        return self.diParameters, overall

    def fnLoadStatData(self, dSparse=0, rescale_small_numbers=True, skip_entries=None):
        if skip_entries is None:
            skip_entries = []

        if path.isfile(os.path.join(self.spath, self.mcmcpath, "sErr.dat")) or \
                path.isfile(os.path.join(self.spath, self.mcmcpath, "isErr.dat")):
            diStatRawData = self.fnLoadsErr()
        else:
            points, lParName, logp = self.fnLoadMCMCResults()
            diStatRawData = {"Parameters": {}}
            diStatRawData["Parameters"]["Chisq"] = {}
            # TODO: Work on better chisq handling
            diStatRawData["Parameters"]["Chisq"]["Values"] = []
            for parname in lParName:
                diStatRawData["Parameters"][parname] = {}
                diStatRawData["Parameters"][parname]["Values"] = []

            seed()
            for j in range(len(points[:, 0])):
                if dSparse == 0 or (dSparse > 1 and j < dSparse) or (1 > dSparse > random()):
                    diStatRawData["Parameters"]["Chisq"]["Values"].append(logp[j])
                    for i, parname in enumerate(lParName):
                        # TODO: this is a hack because Paul does not scale down after scaling up
                        # Rescaling disabled for bumps/refl1d analysis to achieve consistency
                        # if ('rho_' in parname or 'background' in parname) and rescale_small_numbers:
                        #     points[j, i] *= 1E-6
                        diStatRawData["Parameters"][parname]["Values"].append(points[j, i])

            self.fnSaveSingleColumnsFromStatDict(os.path.join(self.spath, self.mcmcpath, "sErr.dat"),
                                                 diStatRawData["Parameters"], skip_entries)

        return diStatRawData

    # deletes the backup directory after restoring it
    def fnRemoveBackup(self, origin=None, target=None):
        if origin is None:
            origin = self.spath
        if target is None:
            target = self.spath + '/rsbackup'
        if path.isdir(target):
            self.fnRestoreBackup(origin, target)
            shutil.rmtree(target)

    def fnReplaceParameterLimitsInSetup(self, sname, flowerlimit, fupperlimit):
        """
        Scans self.runfile file for parameter with name sname and replaces the
        lower and upper fit limits by the given values.
        If a initialization value is given as part of the Parameter() function, it will be replaced as well.
        The value= argument should follow the name= argument in Parameter()
        """

        file = open(os.path.join(self.spath, self.runfile) + '.py', 'r+')
        data = file.readlines()
        file.close()
        smatch = compile(r"(.*?Parameter.*?name=\'"+sname+"[\"\'].+?=).+?(\).+?range\().+?(,).+?(\).*)",
                         IGNORECASE | VERBOSE)
        # version when .range() is present but no parameter value is provided
        smatch2 = compile(r"(.*?\."+sname+'\.range\().+?(,).+?(\).*)', IGNORECASE | VERBOSE)
        newdata = []
        for line in data:
            # apply version 1 for general case
            newline = smatch.sub(r'\1 ' + str(0.5*(flowerlimit+fupperlimit)) + r'\2 ' + str(flowerlimit) + r'\3 ' +
                                 str(fupperlimit) + r'\4', line)
            # apply version 2 to catch both cases, potentially redundant for limits
            newline = smatch2.sub(r'\1 ' + str(flowerlimit) + r'\2 ' + str(fupperlimit) + r'\3', newline)
            newdata.append(newline)

        file = open(os.path.join(self.spath, self.runfile) + '.py', 'w')
        file.writelines(newdata)
        file.close()

    # copies all files from the backup directory (target) to origin
    def fnRestoreBackup(self, origin=None, target=None):
        if origin is None:
            origin = self.spath
        if target is None:
            target = self.spath + '/rsbackup'
        if path.isdir(target):
            for file in glob.glob(target + r'/*.*'):
                shutil.copy(file, origin)

    def fnRestoreFit(self):
        self.problem = self.fnRestoreFitProblem()
        self.state = self.fnRestoreState()
        if self.state is not None:
            # repopulate state with best fit
            p = self.state.best()[0]
            self.problem.setp(p)

    def fnRestoreFitProblem(self):
        from bumps.fitproblem import load_problem

        if path.isfile(os.path.join(self.spath, self.runfile + ".py")):
            problem = load_problem(os.path.join(self.spath, self.runfile + ".py"))
        else:
            print("No problem to reload.")
            problem = None
        return problem

    def fnRestoreState(self):
        import bumps.dream.state
        fulldir = os.path.join(self.spath, self.mcmcpath)
        if path.isfile(os.path.join(fulldir, self.runfile) + '.py') and path.isfile(os.path.join(fulldir, self.runfile)
                                                                                    + '-chain.mc.gz'):
            state = bumps.dream.state.load_state(os.path.join(fulldir, self.runfile))
            state.mark_outliers()  # ignore outlier chains
        else:
            print("No file: " + os.path.join(fulldir, self.runfile) + '.py')
            print("No state to reload.")
            state = None
        return state

    def fnRestoreMolgroups(self, problem):
        # Populates the diMolgroups dictionary based from a saved molgroups.dat file
        diMolgroups = self.fnLoadMolgroups(problem=problem)
        return diMolgroups

    def fnRestoreSmoothProfile(self, M):
        # TODO: Decide what and if to return SLD profile for Bumps fits
        # Returns currently profile for zeroth model if multiproblem fit
        z, rho, irho = M.sld, [], []
        return z, rho, irho

    def fnRunMCMC(self, burn=8000, steps=500, batch=False, fitter='MCMC', reload_problem=True):
        """
        Runs fit for Bumps object.
        Default is 'MCMC', but 'LM' is supported, as well.
        'reload_problem' determines whether the problem is reloaded from disk or whether the internally stored problem
        is used, including any potential best-fit parameters from a previous run or restore.
        """

        # Original Method of Calling the Shell
        '''
        lCommand = ['refl1d', os.path.join(self.spath, self.runfile)+'.py', '--fit=dream', '--parallel', '--init=lhs']
        if batch:
           lCommand.append('--batch')
        lCommand.append('--store=' + self.mcmcpath)
        lCommand.append('--burn=' + str(burn))
        lCommand.append('--steps=' + str(steps))
        lCommand.append('--overwrite')
        call(lCommand)
        '''

        # Command line Python implementation
        '''
        from refl1d import main
        # There is a bug in bumps that prevents running sequential fits because the pool object is terminated
        # from a previous fit but not None. Bumps expects it to be None or will not initiate a new pool.        
        if none_pool:
            from bumps.mapper import MPMapper
            # MPMapper.pool.terminate()  # not always required
            MPMapper.pool = None

        sys.argv = ['refl1d', os.path.join(self.spath, self.runfile)+'.py', '--fit=dream', '--parallel', '--init=lhs']
        if batch:
            sys.argv.append('--batch')
        sys.argv.append('--store=' + self.mcmcpath)
        sys.argv.append('--burn=' + str(burn))
        sys.argv.append('--steps=' + str(steps))
        sys.argv.append('--overwrite')

        main.cli()
        '''

        # Calling bumps functions directly

        model_file = os.path.join(self.spath, self.runfile) + '.py'
        mcmcpath = os.path.join(self.spath, self.mcmcpath)

        if not os.path.isdir(mcmcpath):
            os.mkdir(mcmcpath)

        # save model file in output directory
        shutil.copy(model_file, mcmcpath)

        if reload_problem or self.problem is None:
            self.problem = self.fnRestoreFitProblem()

        mapper = MPMapper.start_mapper(self.problem, None, cpus=0)
        monitors = None if not batch else []

        if fitter == 'MCMC':
            driver = FitDriver(fitclass=DreamFit, mapper=mapper, problem=self.problem, init='lhs', steps=steps,
                               burn=burn, monitors=monitors, xtol=1e-6, ftol=1e-8)
        elif fitter == 'LM':
            driver = FitDriver(fitclass=MPFit, mapper=mapper, problem=self.problem, monitors=monitors,
                               steps=steps, xtol=1e-6, ftol=1e-8)
        else:
            driver = None

        x, fx = driver.fit()

        # try to deal with matplotlib memory leaks
        matplotlib.interactive(False)

        # .err and .par files
        self.problem.output_path = os.path.join(mcmcpath, self.runfile)
        save_best(driver, self.problem, x)

        # try to deal with matplotlib cache issues by deleting the cache
        fig = plt.figure()
        plt.figure().clear()
        plt.cla()
        plt.clf()
        plt.close("all")

        # don't know what files
        if 'models' in dir(self.problem):
            for M in self.problem.models:
                M.fitness.save(os.path.join(mcmcpath, self.runfile))
                break
        else:
            self.problem.fitness.save(os.path.join(mcmcpath, self.runfile))

        # .mcmc and .point files
        driver.save(os.path.join(mcmcpath, self.runfile))

        if not batch:
            # stat table and yet other files
            driver.show()
            # plots and other files
            driver.plot(os.path.join(mcmcpath, self.runfile))

    def fnSaveMolgroups(self, problem):
        # saves bilayer and protein information from a bumps / refl1d problem object into a mol.dat file
        # sequentially using the methods provided in molgroups
        fp = open(self.spath + '/mol.dat', "w")
        z = numpy.linspace(0, problem.dimension * problem.stepsize, problem.dimension, endpoint=False)
        try:
            problem.extra[0].fnWriteGroup2File(fp, 'bilayer', z)
            problem.extra[1].fnWriteGroup2File(fp, 'protein', z)
        except:
            problem.extra.fnWriteGroup2File(fp, 'bilayer', z)
        fp.close()
        stdout.flush()

    def fnUpdateModelPars(self, diNewPars):
        liParameters = list(self.diParameters.keys())
        # sort by number of appereance in runfile
        liParameters = sorted(liParameters, key=lambda keyitem: self.diParameters[keyitem]['number'])
        for element in liParameters:
            if element not in list(diNewPars.keys()):
                print('Parameter ' + element + ' not specified.')
                # check failed -> exit method
                return
            # else:
                # print(element + ' ' + str(diNewPars[element]))

        p = [diNewPars[parameter] for parameter in liParameters]
        self.problem.setp(p)
        self.problem.model_update()


