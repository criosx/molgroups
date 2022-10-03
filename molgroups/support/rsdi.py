from __future__ import print_function
from math import fabs, pow, sqrt
from os import path
from random import seed, normalvariate, random
from re import VERBOSE, IGNORECASE, compile
from sys import exit, stdout
from subprocess import call, Popen
from time import sleep
import numpy
import pandas
import shutil
import glob
import os

from molgroups.support import general

class CDataInteractor:
    def __init__(self, spath='.', mcmcpath='.', runfile=''):
        # script path, MCMC path, runfile (script file)
        self.spath = spath
        self.mcmcpath = mcmcpath
        self.runfile = runfile
        self.diParameters = {}

    def fnBackup(self):
        raise NotImplementedError()

    def fnBackupSimdat(self):
        raise NotImplementedError()

    def fnLoadMolgroups(self, problem=None):
        diMolgroups = {}
        moldict = problem.moldat

        for group in moldict:
            tdata = (moldict[group]['header']).split()  # read header that contains molgroup data
            diMolgroups[tdata[1]] = {}
            diMolgroups[tdata[1]].update({'headerdata': {}})
            diMolgroups[tdata[1]]['headerdata'].update({'Type': tdata[0]})
            diMolgroups[tdata[1]]['headerdata'].update({'ID': tdata[1]})
            for j in range(2, len(tdata), 2):
                diMolgroups[tdata[1]]['headerdata'].update({tdata[j]: tdata[j + 1]})

            zax = moldict[group]['zaxis']
            areaax = moldict[group]['area']
            nslax = moldict[group]['nsl']
            diMolgroups[tdata[1]].update({'zaxis': zax, 'areaaxis': areaax, 'nslaxis': nslax})

        return diMolgroups

    def fnLoadsErr(self):
        sFileName = self.mcmcpath + '/isErr.dat'
        try:
            if path.isfile(self.mcmcpath+'/isErr.dat') and path.isfile(self.mcmcpath+'/sErr.dat'):
                print('-------------------------------')
                print('Found isErr.dat and sErr.dat ?!')
                print('Load isErr.dat as default.')
                print('-------------------------------')
            elif path.isfile(self.mcmcpath+'/isErr.dat'):  # check which type of MC output present
                print('Found isErr.dat\n')
            elif path.isfile(self.mcmcpath+'/sErr.dat'):
                print('Found sErr.dat\n')
                sFileName = self.mcmcpath+'/sErr.dat'

            diStatRawData = {'Parameters': self.fnLoadSingleColumnsIntoStatDict(sFileName)}
            return diStatRawData

        except IOError:
            print('Could not load ' + sFileName + '. \n')
            exit(1)

    # The sparse parameter loads only a fraction of the data, if sparse is larger or equal than 1 this equates to the
    # number of lines loaded. If sparse is smaller than one this equates to a probability that any line will be stored
    @staticmethod
    def fnLoadSingleColumns(sFileName, data=None, exceptions=None, header=True, headerline=None, LoadList=None,
                            sparse=0):

        File = open(sFileName, "r")
        content = File.readlines()
        File.close()

        if data is None:
            data = {}

        if header:  # if headerline in file, read it from there
            splitheaderline = content[0].split()
            content = content[1:]
        else:  # else use the one given as an attribute
            splitheaderline = headerline

        for i, columnname in enumerate(splitheaderline):
            if LoadList is None or (columnname in LoadList):
                data[columnname] = []

        j = 0
        seed()
        for line in content:
            if sparse == 0 or (sparse >= 1 and j < sparse) or (1 > sparse > random()):
                splitline = line.split()
                if splitline:
                    for i, column in enumerate(splitline):
                        if LoadList is None or (splitheaderline[i] in LoadList):
                            try:
                                data[splitheaderline[i]].append(float(splitline[i]))
                            except:
                                data[splitheaderline[i]].append(splitline[i])
                j += 1

        return data

    # Loads single row data into a dictionary ["rowname",'Values',[data]]
    # it appends the data found in the file to the provided list 'data'
    # it will skip rows with names which are either in the exception list
    # or appends data with name extension "_2nd" that are already present in the data list
    # The sparse parameter loads only a fraction of the data, if sparse is larger or equal than 1 this equates to the
    # number of lines loaded. If sparse is smaller than one this equates to a probability that any line will be stored
    @staticmethod
    def fnLoadSingleColumnsIntoStatDict(sFileName, data=None, exceptions=None, header=True, headerline=None,
                                        LoadList=None, sparse=0):

        if not data: data = {}
        file = open(sFileName, "r")
        content = file.readlines()
        file.close()

        if data is None:
            data = {}

        if header:  # if headerline in file, read it from there
            splitheaderline = content[0].split()
            content = content[1:]
        else:  # else use the one given as an attribute
            splitheaderline = headerline

        for i, columnname in enumerate(splitheaderline):
            if LoadList is None or (columnname in LoadList):
                data[columnname] = {}
                data[columnname]['Values'] = []

        j = 0
        seed()
        for line in content:
            if sparse == 0 or (sparse >= 1 and j < sparse) or (sparse < 1 and random() < sparse):
                splitline = line.split()
                if splitline:
                    for i, column in enumerate(splitline):
                        if LoadList is None or (splitheaderline[i] in LoadList):
                            data[splitheaderline[i]]['Values'].append(float(splitline[i]))
                j += 1

        return data

    @staticmethod
    def fnLoadSingleRows(sFileName, data=None, exceptions=None):

        file = open(sFileName, "r")
        content = file.readlines()
        file.close()

        if data is None:
            data = {}

        for line in content:
            splitline = line.split()
            if splitline[0] not in exceptions:
                data[splitline[0]] = []
                for entry in range(1, len(splitline)):
                    data[splitline[0]].append(splitline[entry])
        return data

    @staticmethod
    def fnMCModifyFile(filelist):
        # reads original reflectivity files
        # and modifies them with normal deviates
        # works also with any other 3 or 4
        # column file types

        seed()  # initialize random number generator

        for reflfile in filelist:  # iterate over all refl files
            file = open(reflfile)
            data = file.readlines()
            file.close()
            newdata = []
            iFileType = 0
            for line in data:
                if '#' in line:  # do not modify comment lines
                    newdata.append(line)
                else:
                    columns = line.split()  # access columns
                    if not iFileType:  # Stuart's Mod for 3 Column i.e. NIST
                        if len(columns) == 3:
                            iFileType = 1
                        elif len(columns) == 4:
                            iFileType = 2  # Stuart's Mod for 4 Column i.e. ISIS
                        else:
                            iFileType = 99  # unrecognized format
                    if not len(columns):
                        print('Found empty line in data file.')
                    if len(columns) == 3 and iFileType == 1:
                        try:
                            fvalue = float(columns[1])  # reflectivity
                            ferror = float(columns[2])  # error
                            columns[1] = str(fvalue + normalvariate(0, 1) * ferror)  # modify reflectivity
                            newline = ''
                            for column in columns:  # glue columns together to one line
                                newline = newline + column + ' '
                            newline = newline[:-1] + '\n'  # add newline to end of line
                            newdata.append(newline)  # add line to new data file
                        except:
                            print('-----------------------------------')
                            print('Data file %s corrupt.' % (reflfile))
                            print('File was identified being NIST type')
                            print('-----------------------------------')
                            raise RuntimeError('Corrupt Data file.')
                    elif len(columns) == 4 and iFileType == 2:
                        try:
                            fvalue = float(columns[2])  # reflectivity
                            ferror = float(columns[3])  # error
                            columns[2] = str(fvalue + normalvariate(0, 1) * ferror)  # modify reflectivity
                            newline = ''
                            for column in columns:  # glue columns together to one line
                                newline = newline + column + ' '
                            newline = newline[:-1] + '\n'  # add newline to end of line
                            newdata.append(newline)  # add line to new data file
                        except:
                            print('-----------------------------------')
                            print('Data file %s corrupt.' % (reflfile))
                            print('File was identified being ISIS type.')
                            print('-----------------------------------')
                            raise RuntimeError('Corrupt Data file.')

                    else:
                        print('-----------------------------------------------------')
                        print('Filetype not recognized or contains errors: %s' % (reflfile))
                        print('-----------------------------------------------------')
                        raise RuntimeError('Corrupt Data file.')

            file = open(path.split(reflfile)[-1] + '.mce', "w")  # write modified file into .mce file in
            file.writelines(newdata)  # working directory
            file.close()

    def fnRemoveBackup(self):  # deletes the backup directory
        raise NotImplementedError()

    # saves all data out to a file
    @staticmethod
    def fnSaveSingleColumns(sFilename, data):

        file = open(sFilename, "w")

        for element in data:
            file.write(element + " ")
        file.write("\n")

        for i in range(len(data[list(data.keys())[0]])):
            for element in data:
                file.write(str(data[element][i]) + " ")
            file.write("\n")

        file.close()
        # saves all data out to a file

    @staticmethod
    def fnSaveSingleColumnsFromStatDict(sFilename, data, skipentries=None):

        if skipentries is None:
            skipentries = []
        File = open(sFilename, "w")

        for element in data:
            if element not in skipentries:
                File.write(element + " ")
        File.write("\n")

        for i in range(len(data[list(data.keys())[0]]['Values'])):
            for element in data:
                if element not in skipentries:
                    File.write(str(data[element]['Values'][i]) + " ")
            File.write("\n")

        File.close()

    def fnSimulateData(self, liExpression):
        raise NotImplementedError()

    def fnWriteConstraint2Runfile(self, liExpression):
        raise NotImplementedError()


# Garefl methods will be used if a storage directory for a Markov Chain Monte Carlo (MCMC)
# error analysis cannot be found.
# The MCMC directory is called 'MCMC'
class CBumpsInteractor(CDataInteractor):
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
            # TODO: Not sure this is for bumps(diffraction) or strictly for single model fitproblems
            pnamekeys = list(self.problem.model_parameters().keys())

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

    # copies all files from the backup directory (target) to origin
    def fnRestoreBackup(self, origin=None, target=None):
        if origin is None:
            origin = self.spath
        if target is None:
            target = self.spath + '/rsbackup'
        if path.isdir(target):
            for file in glob.glob(target + r'/*.*'):
                shutil.copy(file, origin)

    def fnRestoreFitProblem(self):
        from bumps.fitproblem import load_problem

        if path.isfile(self.spath + "/" + self.runfile + ".py"):
            problem = load_problem(self.spath + "/" + self.runfile + ".py")
        else:
            print("No problem to reload.")
            problem = None
        return problem

    def fnRestoreState(self):
        import bumps.dream.state
        fulldir = os.path.join(self.spath, self.mcmcpath)
        if path.isfile(os.path.join(fulldir, self.runfile) + '.py'):
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

    def fnRunMCMC(self, burn, steps, batch=False):
        # Original Method of Calling the Shell
        """
        lCommand = ['refl1d', os.path.join(self.spath, self.runfile)+'.py', '--fit=dream', '--parallel', '--init=lhs']
        if batch:
           lCommand.append('--batch')
        lCommand.append('--store=' + self.mcmcpath)
        lCommand.append('--burn=' + str(burn))
        lCommand.append('--steps=' + str(steps))
        lCommand.append('--overwrite')
        call(lCommand)
        """

        # Command line Python implementation
        """
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
        """

        # Calling refl1d functions directly
        from bumps.cli import load_model, save_best
        from bumps.mapper import MPMapper
        from bumps.fitters import fit, FitDriver, DreamFit

        model_file = os.path.join(self.spath, self.runfile) + '.py'
        mcmcpath = os.path.join(self.spath, self.mcmcpath)

        if not os.path.isdir(mcmcpath):
            os.mkdir(mcmcpath)

        # save model file in output directory
        shutil.copy(model_file, mcmcpath)

        problem = load_model(model_file)
        mapper = MPMapper.start_mapper(problem, None, cpus=0)
        monitors = None if not batch else []
        driver = FitDriver(fitclass=DreamFit, mapper=mapper, problem=problem, init='lhs', steps=steps, burn=burn,
                           monitors=monitors, xtol=1e-6, ftol=1e-8)
        x, fx = driver.fit()

        # .err and .par files
        problem.output_path = os.path.join(mcmcpath, self.runfile)
        save_best(driver, problem, x)

        # don't know what files
        if 'models' in dir(problem):
            for M in problem.models:
                M.fitness.save(os.path.join(mcmcpath, self.runfile))
                break
        else:
            problem.fitness.save(os.path.join(mcmcpath, self.runfile))

        # .mcmc and .point files
        driver.save(os.path.join(mcmcpath, self.runfile))

        # stat table and yet other files
        driver.show()

        # plots and other files
        if not batch:
            driver.plot(os.path.join(mcmcpath, self.runfile))

    def fnSaveMolgroups(self, problem):
        # saves bilayer and protein information from a bumps / refl1d problem object into a mol.dat file
        # sequentially using the methods provided in molgroups
        fp = open(self.spath + '/mol.dat', "w")
        z = numpy.linspace(0, problem.dimension * problem.stepsize, problem.dimension, endpoint=False)
        try:
            problem.extra[0].fnWritePar2File(fp, 'bilayer', z)
            problem.extra[1].fnWritePar2File(fp, 'protein', z)
        except:
            problem.extra.fnWritePar2File(fp, 'bilayer', z)
        fp.close()
        stdout.flush()


# Refl1D methods will be used if a storage directory for a Markov Chain Monte Carlo (MCMC)
# error analysis are found.
# The MCMC directory is called 'MCMC'
# The refl1d script name has to be run.py.
class CRefl1DInteractor(CBumpsInteractor):
    def __init__(self, spath='.', mcmcpath='.', runfile='', load_state=True):
        super().__init__(spath, mcmcpath, runfile, load_state=load_state)

    def fnReplaceParameterLimitsInSetup(self, sname, flowerlimit, fupperlimit):
        """
        Scans self.runfile file for parameter with name sname and replaces the
        lower and upper fit limits by the given values.
        Currently, expects the parameter to be defined using the Parameter() method and not just .range()
        on any object variable.
        """

        file = open(os.path.join(self.spath, self.runfile) + '.py', 'r+')
        data = file.readlines()
        file.close()
        smatch = compile(r"(.*?Parameter.*?name=\'"+sname+".+?=).+?(\).+?range\().+?(,).+?(\).*)", IGNORECASE | VERBOSE)
        newdata = []
        for line in data:
            newdata.append(smatch.sub(r'\1 ' + str(0.5*(flowerlimit+fupperlimit)) + r'\2 ' + str(flowerlimit) + r'\3 '
                                      + str(fupperlimit) + r'\4', line))

        file = open(os.path.join(self.spath, self.runfile) + '.py', 'w')
        file.writelines(newdata)
        file.close()

    @staticmethod
    def fnRestoreSmoothProfile(M):
        z, rho, irho = M.fitness.smooth_profile()
        return z, rho, irho

    def fnSimulateData(self, diNewPars):
        liParameters = list(self.diParameters.keys())
        # sort by number of appereance in runfile
        liParameters = sorted(liParameters, key=lambda keyitem: self.diParameters[keyitem]['number'])
        bConsistency = True
        for element in liParameters:
            if element not in list(diNewPars.keys()):
                bConsistency = False
                print('Parameter '+element+' not specified.')
            else:
                print(element + ' ' + str(diNewPars[element]))

        if bConsistency:
            p = [diNewPars[parameter] for parameter in liParameters]
            self.problem.setp(p)
            self.problem.model_update()
            # TODO: By calling .chisq() I currently force an update of the cost function. There must be a better way

            i = 0
            if 'models' in dir(self.problem):
                for M in self.problem.models:
                    M.chisq()
                    qvec, refl = M.fitness.reflectivity()
                    comments = general.extract_comments_from_file(self.spath + '/sim' + str(i) + '.dat', "#")
                    simdata = pandas.read_csv(self.spath + '/sim' + str(i) + '.dat', sep='\s+', header=None,
                                              names=['Q', 'R', 'dR', 'dQ'], skip_blank_lines=True, comment='#')
                    simdata['R'] = refl
                    simdata['Q'] = qvec
                    simdata.to_csv('sim' + str(i) + '.dat', sep=' ', index=None, header=None)
                    general.add_comments_to_start_of_file(self.spath + '/sim' + str(i) + '.dat', comments)
                    i += 1

            else:
                self.problem.chisq()
                qvec, refl = self.problem.fitness.reflectivity()
                comments = general.extract_comments_from_file(self.spath + '/sim' + str(i) + '.dat', "#")
                simdata = pandas.read_csv(self.spath + '/sim' + str(i) + '.dat', sep='\s+', header=None,
                                          names=['Q', 'R', 'dR', 'dQ'], skip_blank_lines=True, comment='#')
                simdata['R'] = refl
                simdata['Q'] = qvec
                simdata.to_csv('sim' + str(i) + '.dat', sep=' ', index=None, header=None)
                general.add_comments_to_start_of_file(self.spath + '/sim' + str(i) + '.dat', comments)
                i += 1

    def fnSimulateErrorBars(self, simpar, qmin=0.008, qmax=0.325, s1min=0.108, s1max=4.397, s2min=0.108, s2max=4.397,
                            tmin=18, tmax=208, nmin=11809, rhomin=-0.56e-6, rhomax=6.34e-6, cbmatmin=1.1e-5,
                            cbmatmax=1.25e-6, mode='water', pre=1):
        i = 0
        last_defined_rho_solv = 0
        while path.isfile(self.spath+'/sim' + str(i) + '.dat'):
            # add error bars
            c1 = s1max / qmax
            c2 = s2max / qmax
            c4 = (tmax - tmin) / (qmax ** 2 - qmin ** 2)
            c3 = tmax - c4 * qmax ** 2
            I = nmin / s1min / s2min / tmin * pow(2, pre)

            if mode == 'air':
                # currently hardwire background for air
                cbmat = 1e-7
            else:
                # dataframe can be empty in case that rho_solv_i is not specified
                # for example in magnetic fits, in this case, use last defined rho_solv
                if simpar[simpar.par == ('rho_solv_' + str(i))].empty:
                    rhosolv = simpar[simpar.par == ('rho_solv_' + str(last_defined_rho_solv))].iloc[0][1]
                else:
                    rhosolv = simpar[simpar.par == ('rho_solv_' + str(i))].iloc[0][1]
                    last_defined_rho_solv = i
                if fabs(rhosolv) > 1E-4:
                    rhosolv = rhosolv * 1E-6

                cbmat = (rhosolv - rhomin) / (rhomax - rhomin) * (cbmatmax - cbmatmin) + cbmatmin

            comments = general.extract_comments_from_file(self.spath + '/sim' + str(i) + '.dat', "#")
            simdata = pandas.read_csv(self.spath+'/sim' + str(i) + '.dat', sep='\s+', skip_blank_lines=True,
                                      comment='#', header=None)

            if len(simdata.columns) == 3:
                simdata.columns = ['Q', 'R', 'dR']
            elif len(simdata.columns) == 4:
                simdata.columns = ['Q', 'R', 'dR', 'dQ']
            simdata['dR'] = 0.0

            for index in simdata.index:
                ns = I * simdata['R'][index] * c1 * c2 * (simdata['Q'][index]) ** 2 * (
                            c3 + c4 * (simdata['Q'][index]) ** 2)
                nb = I * cbmat * c1 * c2 * (simdata['Q'][index]) ** 2 * (c3 + c4 * (simdata['Q'][index]) ** 2)
                dRoR = sqrt(ns + 2 * nb) / ns
                dR = simdata.iloc[index, 1] * dRoR  # see manuscript for details on calculation
                simdata['dR'][index] = dR
                # modify reflectivity within error bars
                simdata['R'][index] = simdata['R'][index] + normalvariate(0, 1) * dR

            simdata.to_csv(self.spath+'/sim' + str(i) + '.dat', sep=' ', index=None, header=None)
            general.add_comments_to_start_of_file(self.spath + '/sim' + str(i) + '.dat', comments)
            i += 1


class CGaReflInteractor(CRefl1DInteractor):
    def __init__(self, spath='.', mcmcpath='.', runfile='', load_state=True):
        super().__init__(spath, mcmcpath, runfile, load_state=load_state)

    def fnBackup(self, origin=None, target=None):
        if origin is None:
            origin = self.spath
        if target is None:
            target = self.spath + '/rsbackup'
        if not path.isdir(target):
            os.mkdir(target)
        for file in glob.glob(origin + r'/*.dat'):
            shutil.copy(file, target)
        for file in glob.glob(origin + r'/*.cc'):
            shutil.copy(file, target)
        for file in glob.glob(origin + r'/*.h'):
            shutil.copy(file, target)
        for file in glob.glob(origin + r'/*.o'):
            shutil.copy(file, target)
        for file in glob.glob(origin + r'/*.py'):
            shutil.copy(file, target)
        for file in glob.glob(origin + r'/*.pyc'):
            shutil.copy(file, target)
        shutil.copy(origin + '/Makefile', target)
        shutil.copy(origin + '/fit', target)

    def fnBackupSimdat(self):
        def backup_simdat(i):
            if not path.isfile('simbackup' + str(i) + '.dat'):
                pr = Popen(["cp", 'sim' + str(i) + '.dat', 'simbackup' + str(i) + '.dat'])
                pr.wait()
            else:
                pr = Popen(["cp", 'simbackup' + str(i) + '.dat', 'sim' + str(i) + '.dat'])
                pr.wait()

        i = 0
        while path.isfile('fit' + str(i) + '.dat'):
            backup_simdat(i)
            i += 1

    def fnLoadFileListAndChangeToLocal(self):
        # scans the setup.c file and creates a list of the filenames of the reflectivity data files, it also
        # modifies setup.c in this way that it loads copies of reflectivity dat files located in the working directory
        # with the file ending .mce, those files are the modified files by the statistical error analysis

        File = open(self.spath+'/setup.cc', "r")
        data = File.readlines()
        File.close()
        newdata = []
        filelist = []
        smatch1 = compile(r'fit_data.+?\"(.+?)\"', IGNORECASE | VERBOSE)
        smatch2 = compile(r'(fit_data.+?\").+?(\")', IGNORECASE | VERBOSE)
        smatch3 = compile(r'\s*//', IGNORECASE | VERBOSE)
        for line in data:
            if ('fit_data' in line) and (not smatch3.match(line)):  # scan for file loading but not comment lines
                filelist.append(smatch1.search(line).group(1))  # append path+filename to filelist
                newdata.append(smatch2.sub(r'\1' +
                                           path.split(filelist[-1])[-1] + r'.mce\2',
                                           line))  # modifiy setup.c line with new filename
            else:
                newdata.append(line)  # no file loading -> do not modify setup.c line
        file = open(self.spath+'/setup.cc', "w")  # write setup.c
        file.writelines(newdata)
        file.close()
        return filelist

    def fnLoadMolgroups(self, problem=None):
        diMolgroups = {}
        li = []

        file = open(self.spath + '/mol.dat')
        data = file.readlines()
        file.close()

        i = 0
        while i < len(data):
            tdata = (data[i]).split()  # read header that contains molgroup data
            diMolgroups[tdata[1]] = {}
            diMolgroups[tdata[1]].update({'headerdata': {}})
            diMolgroups[tdata[1]]['headerdata'].update({'Type': tdata[0]})
            diMolgroups[tdata[1]]['headerdata'].update({'ID': tdata[1]})
            for j in range(2, len(tdata), 2):
                diMolgroups[tdata[1]]['headerdata'].update({tdata[j]: tdata[j + 1]})

            i += 2  # skip header line for data columns
            zax = li[:]
            areaax = li[:]
            nslax = li[:]
            diMolgroups[tdata[1]].update({'zaxis': zax, 'areaaxis': areaax, 'nslaxis': nslax})

            while i < len(data):
                tline = (data[i]).split()
                if tline:
                    diMolgroups[tdata[1]]['zaxis'].append(float(tline[0]))
                    diMolgroups[tdata[1]]['areaaxis'].append(float(tline[1]))
                    diMolgroups[tdata[1]]['nslaxis'].append(float(tline[2]))
                    i += 1
                else:
                    break
            i += 1

        return diMolgroups

    def fnGetNumberOfModelsFromSetupC(self):
        file = open(self.runfile, "r")  # open setup.c
        data = file.readlines()
        file.close()
        smatch = compile(r'define\s+MODELS\s+(.+?)\n', IGNORECASE | VERBOSE)
        for line in data:  # search through setup.c
            if smatch.search(line):  # searching for MODELS constant
                i = smatch.search(line).group(1)
                return int(i)
        return 0

    def fnGetTaggedParameters(self):  # returns a list of the name and the
        file = open(self.spath + '/setup.cc')  # range + stepsize information of parameters
        data = file.readlines()  # which are tagged for displacement error
        file.close()  # analysis
        output = []
        for line in data:
            if '!rstag' in line:
                smatch = compile(r'pars_add\(pars.*?\"(.+?)\"\s*,.+?,(.+?),(.+?)\).+!rstag\s+(.+?)\s+(.+?)\s+(.+?)\s+!',
                                 IGNORECASE | VERBOSE)
                output.append(smatch.search(line).groups())
        return output

    def fnImportMCMCBestFit(self):  # imports best-fit from MCMC

        call(['cp', 'setup.c', 'setup.back'])
        call(['cp', 'setup.cc', 'setup.backcc'])

        self.fnLoadParameters()                             # Load Parameters and modify setup.cc

        # get list of parameters from setup.c/par.dat and sort by number of appearance in setup.cc
        li_parameters = list(self.diParameters.keys())
        li_parameters = sorted(li_parameters, key=lambda keyitem: self.diParameters[keyitem]['number'])

        # change setup.c to quasi fix all parameters
        li_addition = []
        for parameter in li_parameters:
            li_addition.append(('%s = %s;\n' %
                                (self.diParameters[parameter]['variable'], self.diParameters[parameter]['value'])))
        self.fnWriteConstraint2Runfile(li_addition)
        call(["rm", "-f", "setup.o"])
        call(["sync"])  # synchronize file system
        sleep(1)  # wait for system to clean up
        self.fnMake()  # compile changed setup.c
        call(['nice', './fit', '-el'])  # initiate write out par.dat via LM fit
        call(['cp', 'pop_bak.dat', 'pop.dat'])
        call(['mv', 'setup.back', 'setup.c'])
        call(['mv', 'setup.backcc', 'setup.cc'])
        call(["sync"])  # synchronize file system
        sleep(1)  # wait for system to clean up

    def fnLoadParameters(self):
        filename = self.spath + '/setup.cc'
        self.diParameters = {}

        # check wether an existing MCMC exists
        if path.isfile(self.mcmcpath + '/run.par'):
            print('Found ' + self.mcmcpath + '/run.par \n')
            print('Loading MCMC best-fit parameters ...')
            File = open(self.mcmcpath + '/run.par')
            data = File.readlines()
            File.close()
            tParValues = []
            for line in data:
                tParValues.append(float(line.split()[-1]))
            chisq = 0  # for the moment I cannot import chisq
        else:
            if path.isfile(self.spath + '/par.dat'):
                file = open(self.spath + '/par.dat')
                data = file.readlines()
                file.close()
            else:
                print('--------------------------------------')
                print('No par.dat found - Initializing fit.')
                print('--------------------------------------')
                self.fnMake()
                pr = Popen(["nice", "./fit", "-S", "-n", "1"])  # initial genetic run
                pr.wait()
                if path.isfile(self.spath + '/par.dat'):
                    file = open(self.spath + '/par.dat')
                    data = file.readlines()
                    file.close()
                else:
                    print('--------------------------------------')
                    print('Could not start up fit.')
                    print('--------------------------------------')
                    raise RuntimeError('Could not start up fit.')
            while data[-1][0] == '#':  # delete comments in last lines
                data = data[:-1]
            tParValues = (data[-1]).split()[2:]  # best fit values from last row
            chisq = float((data[-1]).split()[1])  # chi squared from second column last row

        File = open(filename)
        setupdata = File.readlines()
        File.close()
        # reminder .+? match all character
        smatch = compile(r'pars_add\(pars.*?\"(.+?)\"\s*,.+?\((.+?)\)\s*,(.+?),(.+?)\)', IGNORECASE | VERBOSE)
        smatch2 = compile(r'\s*//', IGNORECASE | VERBOSE)
        for i in range(len(tParValues)):  # iterate over read parameters
            for j in range(len(setupdata)):  # scan setup.c for the first parameter
                setupline = setupdata[j]  # initialization and fetch that name
                if smatch.search(setupline) and (not smatch2.match(setupline)):  # plus its fit ranges
                    sParName = smatch.search(setupline).group(1)
                    sVarName = smatch.search(setupline).group(2)
                    flowerlimit = float(smatch.search(setupline).group(3))
                    fupperlimit = float(smatch.search(setupline).group(4))
                    del setupdata[j]  # delete the line where parameter found
                    break
            else:
                print('--------------------------------------')
                print('Parameters do not match in setup and parameter file ')  # no parameter in setup.c left!
                print('--------------------------------------')
                raise RuntimeError('Mismatch between setup file and par.dat.')
            if sParName in self.diParameters:
                print('--------------------------------------')
                print('The parameter %s is defined twice in the garefl setup file.' % sParName)
                print('This is not supported by rs.py.')
                print('--------------------------------------')
                raise RuntimeError('Doubled parameter in setup file.')
            ipar = int(i)  # parameter number
            frelval = float(tParValues[i])  # relative par value is stored in file
            fvalue = (fupperlimit - flowerlimit) * frelval + flowerlimit  # calculate absolute par value
            self.diParameters[sParName] = {'number': ipar, 'variable': sVarName, 'lowerlimit': flowerlimit,
                                      'upperlimit': fupperlimit, 'relval': frelval, 'value': fvalue}

        # There was a time when all parameters were given as relative numbers between the bounds
        # This tests if any parameter is out of bounds to determine whether we have relative or
        # absolute parameter values
        # Code disabled, as this should not occur anymore
        # bAbsoluteParameters = False
        # for sParName in list(diParameters.keys()):
        #     if diParameters[sParName]['relval'] > 1:
        #         bAbsoluteParameters = True
        # if bAbsoluteParameters:
        #     for sParName in list(diParameters.keys()):
        #         diParameters[sParName]['value'] = diParameters[sParName]['relval']
        #         diParameters[sParName]['relval'] = (diParameters[sParName]['value'] -
        #                                             diParameters[sParName]['lowerlimit']) / \
        #                                            (diParameters[sParName]['upperlimit'] -
        #                                             diParameters[sParName]['lowerlimit'])

        # redo refl1d convention to multiply small parameters (nSLD, background) by 1E6
        # in case bounds were specified in 1E-6 units (garefl)
        for sParName in list(self.diParameters.keys()):
            if self.diParameters[sParName]['value'] < self.diParameters[sParName]['lowerlimit'] or \
                    self.diParameters[sParName]['value'] > self.diParameters[sParName]['upperlimit']:
                self.diParameters[sParName]['value'] /= 1E6

        return self.diParameters, chisq

    def fnLoadStatData(self, dSparse=0., rescale_small_numbers=True, skip_entries=None):
        if skip_entries is None:
            skip_entries = []
        # Load a file like iSErr or SErr into
        # the self.liStatResult object
        sStatResultHeader = ''
        liStatResult = []
        if not (path.isfile(self.mcmcpath+'/sErr.dat') or path.isfile(self.mcmcpath+'/isErr.dat')):
            # call super method to recreate sErr.dat from MCMC path, default MCMC runfile for garefl is 'run'
            store_runfile = self.runfile
            self.runfile = 'run'
            super(CRefl1DInteractor, self).fnLoadStatData(dSparse)
            self.runfile = store_runfile

        return self.fnLoadsErr()

    def fnMake(self):
        # make setup.c and print sys output
        call(["rm", "-f", self.spath + "/setup.o"])
        call(["sync"])  # synchronize file system
        sleep(1)  # wait for system to clean up
        call(['make', '-C', self.spath])

    def fnReplaceParameterLimitsInSetup(self, sname, flowerlimit, fupperlimit):
        # scans setup.c file for parameter with name sname and replaces the
        # lower and upper fit limits by the given values
        file = open(self.spath+'/setup.cc', 'r+')
        data = file.readlines()
        file.close()
        smatch = compile(r'(pars_add\(pars.*?\"' + sname + '.+?,.+?,).+?(,).+?(\))', IGNORECASE | VERBOSE)
        newdata = []
        for line in data:
            newdata.append(smatch.sub(r'\1 ' + str(flowerlimit) + r'\2 '
                                      + str(fupperlimit) + r'\3', line))

        file = open(self.spath+'/setup.cc', 'w')
        file.writelines(newdata)
        file.close()

    def fnRestoreFileList(self, filelist):  # not used
        file = open(self.spath+'/setup.cc', "r")
        data = file.readlines()
        file.close()
        newdata = []
        smatch = compile(r'(fit_data.+?\").+?(\")', IGNORECASE | VERBOSE)
        for line in data:
            if 'void constr_models' in line:
                newdata.append(smatch.sub(r'\1' + filelist[0] + r'\2', line))
                pr = Popen(['rm', '-f', path.split(filelist[0])[-1] + '.mce'])
                pr.wait()
                del filelist[0]
            else:
                newdata.append(line)
        file = open(self.spath+'/setup.cc', "w")
        file.writelines(newdata)
        file.close()

    def fnRestoreFitProblem(self):
        from refl1d import garefl
        # iNumberOfModels = self.fnGetNumberOfModelsFromSetupC()  # how many profiles to analyze?
        if path.isfile(self.spath + '/model.so'):
            self.fnMake()
            problem = garefl.load(self.spath + '/model.so')
        else:
            print('No problem to reload.')
            problem = None
        return problem

    def fnRestoreMolgroups(self, problem):
        problem.active_model.fitness.output_model()
        diMolgroups = self.fnLoadMolgroups()
        return diMolgroups

    @staticmethod
    def fnRestoreSmoothProfile(M):
        z, rho, irho = M.fitness.smooth_profile()
        return z, rho, irho

    def fnRunMCMC(self, burn, steps, batch=False, compile_setup=True):
        if compile_setup:
            self.fnMake()
        CRefl1DInteractor.fnRunMCMC(self, burn=burn, steps=steps, batch=batch)

    def fnSaveMolgroups(self, problem):
        # should call the setup.cc function that saves mol.dat
        # TODO: Needs testing
        problem.active_model.fitness.output_model()

    def fnSimulateData(self, diAddition):
        # change dir of parname:parvalue into list of expressions to be inserted into setup.cc
        liExpression = []
        for parameter in diAddition:
            parvalue = diAddition[parameter]
            strchange = ''
            parvaluefinal = parvalue
            if ('rho' in parameter) and fabs(parvalue) > 1E-4:
                parvaluefinal = parvalue * 1E-6
                strchange = ' => changed to ' + str(parvaluefinal)
            elif ('background' in parameter) and parvalue > 1E-4:
                # The issue with the background is two different uses
                # some setup.cc use negative numbers as parameter values and then
                # compute background=10^value, others use positive values and then
                # compute background=value *1E-6. This should catch it all
                parvaluefinal = parvalue * 1E-6
                strchange = ' => changed to ' + str(parvaluefinal)

            print(str(parameter) + ' ' + str(parvalue) + strchange)
            liExpression.append(('%s = %s;\n' % (self.diParameters[parameter]['variable'], parvaluefinal)))

        shutil.copy(self.spath + '/setup.cc', self.spath + '/setup_bak.cc')
        self.fnWriteConstraint2Runfile(liExpression)                        # change runfile to quasi-fix all parameters
        self.fnMake()                                                       # compile changed setup.c
        call(["rm", "-f", self.spath + "/mol.dat"])
        call([self.spath + "/fit", "-g"])                                   # write out profile.dat and fit.dat
        call(["sync"])                                                      # synchronize file system
        sleep(1)
        # restore setup file without quasi-fixed parameters
        shutil.copy(self.spath + '/setup_bak.cc', self.spath + '/setup.cc')
        call(["rm", "-f", self.spath + '/setup_bak.cc'])

        i = 0
        while path.isfile(self.spath+'/fit' + str(i) + '.dat'):
            simdata = pandas.read_csv(self.spath+'/fit' + str(i) + '.dat', sep='\s+', header=None,
                                      names=['Q', 'dQ', 'R', 'dR', 'fit'],
                                      skip_blank_lines=True, comment='#')
            del simdata['dQ']
            del simdata['R']
            simdata = simdata.rename(columns={'fit': 'R'})
            simdata = simdata[['Q', 'R', 'dR']]
            simdata.to_csv(self.spath+'/sim' + str(i) + '.dat', sep=' ', index=None, header=False)
            i += 1

    @staticmethod
    def fnTruncateParDat():
        if path.exists("par.dat"):  # check if something is there to truncate
            file = open("par.dat", 'r')
            data = file.readlines()
            file.close()

            iFirstDataLine = -1
            iLastDataLine = -1

            for i in range(int(len(data) / 2)):  # find position of first and last line
                if (data[i][0] != '#') and (iFirstDataLine == -1):  # with data
                    iFirstDataLine = i
                if (data[len(data) - 1 - i][0] != '#') and (iLastDataLine == -1):
                    iLastDataLine = len(data) - 1 - i
                if (iFirstDataLine != -1) and (iLastDataLine != -1):
                    break

            if (iFirstDataLine != -1) and (iLastDataLine != -1):  # cut out everything between those two
                data1 = data[:iFirstDataLine + 1]  # lines
                data2 = data[iLastDataLine:]
                file = open("par.dat", 'w')
                file.writelines(data1)
                file.writelines(data2)
                file.close()

    @staticmethod
    def fnWriteOutGareflModel():
        pr = Popen(["./fit", "-g"])  # write out profile.dat and fit.dat
        pr.wait()

    def fnWriteConstraint2Runfile(self, liExpression):
        # Writes a list of expressions at the beginning of
        # the constraint section in setup.c
        File = open(self.spath+'/setup.cc', "r")  # open setup.c
        data = File.readlines()
        File.close()
        newdata = []

        for i in range(len(data)):  # need integer index for slicing
            if 'void constr_models(' in data[i]:  # searching for function declaration
                newdata = data[:i + 2] + liExpression + data[i + 2:]  # inserting Expression two lines later,
                # after the initial '{' of the function

                break  # break loop, only one insertion
        File = open(self.spath+'/setup.cc', "w")  # write changed setup.cc
        File.writelines(newdata)
        File.close()


class CSASViewInteractor(CBumpsInteractor):
    def __init__(self, spath='.', mcmcpath='.', runfile='', load_state=True):
        super().__init__(spath, mcmcpath, runfile, load_state=load_state)
