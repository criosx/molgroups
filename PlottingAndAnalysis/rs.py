# issues: ISIS data files cannot contain headerline without specifier as they normally do, this will
#        lead to MC not starting up

from __future__ import print_function
from math import fabs, pow, floor, ceil, sqrt, log10
from numpy import subtract, minimum, maximum, average, array
from operator import itemgetter
from os import getcwd, remove, rename, replace, path, kill, devnull
from random import seed, normalvariate, random
from re import VERBOSE, IGNORECASE, compile
from scipy import stats, special
from shutil import rmtree
from sys import argv, exit
from subprocess import call, Popen
from time import time, gmtime, sleep, ctime
from bumps.dream.stats import parse_var
import zipfile
import numpy
import pandas
import re


class CDataInteractor():
    def __init__(self, spath='.', mcmcpath='.', runfile=''):
        # script path, MCMC path, runfile (script file)
        self.spath = spath
        self.mcmcpath = mcmcpath
        self.runfile = runfile

    def fnLoadMolgroups(self):
        diMolgroups = {}
        li = []

        file = open(self.spath + '/mol.dat')
        data = file.readlines()
        file.close()

        i = 0
        while 1:
            tdata = (data[i]).split()  # read header that contains molgroup data
            diMolgroups[tdata[1]] = {}
            diMolgroups[tdata[1]].update({'headerdata': {}})
            diMolgroups[tdata[1]]['headerdata'].update({'Type': tdata[0]})
            diMolgroups[tdata[1]]['headerdata'].update({'ID': tdata[1]})
            j = 2
            while 1:
                diMolgroups[tdata[1]]['headerdata'].update({tdata[j]: tdata[j + 1]})
                j += 2
                if j == len(tdata):
                    break

            i += 2  # skip header line for data columns
            zax = li[:]
            areaax = li[:]
            nslax = li[:]
            diMolgroups[tdata[1]].update({'zaxis': zax, 'areaaxis': areaax, 'nslaxis': nslax})

            while 1:
                if i >= len(data):
                    break
                tline = (data[i]).split()
                if tline:
                    diMolgroups[tdata[1]]['zaxis'].append(float(tline[0]))
                    diMolgroups[tdata[1]]['areaaxis'].append(float(tline[1]))
                    diMolgroups[tdata[1]]['nslaxis'].append(float(tline[2]))
                    i += 1
                else:
                    break

            i += 1
            if i >= len(data):
                break
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
    def fnSaveSingleColumnsFromStatDict(sFilename, data, skipentries=[]):

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


# Garefl methods will be used if a storage directory for a Markov Chain Monte Carlo (MCMC)
# error analysis cannot be found.
# The MCMC directory is called 'MCMC'

class CBumpsInteractor(CDataInteractor):
    def __init__(self, spath='.', mcmcpath='.', runfile=''):
        super().__init__(spath, mcmcpath, runfile)

        # patterns to extract parmeter information from .err results files
        self.VAR_PATTERN1 = re.compile(r"""
            ^.(?P<parname>.+?)\ =\ 
            (?P<best>[0-9.eE+-]+?)\ in\ \[           
            (?P<lowbound>[0-9.eE+-]+?),
            (?P<highbound>[0-9.eE+-]+?)\]
            .*?
            $""", re.VERBOSE)

        self.VAR_PATTERN2 = re.compile(r"""
           ^\ *
           (?P<parnum>[0-9]+)\ +
           (?P<parname>.+?)\ +
           (?P<mean>[0-9.-]+?)
           \((?P<err>[0-9]+)\)
           (e(?P<exp>[+-]?[0-9]+))?\ +
           (?P<median>[0-9.eE+-]+?)\ +
           (?P<best>[0-9.eE+-]+?)\ +
           \[\ *(?P<lo68>[0-9.eE+-]+?)\ +
           (?P<hi68>[0-9.eE+-]+?)\]\ +
           \[\ *(?P<lo95>[0-9.eE+-]+?)\ +
           (?P<hi95>[0-9.eE+-]+?)\]
           \ *$
           """, re.VERBOSE)

        self.VAR_PATTERN_CHISQ = re.compile(r"""
            ^\[chisq= 
            (?P<value>[0-9.eE+-]+?)\(
            .*?
            $""", re.VERBOSE)

        self.VAR_PATTERN_OVERALLCHISQ = re.compile(r"""
            ^\[overall\ chisq= 
            (?P<value>[0-9.eE+-]+?)\(
            .*?
            $""", re.VERBOSE)

    def fnLoadMCMCResults(self):

        import bumps.dream.state

        state = bumps.dream.state.load_state(self.mcmcpath + '/' + self.runfile)
        state.mark_outliers()  # ignore outlier chains

        # load Parameter
        data = self.fnLoadSingleColumns(self.mcmcpath + '/' + self.runfile + '.par', header=False,
                                        headerline=['parname', 'bestfitvalue'])
        lParName = []
        vars = []
        i = 0
        for parname in data['parname']:
            lParName.append(parname)
            vars.append(i)
            i += 1

        draw = state.draw(1, vars, None)
        points = draw.points
        logp = draw.logp

        return points, lParName, logp

    # LoadStatResults returns a list of variable names, a logP array, and a numpy.ndarray
    # [values,var_numbers].
    def fnLoadParameters(self):
        def parse_var(line):
            """
            Parse a line returned by format_vars back into the statistics for the
            variable on that line.
            """

            m1 = self.VAR_PATTERN1.match(line)
            m2 = self.VAR_PATTERN2.match(line)
            results = None
            if m2:
                exp = int(m2.group('exp')) if m2.group('exp') else 0
                results = {'index': int(m2.group('parnum')), 'name': m2.group('parname'),
                           'mean': float(m2.group('mean')) * 10 ** exp, 'median': float(m2.group('median')),
                           'best': float(m2.group('best')), 'p68': (float(m2.group('lo68')), float(m2.group('hi68'))),
                           'p95': (float(m2.group('lo95')), float(m2.group('hi95')))}
            elif m1:
                results = {'name': m1.group('parname'), 'best': float(m1.group('best')),
                           'bounds': (float(m1.group('lowbound')), float(m1.group('highbound')))}
            return results

        # this code is from bumps.util parse_errfile with modifications
        diParameters = {}
        chisq = []          # not used here
        overall = None
        filename = self.mcmcpath + '/' + self.runfile + '.err'
        with open(filename) as fid:
            for line in fid:
                if line.startswith("[overall"):
                    m = self.VAR_PATTERN_OVERALLCHISQ.match(line)
                    overall = float(m.group('value'))
                    continue

                if line.startswith("[chisq"):
                    m = self.VAR_PATTERN_CHISQ.match(line)
                    chisq.append(float(m.group('value')))
                    continue

                p = parse_var(line)
                if p is not None:
                    parametername = p['name']
                    if parametername not in diParameters:
                        diParameters[parametername] = {}
                    if 'bounds' in p:
                        diParameters[parametername]['lowerlimit'] = p['bounds'][0]
                        diParameters[parametername]['upperlimit'] = p['bounds'][1]
                        diParameters[parametername]['value'] = p['best']
                        diParameters[parametername]['relval'] = (p['best'] - p['bounds'][0]) / \
                                                                (p['bounds'][1] - p['bounds'][0])
                        diParameters[parametername]['variable'] = p['name']
                    else:
                        diParameters[parametername]['number'] = p['index']
                        diParameters[parametername]['error'] = (p['p68'][1] - p['p68'][0])/2

        if overall is None:
            overall = chisq[0]

        return diParameters, overall

    def fnLoadStatData(self, dSparse=0, rescale_small_numbers=True, skip_entries=[]):
        if path.isfile(self.mcmcpath + '/sErr.dat') or path.isfile(self.mcmcpath + '/isErr.dat'):
            diStatRawData = self.fnLoadsErr()
        else:
            points, lParName, logp = self.fnLoadMCMCResults()

            diStatRawData = {'Parameters': {}}
            diStatRawData['Parameters']['Chisq'] = {}  # TODO: Work on better chisq handling
            diStatRawData['Parameters']['Chisq']['Values'] = []
            for parname in lParName:
                diStatRawData['Parameters'][parname] = {}
                diStatRawData['Parameters'][parname]['Values'] = []

            seed()
            for j in range(len(points[:, 0])):
                if dSparse == 0 or (dSparse > 1 and j < dSparse) or (1 > dSparse > random()):
                    diStatRawData['Parameters']['Chisq']['Values'].append(logp[j])
                    for i, parname in enumerate(lParName):
                        # TODO: this is a hack because Paul does not scale down after scaling up
                        # Rescaling disabled for bumps/refl1d analysis to achieve consistency
                        # if ('rho_' in parname or 'background' in parname) and rescale_small_numbers:
                        #     points[j, i] *= 1E-6
                        diStatRawData['Parameters'][parname]['Values'].append(points[j, i])

            self.fnSaveSingleColumnsFromStatDict(self.mcmcpath + '/sErr.dat', diStatRawData['Parameters'], skip_entries)

        return diStatRawData

    def fnRestoreFitProblem(self):
        from bumps.fitproblem import load_problem
        problem = load_problem(self.spath + '/' + self.runfile + '.py')
        return problem

    def fnRestoreMolgroups(self, problem):
        diMolgroups = self.fnLoadMolgroups()
        return diMolgroups

    @staticmethod
    def fnRestoreSmoothProfile(M):
        # TODO: Decide what and if to return SLD profile for Bumps fits
        z, rho, irho = [], [], []
        return z, rho, irho


# Refl1D methods will be used if a storage directory for a Markov Chain Monte Carlo (MCMC)
# error analysis are found.
# The MCMC directory is called 'MCMC'
# The refl1d script name has to be run.py.
class CRefl1DInteractor(CBumpsInteractor):
    def __init__(self, spath='.', mcmcpath='.', runfile=''):
        super().__init__(spath, mcmcpath, runfile)

        # patterns to extract parmeter information from .err results files
        self.VAR_PATTERN1 = re.compile(r"""
            ^\[(?P<parnum>[0-9]+)\]\ =\ Parameter\(
            (?P<best>[0-9.eE+-]+?),\ name='
            (?P<parname>.+?)',\ bounds=\(
            (?P<lowbound>[0-9.eE+-]+?),
            (?P<highbound>[0-9.eE+-]+?)\)
            .*?
            $""", re.VERBOSE)

    @staticmethod
    def fnRestoreSmoothProfile(M):
        z, rho, irho = M.fitness.smooth_profile()
        return z, rho, irho


class CGaReflInteractor(CRefl1DInteractor):
    def __init__(self, spath='.', mcmcpath='.', runfile=''):
        super().__init__(spath, mcmcpath, runfile)

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

    def fnLoadParameters(self):
        filename = self.spath + '/' + self.runfile
        diParameters = {}

        # check wether an existing MCMC exists
        if path.isfile(self.spath + '/run.par'):
            print('Found ' + self.spath + '/run.par \n')
            print('Using MCMC best-fit parameters ...')
            File = open(self.spath + '/run.par')
            data = File.readlines()
            File.close()
            tParValues = []
            for line in data:
                tParValues.append(float(line.split()[-1]))
            chisq = 0  # for the moment I cannot import chisq
        elif zipfile.is_zipfile(self.spath):  # unpack from zip-file
            File = zipfile.ZipFile(self.spath, 'r')
            File.extractall(path.dirname(self.spath))
            File.close()
            spath2 = path.dirname(self.spath)
            if spath2 == '':
                spath2 = '.'
            File = open(spath2 + '/' + path.basename(self.spath)[:-4] + '/run.par')
            data = File.readlines()
            File.close()
            tParValues = []
            for line in data:
                tParValues.append(float(line.split()[-1]))
            chisq = 0  # for the moment I cannot import chisq
            rmtree(spath2 + '/' + path.basename(self.spath)[:-4])
            if path.isdir(spath2 + '/__MACOSX'):
                rmtree(spath2 + '/__MACOSX')
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
            if sParName in diParameters:
                print('--------------------------------------')
                print('The parameter %s is defined twice in the garefl setup file.' % sParName)
                print('This is not supported by rs.py.')
                print('--------------------------------------')
                raise RuntimeError('Doubled parameter in setup file.')
            ipar = int(i)  # parameter number
            frelval = float(tParValues[i])  # relative par value is stored in file
            fvalue = (fupperlimit - flowerlimit) * frelval + flowerlimit  # calculate absolute par value
            diParameters[sParName] = {'number': ipar, 'variable': sVarName, 'lowerlimit': flowerlimit,
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
        for sParName in list(diParameters.keys()):
            if diParameters[sParName]['value'] < diParameters[sParName]['lowerlimit'] or \
                    diParameters[sParName]['value'] > diParameters[sParName]['upperlimit']:
                diParameters[sParName]['value'] /= 1E6

        return diParameters, chisq

    def fnLoadStatData(self, dSparse=0., rescale_small_numbers=True, skip_entries=[]):

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

    def fnRestoreFitProblem(self):
        from refl1d import garefl
        # iNumberOfModels = self.fnGetNumberOfModelsFromSetupC()  # how many profiles to analyze?
        self.fnMake()
        problem = garefl.load(self.spath + '/model.so')
        return problem

    def fnRestoreMolgroups(self, problem):
        problem.active_model.fitness.output_model()
        diMolgroups = self.fnLoadMolgroups()
        return diMolgroups

    @staticmethod
    def fnRestoreSmoothProfile(M):
        z, rho, irho = M.fitness.smooth_profile()
        return z, rho, irho


class CMolStat:
    def __init__(self, fitsource='refl1d', spath='.', mcmcpath='.', runfile='run'):
        """

        """
        self.diParameters = {}
        # Dictionary with all the parameters
        # self.diParameters data structure:
        # dictionary: sparameter : dictionary
        #                          'number'    : int       # order of initialization in ga_refl
        #                          'lowerlimit'  : float   # constraint range lower limit
        #                          'upperlimit'  : float   # constraint range upper limit
        #                          'value'  : float        # absolute value
        #                          'error'  : float        # error derived from covariance matrix
        #                          'relval' : float        # relative value between 0 and 1 in terms of the constraints
        #                          'variable: string       # associated ga_refl variable
        self.diMolgroups = {}
        self.diStatResults = {}
        # dictionary of statistical results for various routines
        # from fnAnalyzeStatFile
        #   'Parameters' key contains dictionary
        #       parameter names are keys contain dictionary
        #           'Values' key contains ordered list of all MC parameters
        #           'LowPerc' contains lower percentile
        #           'Median' contains median
        #           'HighPerc' contains higher percentile
        #   'MaxParameterLength' contains length of longest parameter name
        # from fnRecreateStatFile
        #   'NumberOfStatValues' number of MC points
        #   'nSLDProfiles' key contains list of lists of lists
        #      of all profiles 'nSLDProfiles'[MCiteration][model][z,rho]
        # from fnLoadMolgroups
        #   'Molgroups' contains information about molecular groups
        #     of all iterations 'Molgroups'[MCiteration] {molgroupname {'zaxis' [z]
        #                                                               'areaxis' [area]
        #                                                               'nslaxis' [nsl]
        #                                                               'headerdata' {'property' value}
        #                                                               'property1' value
        #                                                               'property2' value}}

        self.fitsource = fitsource                  # define fitting software
        self.spath = spath
        self.mcmcpath = mcmcpath
        self.runfile = runfile

        self.liStatResult = []                      # a list that contains the isErr.dat or sErr.dat file line by line
        self.sStatResultHeader = ''                 # headerline from isErr.dat or sErr.dat
        self.sStatFileName = ''                     # Name of the statistical File
        self.chisq = 0.
        self.fMolgroupsStepsize = 0.
        self.iMolgroupsDimension = 0
        self.fMolgroupsNormArea = 0
        # check for system and type of setup file

        self.Interactor = None
        if self.fitsource == 'bumps':
            self.Interactor = CBumpsInteractor(spath, mcmcpath, runfile)
        elif self.fitsource == 'refl1d':
            self.Interactor = CRefl1DInteractor(spath, mcmcpath, runfile)
        elif self.fitsource == 'garefl':
            self.Interactor = CGaReflInteractor(spath, mcmcpath, runfile)

    def fnAnalyzeStatFile(self, fConfidence=-1, sparse=0):  # summarizes stat file

        self.fnLoadStatData(sparse)                     # Load data from file into list
        self.fnLoadParameters()                         # Load Parameters for limits

        if fConfidence > 1:
            fConfidence = 1
        if fConfidence < 0:
            fConfidence = special.erf(-1 * fConfidence / sqrt(2))

        fLowerPercentileMark = 100.0 * (1 - fConfidence) / 2
        fHigherPercentileMark = (100 - fLowerPercentileMark)

        iNumberOfMCIterations = self.diStatResults['NumberOfStatValues']  # how many iterations already done
        iMaxParameterNameLength = self.diStatResults['MaxParameterLength']
        print('Analysis of current MC simulation ...')
        print('Number of iterations: %(ni)d' % {'ni': iNumberOfMCIterations})

        fMaxConvergence = 0
        iHistoryLength = 5

        for element in sorted(list(self.diParameters.keys()),
                              key=lambda sParameter: self.diParameters[sParameter]['number']):

            fLowPerc = stats.scoreatpercentile \
                (self.diStatResults['Parameters'][element]['Values'], fLowerPercentileMark)  # Calculate Percentiles
            fMedian = stats.scoreatpercentile \
                (self.diStatResults['Parameters'][element]['Values'], 50.)
            fHighPerc = stats.scoreatpercentile \
                (self.diStatResults['Parameters'][element]['Values'], fHigherPercentileMark)

            self.diStatResults['Parameters'][element]['LowPerc'] = fLowPerc
            self.diStatResults['Parameters'][element]['Median'] = fMedian
            self.diStatResults['Parameters'][element]['HighPerc'] = fHighPerc

            flowerlimit = self.diParameters[element]['lowerlimit']
            fupperlimit = self.diParameters[element]['upperlimit']
            temp = abs(fupperlimit - flowerlimit)

            sGraphOutput = '['
            itemp1 = int((fLowPerc - flowerlimit) / temp * 10 + 0.5)
            itemp2 = int((fMedian - flowerlimit) / temp * 10 + 0.5)
            itemp3 = int((fHighPerc - flowerlimit) / temp * 10 + 0.5)

            for i in range(11):
                s1 = ' '
                if itemp1 == i or itemp3 == i: s1 = '|'
                if itemp2 == i:
                    if s1 == '|':
                        s1 = '+'
                    else:
                        s1 = '-'
                sGraphOutput += s1
            sGraphOutput += ']'

            if (fLowPerc - flowerlimit) < temp * 0.01:
                self.diStatResults['Parameters'][element]['LowerLimitCollision'] = True
                sGraphOutput = '#' + sGraphOutput[1:]
            else:
                self.diStatResults['Parameters'][element]['LowerLimitCollision'] = False

            if (fupperlimit - fHighPerc) < temp * 0.01:
                self.diStatResults['Parameters'][element]['UpperLimitCollision'] = True
                sGraphOutput = sGraphOutput[:-1] + '#'
            else:
                self.diStatResults['Parameters'][element]['UpperLimitCollision'] = False

            if iNumberOfMCIterations > iHistoryLength:  # at least five iterations

                fLowPercHistory = []
                fMedianHistory = []
                fHighPercHistory = []

                for i in range((iHistoryLength - 1) * (-1), 0):
                    fLowPercHistory.append(
                        stats.scoreatpercentile(self.diStatResults['Parameters'][element]['Values'][:i],
                                                fLowerPercentileMark))
                    fMedianHistory.append(
                        stats.scoreatpercentile(self.diStatResults['Parameters'][element]['Values'][:i], 50))
                    fHighPercHistory.append(
                        stats.scoreatpercentile(self.diStatResults['Parameters'][element]['Values'][:i],
                                                fHigherPercentileMark))

                fLowPercHistory.append(fLowPerc)
                fMedianHistory.append(fMedian)
                fHighPercHistory.append(fHighPerc)

                fLowPercAverage = average(fLowPercHistory)
                fMedianAverage = average(fMedianHistory)
                fHighPercAverage = average(fHighPercHistory)

                fSigma = []
                for i in range(iHistoryLength):
                    fSigma.append(abs(fLowPercHistory[i] - fLowPercAverage) / temp)
                    fSigma.append(abs(fMedianHistory[i] - fMedianAverage) / temp)
                    fSigma.append(abs(fHighPercHistory[i] - fHighPercAverage) / temp)

                fMaxConvergence = max(fMaxConvergence, max(fSigma))

                fLowPercConv = abs(
                    fLowPerc - stats.scoreatpercentile(self.diStatResults['Parameters'][element]['Values'][:-1],
                                                       fLowerPercentileMark)) / temp
                fMedianConv = abs(
                    fMedian - stats.scoreatpercentile(self.diStatResults['Parameters'][element]['Values'][:-1],
                                                      50)) / temp
                fHighPercConv = abs(
                    fHighPerc - stats.scoreatpercentile(self.diStatResults['Parameters'][element]['Values'][:-1],
                                                        fHigherPercentileMark)) / temp

                self.diStatResults['Parameters'][element]['LowPercConv'] = fLowPercConv
                self.diStatResults['Parameters'][element]['MedianConv'] = fMedianConv
                self.diStatResults['Parameters'][element]['HighPercConv'] = fHighPercConv

            else:
                fMaxConvergence = 1e6
                fLowPercConv = fMedianConv = fHighPercConv = 0

            sPrintString = '%(el)'
            sPrintString += str(iMaxParameterNameLength)
            sPrintString += 's  %(sg)s  [%(ll)10.4g,%(ul)10.4g]  [%(lp)10.4g(%(lpc).3f), %(m)10.4g(%(mc).3f), ' \
                            '%(hp)10.4g(%(hpc).3f)] (-%(ld)10.4g, +%(hd)10.4g)'

            print(sPrintString %
                  {'el': element, 'll': flowerlimit, 'ul': fupperlimit, 'lp': fLowPerc,
                   'lpc': fLowPercConv, 'ld': (fMedian - fLowPerc), 'm': fMedian,
                   'mc': fMedianConv, 'hd': (fHighPerc - fMedian), 'hp': fHighPerc,
                   'hpc': fHighPercConv, 'sg': sGraphOutput})

        self.diStatResults['Convergence'] = fMaxConvergence
        print('Maximum deviation from average over last %(iHL)d iterations: %(maxc).4f' %
              {'iHL': iHistoryLength, 'maxc': fMaxConvergence})
        print('Confidence level: %(fCL).4f' % {'fCL': fConfidence})

    # -------------------------------------------------------------------------------

    def fnBackup(self):
        if not path.isdir('rsbackup'):  # create backup dir
            pr = Popen(["mkdir", "rsbackup"])
            pr.wait()
        pr = Popen(["cp", "pop.dat", "rsbackup/"])  # save all data that will
        pr.wait()
        pr = Popen(["cp", "par.dat", "rsbackup/"])  # change during fit
        pr.wait()
        pr = Popen(["cp", "covar.dat", "rsbackup/"])
        pr.wait()
        call("cp fit*.dat rsbackup/", shell=True)
        pr = Popen(["cp", "fit", "rsbackup/"])
        pr.wait()
        call("cp model*.dat rsbackup/", shell=True)
        pr = Popen(["cp", "pop_bak.dat", "rsbackup/"])
        pr.wait()
        call("cp profile*.dat rsbackup/", shell=True)
        call("cp " + self.setupfilename + " rsbackup/", shell=True)
        pr = Popen(["cp", "setup.o", "rsbackup/"])
        pr.wait()

    # -------------------------------------------------------------------------------

    def fnCalculateMolgroupProperty(self, fConfidence):
        def fnFindMaxFWHM(lList):
            imax = 0
            maxvalue = 0
            for i in range(len(lList)):
                if lList[i] > maxvalue:
                    maxvalue = lList[i]
                    imax = i

            ifwhmplus = imax
            ifwhmminus = imax
            while lList[ifwhmplus] > (maxvalue / 2):
                ifwhmplus = ifwhmplus + 1
                if ifwhmplus == len(lList):
                    ifwhmplus = ifwhmplus - 1
                    break
            while lList[ifwhmminus] > (maxvalue / 2):
                ifwhmminus = ifwhmminus - 1
                if ifwhmminus == (-1):
                    ifwhmplus = 0
                    break
            return imax, maxvalue, ifwhmminus, ifwhmplus

        try:
            self.diStatResults = self.fnLoadObject(self.mcmcpath + '/StatDataPython.dat')
            print('Loaded statistical data from StatDataPython.dat')
        except IOError:
            print('Failure to load StatDataPython.dat.')
            print('Recreate statistical data from sErr.dat.')
            self.fnRecreateStatistical()

        # Try to import any fractional protein profiles stored in envelopefrac1/2.dat
        # after such an analysis has been done separately
        try:
            pdf_frac1 = pandas.read_csv(self.mcmcpath + '/envelopefrac1.dat', sep=' ')
            pdf_frac2 = pandas.read_csv(self.mcmcpath + '/envelopefrac2.dat', sep=' ')

            for mcmc_iter in range(len(self.diStatResults['Molgroups'])):
                striter = 'iter' + str(mcmc_iter)
                self.diStatResults['Molgroups'][mcmc_iter]['frac1'] = {}
                self.diStatResults['Molgroups'][mcmc_iter]['frac2'] = {}
                self.diStatResults['Molgroups'][mcmc_iter]['frac1']['zaxis'] = pdf_frac1['zaxis'].tolist()
                self.diStatResults['Molgroups'][mcmc_iter]['frac1']['areaaxis'] = pdf_frac1[striter].tolist()
                self.diStatResults['Molgroups'][mcmc_iter]['frac2']['zaxis'] = pdf_frac2['zaxis'].tolist()
                self.diStatResults['Molgroups'][mcmc_iter]['frac2']['areaaxis'] = pdf_frac2[striter].tolist()
        except IOError:
            print('Did not find any fractional envelopes ...')

        diResults = {}
        l_molgroups = list(self.diStatResults['Molgroups'][0].keys())
        for s_molgroup in l_molgroups:  # create results for individual molgroups
            diResults[s_molgroup + '_COM'] = []
            diResults[s_molgroup + '_INT'] = []
            diResults[s_molgroup + '_AVG'] = []

        for mcmc_iter in range(len(self.diStatResults['Molgroups'])):  # cycle over all MC iterations
            mgdict = self.diStatResults['Molgroups'][mcmc_iter]
            for s_molgroup in list(mgdict.keys()):  # cycle over all molgroups

                if sum(mgdict[s_molgroup]['areaaxis']):
                    f_com, f_int = average(mgdict[s_molgroup]['zaxis'], weights=mgdict[s_molgroup]['areaaxis'],
                                           returned=True)  # Calculate Center of Mass and Integral
                    f_int = sum(mgdict[s_molgroup]['areaaxis']) * (
                            mgdict[s_molgroup]['zaxis'][1] - mgdict[s_molgroup]['zaxis'][0])
                else:
                    f_com = 1E5
                    f_int = 0  # total integral in volume per surface area (unit: Å)
                    # taking into account grid spacing of z-axis
                f_int = f_int / mgdict['normarea']['areaaxis'][0]

                #                for j in range(len(mgdict[sMolgroup]['areaaxis'])):
                #                    if mgdict[sMolgroup]['areaaxis'][j]:
                #                        fCOM=fCOM-mgdict[sMolgroup]['zaxis'][j]
                #                        # normalize COM to start of molecular group
                #                        break

                f_avg = average(mgdict[s_molgroup]['areaaxis'])

                diResults[s_molgroup + '_COM'].append(f_com)
                diResults[s_molgroup + '_INT'].append(f_int)
                diResults[s_molgroup + '_AVG'].append(f_avg)

            # calculate ratios between frac1 and frac2
            if ('frac1' in l_molgroups) and ('frac2' in l_molgroups):
                if 'ratio_f1f2' not in list(diResults.keys()):
                    diResults['ratio_f1f2'] = []
                if 'ratio_f2f1' not in list(diResults.keys()):
                    diResults['ratio_f2f1'] = []
                diResults['ratio_f1f2'].append(diResults['frac1_INT'][-1] / diResults['frac2_INT'][-1])
                diResults['ratio_f2f1'].append(diResults['frac2_INT'][-1] / diResults['frac1_INT'][-1])

            # percentage water in sub-membrane space and other regions for tBLM and other bilayer
            # get vf_bilayer from molecular group lipid2
            vf_bilayer = float(mgdict['lipid2']['headerdata']['nf'])

            # prepare arrays for summing up molgroups
            total_components = numpy.zeros(len(mgdict['normarea']['areaaxis']))

            if 'headgroup1' in l_molgroups:
                total_components += numpy.array(mgdict['headgroup1']['areaaxis'])
                f_vol_headgroup1 = float(mgdict['headgroup1']['headerdata']['vol']) * \
                    float(mgdict['headgroup1']['headerdata']['nf'])
                # look for other headgroups
                j = 2
                while True:
                    if 'headgroup1_' + str(j) in l_molgroups:
                        total_components += numpy.array(mgdict['headgroup1_' + str(j)]['areaaxis'])
                        f_vol_headgroup1 += float(mgdict['headgroup1_' + str(j)]['headerdata']['vol']) * \
                            float(mgdict['headgroup1_' + str(j)]['headerdata']['nf'])
                    else:
                        break
                    j += 1

                if 'tether' in l_molgroups and 'tetherg' in l_molgroups \
                        and 'normarea' in l_molgroups and 'bME' in l_molgroups:

                    f_vol_submembrane = mgdict['normarea']['areaaxis'][0] * \
                        (float(mgdict['tether']['headerdata']['l']) + float(mgdict['tetherg']['headerdata']['l'])) * \
                        vf_bilayer

                    # sub membrane components
                    f_vol_components = float(mgdict['bME']['headerdata']['vol']) * \
                        float(mgdict['bME']['headerdata']['nf']) + float(mgdict['tether']['headerdata']['vol']) * \
                        float(mgdict['tether']['headerdata']['nf']) + float(mgdict['tetherg']['headerdata']['vol']) * \
                        float(mgdict['tetherg']['headerdata']['nf']) + f_vol_headgroup1

                    f_total_tether_length = float(mgdict['tether']['headerdata']['l']) + \
                        float(mgdict['tetherg']['headerdata']['l'])
                    if 'fTotalTetherLength' not in list(diResults.keys()):
                        diResults['fTotalTetherLength'] = []
                    diResults['fTotalTetherLength'].append(f_total_tether_length)

                    f_tether_density = float(mgdict['tether']['headerdata']['nf']) / mgdict['normarea']['areaaxis'][0]
                    if 'fTetherDensity' not in list(diResults.keys()):
                        diResults['fTetherDensity'] = []
                    diResults['fTetherDensity'].append(f_tether_density)

                    total_components += numpy.array(mgdict['bME']['areaaxis']) + \
                        numpy.array(mgdict['tether']['areaaxis']) + numpy.array(mgdict['tetherg']['areaaxis'])

            if 'lipid1' in l_molgroups and 'methyl1' in l_molgroups:

                f_total_lipid1_length = float(mgdict['lipid1']['headerdata']['l']) + float(
                    mgdict['methyl1']['headerdata']['l'])
                if 'fTotalLipid1Length' not in list(diResults.keys()):
                    diResults['fTotalLipid1Length'] = []
                diResults['fTotalLipid1Length'].append(f_total_lipid1_length)

                total_components += numpy.array(mgdict['lipid1']['areaaxis']) + \
                    numpy.array(mgdict['methyl1']['areaaxis'])

            if 'lipid2' in l_molgroups and 'methyl2' in l_molgroups:

                f_total_lipid2_length = float(mgdict['lipid2']['headerdata']['l']) + float(
                    mgdict['methyl2']['headerdata']['l'])
                if 'fTotalLipid2Length' not in list(diResults.keys()):
                    diResults['fTotalLipid2Length'] = []
                diResults['fTotalLipid2Length'].append(f_total_lipid2_length)

                f_area_per_lipid2 = float(mgdict['lipid2']['headerdata']['vol']) / float(
                    mgdict['lipid2']['headerdata']['l'])
                if 'fAreaPerLipid2' not in list(diResults.keys()):
                    diResults['fAreaPerLipid2'] = []
                diResults['fAreaPerLipid2'].append(f_area_per_lipid2)

                total_components += numpy.array(mgdict['lipid2']['areaaxis']) + \
                    numpy.array(mgdict['methyl2']['areaaxis'])

            if 'headgroup2' in l_molgroups:
                f_vol_headgroup2 = float(mgdict['headgroup2']['headerdata']['vol']) * \
                    float(mgdict['headgroup2']['headerdata']['nf'])
                total_components += numpy.array(mgdict['headgroup2']['areaaxis'])
                # look for other headgroups
                j = 2
                while True:
                    if 'headgroup2_' + str(j) in l_molgroups:
                        total_components += numpy.array(mgdict['headgroup2_' + str(j)]['areaaxis'])
                        f_vol_headgroup2 += float(mgdict['headgroup2_' + str(j)]['headerdata']['vol']) * \
                                          float(mgdict['headgroup2_' + str(j)]['headerdata']['nf'])
                    else:
                        break
                    j += 1

            if 'defect_hc' in l_molgroups:
                total_components = total_components + mgdict['defect_hc']['areaaxis']
                
            if 'defect_hg' in l_molgroups:
                total_components = total_components + mgdict['defect_hg']['areaaxis']

            if 'protein' in l_molgroups:
                total_components = total_components + mgdict['protein']['areaaxis']
                for i in range(len(total_components)):
                    if total_components[i] > vf_bilayer * mgdict['normarea']['areaaxis'][i]:
                        total_components[i] = vf_bilayer * mgdict['normarea']['areaaxis'][i]

            # calculate water and protein fractions if bilayer is present
            if 'headgroup1' in l_molgroups:
                f_start_hg1 = float(mgdict['headgroup1']['headerdata']['z']) - 0.5 * float(
                    mgdict['headgroup1']['headerdata']['l'])
                f_start_hc = float(mgdict['lipid1']['headerdata']['z']) - 0.5 * float(
                    mgdict['lipid1']['headerdata']['l'])
                f_start_methyl2 = float(mgdict['methyl2']['headerdata']['z']) - 0.5 * float(
                    mgdict['methyl2']['headerdata']['l'])
                f_start_hg2 = float(mgdict['headgroup2']['headerdata']['z']) - 0.5 * float(
                    mgdict['headgroup2']['headerdata']['l'])
                f_startbulk = float(mgdict['headgroup2']['headerdata']['z']) + 0.5 * float(
                    mgdict['headgroup2']['headerdata']['l'])

                f_step_size = mgdict['normarea']['zaxis'][1] - mgdict['normarea']['zaxis'][0]

                i_start_hg1 = int(floor(f_start_hg1 / f_step_size + 0.5))
                i_start_hc = int(floor(f_start_hc / f_step_size + 0.5))
                i_start_methyl2 = int(floor(f_start_methyl2 / f_step_size + 0.5))
                i_start_hg2 = int(floor(f_start_hg2 / f_step_size + 0.5))
                i_start_bulk = int(floor(f_startbulk / f_step_size + 0.5))

                ref = mgdict['normarea']['areaaxis']
                ratio = 1 - sum(total_components[0:i_start_hg1]) / sum(ref[0:i_start_hg1])
                if 'WaterFracSubMembrane' not in list(diResults.keys()):
                    diResults['WaterFracSubMembrane'] = []
                diResults['WaterFracSubMembrane'].append(ratio)
                ratio = 1 - sum(total_components[i_start_hg1:i_start_hc]) / sum(ref[i_start_hg1:i_start_hc])
                if 'WaterFracHeadgroup1' not in list(diResults.keys()):
                    diResults['WaterFracHeadgroup1'] = []
                diResults['WaterFracHeadgroup1'].append(ratio)
                ratio = 1 - sum(total_components[i_start_hc:i_start_methyl2]) / sum(ref[i_start_hc:i_start_methyl2])
                if 'WaterFracLipid1' not in list(diResults.keys()):
                    diResults['WaterFracLipid1'] = []
                diResults['WaterFracLipid1'].append(ratio)
                ratio = 1 - sum(total_components[i_start_methyl2:i_start_hg2]) / sum(ref[i_start_methyl2:i_start_hg2])
                if 'WaterFracLipid2' not in list(diResults.keys()):
                    diResults['WaterFracLipid2'] = []
                diResults['WaterFracLipid2'].append(ratio)
                ratio = 1 - sum(total_components[i_start_hc:i_start_hg2]) / sum(ref[i_start_hc:i_start_hg2])
                if 'WaterFracHydrocarbon' not in list(diResults.keys()):
                    diResults['WaterFracHydrocarbon'] = []
                diResults['WaterFracHydrocarbon'].append(ratio)
                ratio = 1 - sum(total_components[i_start_hg2:i_start_bulk]) / sum(ref[i_start_hg2:i_start_bulk])
                if 'WaterFracHeadgroup2' not in list(diResults.keys()):
                    diResults['WaterFracHeadgroup2'] = []
                diResults['WaterFracHeadgroup2'].append(ratio)

                # fraction of protein in certain parts of the membrane
                for group in ['protein', 'frac1', 'frac2']:
                    if group in l_molgroups:

                        if sum(mgdict[group]['areaaxis']) == 0.0:
                            f_frac_submembrane = 0
                            f_frac_inner_headgroup = 0
                            f_frac_inner_hydrocarbon = 0
                            f_frac_outer_hydrocarbon = 0
                            f_frac_hydrocarbon = 0
                            f_frac_outer_headgroup = 0
                            f_frac_headgroups = 0
                            f_frac_bulk = 0
                            f_frac_inner_leaflet = 0
                            f_frac_outer_leaflet = 0
                        else:
                            f_frac_submembrane = sum(mgdict[group]['areaaxis'][0:i_start_hg1]) / \
                                                 sum(mgdict[group]['areaaxis'])
                            f_frac_inner_headgroup = sum(mgdict[group]['areaaxis'][i_start_hg1:i_start_hc]) / \
                                                     sum(mgdict[group]['areaaxis'])
                            f_frac_inner_hydrocarbon = sum(mgdict[group]['areaaxis'][i_start_hc:i_start_methyl2]) / \
                                                       sum(mgdict[group]['areaaxis'])
                            f_frac_outer_hydrocarbon = sum(mgdict[group]['areaaxis'][i_start_methyl2:i_start_hg2]) / \
                                                       sum(mgdict[group]['areaaxis'])
                            f_frac_hydrocarbon = sum(mgdict[group]['areaaxis'][i_start_hc:i_start_hg2]) / \
                                                 sum(mgdict[group]['areaaxis'])
                            f_frac_outer_headgroup = sum(mgdict[group]['areaaxis'][i_start_hg2:i_start_bulk]) / \
                                                     sum(mgdict[group]['areaaxis'])
                            f_frac_headgroups = f_frac_inner_headgroup + f_frac_outer_headgroup
                            f_frac_bulk = sum(mgdict[group]['areaaxis'][i_start_bulk:]) / sum(mgdict[group]['areaaxis'])
                            f_frac_inner_leaflet = f_frac_inner_headgroup + f_frac_inner_hydrocarbon
                            f_frac_outer_leaflet = f_frac_outer_headgroup + f_frac_outer_hydrocarbon

                        if 'FracSubmembrane_' + group not in list(diResults.keys()):
                            diResults['FracSubmembrane_' + group] = []
                        diResults['FracSubmembrane_' + group].append(f_frac_submembrane)
                        if 'FracHydrocarbon_' + group not in list(diResults.keys()):
                            diResults['FracHydrocarbon_' + group] = []
                        diResults['FracHydrocarbon_' + group].append(f_frac_hydrocarbon)
                        if 'FracInnerHydrocarbon_' + group not in list(diResults.keys()):
                            diResults['FracInnerHydrocarbon_' + group] = []
                        diResults['FracInnerHydrocarbon_' + group].append(f_frac_inner_hydrocarbon)
                        if 'FracOuterHydrocarbon_' + group not in list(diResults.keys()):
                            diResults['FracOuterHydrocarbon_' + group] = []
                        diResults['FracOuterHydrocarbon_' + group].append(f_frac_outer_hydrocarbon)
                        if 'FracInnerHeadgroup_' + group not in list(diResults.keys()):
                            diResults['FracInnerHeadgroup_' + group] = []
                        diResults['FracInnerHeadgroup_' + group].append(f_frac_inner_headgroup)
                        if 'FracOuterHeadgroup_' + group not in list(diResults.keys()):
                            diResults['FracOuterHeadgroup_' + group] = []
                        diResults['FracOuterHeadgroup_' + group].append(f_frac_outer_headgroup)
                        if 'FracHeadgroups_' + group not in list(diResults.keys()):
                            diResults['FracHeadgroups_' + group] = []
                        diResults['FracHeadgroups_' + group].append(f_frac_headgroups)
                        if 'FracBulk_' + group not in list(diResults.keys()):
                            diResults['FracBulk_' + group] = []
                        diResults['FracBulk_' + group].append(f_frac_bulk)
                        if 'FracInnerLeaflet_' + group not in list(diResults.keys()):
                            diResults['FracInnerLeaflet_' + group] = []
                        diResults['FracInnerLeaflet_' + group].append(f_frac_inner_leaflet)
                        if 'FracOuterLeaflet_' + group not in list(diResults.keys()):
                            diResults['FracOuterLeaflet_' + group] = []
                        diResults['FracOuterLeaflet_' + group].append(f_frac_outer_leaflet)

                        # calculate peak position and FWHM for spline profile
                        imax, maxvalue, ifwhmminus, ifwhmplus = fnFindMaxFWHM(mgdict[group]['areaaxis'])
                        if 'PeakPosition_' + group not in list(diResults.keys()):
                            diResults['PeakPosition_' + group] = []
                        diResults['PeakPosition_' + group].append(mgdict[group]['zaxis'][imax] - f_startbulk)
                        if 'PeakValue_' + group not in list(diResults.keys()):
                            diResults['PeakValue_' + group] = []
                        diResults['PeakValue_' + group].append(mgdict[group]['areaaxis'][imax])
                        if 'FWHMMinusPosition_' + group not in list(diResults.keys()):
                            diResults['FWHMMinusPosition_' + group] = []
                        diResults['FWHMMinusPosition_' + group].append(mgdict[group]['zaxis'][ifwhmminus] - f_startbulk)
                        if 'FWHMPlusPosition_' + group not in list(diResults.keys()):
                            diResults['FWHMPlusPosition_' + group] = []
                        diResults['FWHMPlusPosition_' + group].append(mgdict[group]['zaxis'][ifwhmplus] - f_startbulk)
                        if 'FWHM_' + group not in list(diResults.keys()):
                            diResults['FWHM_' + group] = []
                        diResults['FWHM_' + group].append(
                            mgdict[group]['zaxis'][ifwhmplus] - mgdict[group]['zaxis'][ifwhmminus])

        if fConfidence > 1:
            fConfidence = 1
        if fConfidence < 0:
            fConfidence = special.erf(-1 * fConfidence / sqrt(2))

        fLowerPercentileMark = 100.0 * (1 - fConfidence) / 2
        fHigherPercentileMark = (100 - fLowerPercentileMark)

        File = open(self.mcmcpath + "/CalculationResults.dat", "w")
        for element, value in sorted(diResults.items()):
            fLowPerc = stats.scoreatpercentile(diResults[element], fLowerPercentileMark)  # Calculate Percentiles
            fMedian = stats.scoreatpercentile(diResults[element], 50.)
            fHighPerc = stats.scoreatpercentile(diResults[element], fHigherPercentileMark)

            sPrintString = '%(el)s'
            sPrintString += '  [%(lp)10.4g, %(m)10.4g, %(hp)10.4g] (-%(ld)10.4g, +%(hd)10.4g)'

            soutput = sPrintString % {'el': element, 'lp': fLowPerc, 'ld': (fMedian - fLowPerc), 'm': fMedian,
                                      'hd': (fHighPerc - fMedian), 'hp': fHighPerc}
            File.write(soutput)
            print(soutput)
            File.write('\n')

        File.close()

    # -------------------------------------------------------------------------------
    # calculates most likely, 68% and 94% confidence limits
    def fnCalcConfidenceLimits(self, data, method=1):

        # what follows is a set of three routines, courtesy to P. Kienzle, calculating
        # the shortest confidence interval

        def credible_interval(x, ci=0.95):
            """
            Find the credible interval covering the portion *ci* of the data.
            Returns the minimum and maximum values of the interval.
            *x* are samples from the posterior distribution.
           *ci* is the portion in (0,1], and defaults to 0.95.
            This function is faster if the inputs are already sorted.
            If *ci* is a vector, return a vector of intervals.
            """
            x.sort()

            # Simple solution: ci*N is the number of points in the interval, so
            # find the width of every interval of that size and return the smallest.
            if numpy.isscalar(ci):
                return _find_interval(x, ci)
            else:
                return [_find_interval(x, i) for i in ci]

        def _find_interval(x, ci):
            """
            Find credible interval ci in sorted, unweighted x
            """
            size = int(ci * len(x))
            if size > len(x) - 0.5:
                return x[0], x[-1]
            else:
                width = list(numpy.array(x[size:]) - numpy.array(x[:(-1 * size)]))
                idx = numpy.argmin(width)
                return x[idx], x[idx + size]

        # ------main function starting here------

        # traditional method, taking percentiles of the entire distribution
        if method == 1:
            return [stats.scoreatpercentile(data, percentile) for percentile in [2.3, 15.9, 50, 84.1, 97.7]]
        # alternative method, starting from maximum of distribution
        elif method == 2:
            (histo, low_range, bin_size, _) = stats.histogram(data, numbins=int(len(data) / 5))

            # normalize histogram
            sumlist = sum(histo)
            for i in range(len(histo)):
                histo[i] = float(histo[i]) / sumlist

            maxindex = numpy.argmax(histo)
            print(maxindex, histo[maxindex])
            if histo[maxindex] == 1:
                return [data[0], data[0], data[0], data[0], data[0]]

            # calculate a smoother maximum value
            a = maxindex;
            c = maxindex
            if 0 < maxindex < len(histo) - 1:
                a = maxindex - 1
                c = maxindex + 1
            maxindexsmooth = a * histo[a] + maxindex * histo[maxindex] + c * histo[c]
            maxindexsmooth = maxindexsmooth / (histo[a] + histo[maxindex] + histo[c])
            maxvaluesmooth = low_range + (maxindexsmooth + 0.5) * bin_size

            a = maxindex;
            c = maxindex
            confidence = histo[maxindex]
            while confidence < 0.6827:
                if a > 0:
                    a -= 1
                    confidence += histo[a]
                if c < len(histo) - 1:
                    c += 1
                    confidence += histo[c]
            onesigmam = low_range + (a + 0.5) * bin_size
            onesigmap = low_range + (c + 0.5) * bin_size

            while confidence < 0.9545:
                if a > 0:
                    a -= 1
                    confidence += histo[a]
                if c < len(histo) - 1:
                    c += 1
                    confidence += histo[c]

            twosigmam = low_range + (a + 0.5) * bin_size
            twosigmap = low_range + (c + 0.5) * bin_size

            # print data
            # print histo, twosigmam, onesigmam, maxvaluesmooth, onesigmap, twosigmap, low_range, bin_size

            return [twosigmam, onesigmam, maxvaluesmooth, onesigmap, twosigmap]

        # shortest confidence interval method, NIST recommended
        else:
            twosigmam, twosigmap = credible_interval(data, 0.95)
            onesigmam, onesigmap = credible_interval(data, 0.68)
            reported = 0.5 * (onesigmam + onesigmap)
            return [twosigmam, onesigmam, reported, onesigmap, twosigmap]

    # -------------------------------------------------------------------------------

    def fnCheckFit(self):  # checks fit directory for certain conditions

        if path.isfile('isErr.dat'):  # check which type of MC output present
            self.sStatFileName = 'isErr.dat'
        elif path.isfile('isErr.dat'):
            self.sStatFileName = 'sErr.dat'
        else:
            self.sStatFileName = ''

        try:
            self.fnLoadStatData()  # Load data from file into list
            self.fnLoadParameters()  # Load Parameters for limits

            if sorted(list(self.diParameters.keys())) != sorted(list(self.diStatResults['Parameters'].keys())):
                print('----------------------------------------------------')
                print('----------------------------------------------------')
                print('setup.c and Stat File do not agree -----------------')
                print('backing up Stat File -------------------------------')
                print('----------------------------------------------------')
                print('----------------------------------------------------')

                i = 2
                while True:
                    if path.isfile(self.sStatFileName[:-4] + str(i)[1:] + ".bak"):
                        pass
                    else:
                        pr = Popen(["mv ", self.sStatFileName, self.sStatFileName[:-4] + str(i)[1:] + ".bak"])
                        pr.wait()
                        self.sStatFileName = ''
                        self.liStatResult = []
                        self.sStatResultHeader = ''
                        break

                    i += 1

        except IOError:
            pass

        # -------------------------------------------------------------------------------

    def fnContourData(self, sname, dZGrid, dRhoGrid, dAreaGrid=0.1):

        def fnMap2Array(liArray, liDimensions, dx1, dy1, dx2, dy2):  # map nSLD profile onto array
            # by drawing lines between (dz1,drho1)
            # and (dz2,drho2)

            if not liArray:  # if array is empty create with first point
                liArray.append([0])
                liDimensions[0] = dx1
                liDimensions[1] = dx1
                liDimensions[3] = dy1
                liDimensions[4] = dy1

            while dx1 > liDimensions[1]:  # check if array has to be extended
                for i in range(len(liArray)):  # and extend it if necessary
                    liArray[i].append(0)
                liDimensions[1] = liDimensions[1] + liDimensions[2]
            while dx1 < liDimensions[0]:
                for i in range(len(liArray)):
                    liArray[i].insert(0, 0)
                liDimensions[0] = liDimensions[0] - liDimensions[2]
            while dy1 > liDimensions[4]:
                li = [0] * len(liArray[0])
                liArray.append(li)
                liDimensions[4] = liDimensions[4] + liDimensions[5]
            while dy1 < liDimensions[3]:
                li = [0] * len(liArray[0])
                liArray.insert(0, li)
                liDimensions[3] = liDimensions[3] - liDimensions[5]
            while dx2 > liDimensions[1]:
                for i in range(len(liArray)):
                    liArray[i].append(0)
                liDimensions[1] = liDimensions[1] + liDimensions[2]
            while dx2 < liDimensions[0]:
                for i in range(len(liArray)):
                    liArray[i].insert(0, 0)
                liDimensions[0] = liDimensions[0] - liDimensions[2]
            while dy2 > liDimensions[4]:
                li = [0] * len(liArray[0])
                liArray.append(li)
                liDimensions[4] = liDimensions[4] + liDimensions[5]
            while dy2 < liDimensions[3]:
                li = [0] * len(liArray[0])
                liArray.insert(0, li)
                liDimensions[3] = liDimensions[3] - liDimensions[5]

            iXStart = int((dx1 - liDimensions[0]) / liDimensions[2])  # calculate coordinates from (z,rho)
            iXEnd = int((dx2 - liDimensions[0]) / liDimensions[2])  # data points
            iYStart = int((dy1 - liDimensions[3]) / liDimensions[5])
            iYEnd = int((dy2 - liDimensions[3]) / liDimensions[5])

            diffX = iXEnd - iXStart  # how many indizes do go?
            diffY = iYEnd - iYStart

            if abs(diffX) > abs(diffY):  # which direction more steps?
                diff = abs(diffX)
            else:
                diff = abs(diffY)

            if diff == 0:  # treat single point differently because
                liArray[iYStart][iXStart] = liArray[iYStart][iXStart]  # of division by zero error -> do nothing!
            else:
                fStepX = float(diffX) / float(diff)  # calculate stepsize for each direction
                fStepY = float(diffY) / float(diff)

                i = 0  # draw by increasing field occupancy
                while i < diff:  # and thus create surface/contour plot
                    iX = iXStart + int(round(fStepX * float(i), 0))
                    iY = iYStart + int(round(fStepY * float(i), 0))
                    liArray[iY][iX] += 1
                    i += 1

        self.fnRecreateStatistical()  # Load Parameters and profiles from stat data
        iNumberOfModels = self.fnGetNumberOfModelsFromSetupC()  # how many profiles to analyze?

        liContourArray = []  # initiate array
        liContourArrayDimensions = []  # array dimensions
        i = 0  # result for contour plot of profiles
        while i < iNumberOfModels:  # list of lists for each profile
            liContourArray = liContourArray + [[]]
            liContourArrayDimensions = liContourArrayDimensions + [[0., 0., dRhoGrid, 0., 0., dZGrid]]
            # dimension format ymin, ymax, ystep, xmin, xmax, xstep
            i = i + 1

        for iteration in self.diStatResults['nSLDProfiles']:  # cycle through all individual stat. results
            i = 0
            for model in iteration:  # cycle through all models
                for l in range(len(model[0])):  # extract nSLD profile data point by point
                    dz = round(model[0][l] / dZGrid) * dZGrid  # round to grid precision
                    drho = round(model[1][l] / dRhoGrid) * dRhoGrid  # round nSLD to grid precision
                    if l != 0:
                        fnMap2Array(liContourArray[i],
                                    liContourArrayDimensions[i],
                                    drhoold, dzold, drho, dz)

                    dzold = dz
                    drhoold = drho
                i += 1

        i = 0
        print('Processing data for output ...')
        while i < iNumberOfModels:  # loop through all models

            print('Model %i: %i x %i' % (i, len(liContourArray[i][0]),
                                         len(liContourArray[i])))
            sFileName = 'Cont_nSLD_Array' + str(i) + '.dat'  # write out array
            file = open(sFileName, "w")
            for line in liContourArray[i]:
                sLine = ''
                for element in line:
                    sLine = sLine + str(element) + ' '
                sLine = sLine + '\n'
                file.write(sLine)
            file.close()

            dZMin = liContourArrayDimensions[i][3]
            dZMax = liContourArrayDimensions[i][4]
            dZStep = liContourArrayDimensions[i][5]
            dRhoMin = liContourArrayDimensions[i][0]
            dRhoMax = liContourArrayDimensions[i][1]
            dRhoStep = liContourArrayDimensions[i][2]

            dZ = dZMin  # write out x-dimension wave
            sFileName = 'Cont_nSLD_DimZ' + str(i) + '.dat'  # dimension wave has one point extra for Igor
            file = open(sFileName, "w")
            while (dZ <= dZMax + dZStep):
                sLine = str(round(dZ / dZGrid) * dZGrid) + '\n'
                file.write(sLine)
                dZ = dZ + dZStep
            file.close()

            dRho = dRhoMin  # write out y-dimension wave
            sFileName = 'Cont_nSLD_DimRho' + str(i) + '.dat'  # dimension wave has one point extra for Igor
            file = open(sFileName, "w")
            while (dRho <= dRhoMax + dRhoStep):
                sLine = str(round(dRho / dRhoGrid) * dRhoGrid) + '\n'
                file.write(sLine)
                dRho = dRho + dRhoStep
            file.close()

            i = i + 1

        if 'Molgroups' in self.diStatResults:

            liContourArray = []  # initiate array
            liContourArrayDimensions = []  # array dimensions
            for _ in self.diStatResults['Molgroups'][0]:  # iterate through molecular group names
                liContourArray = liContourArray + [[]]
                liContourArrayDimensions = liContourArrayDimensions + [[0., 0., dAreaGrid, 0., 0., dZGrid]]
                # dimension format ymin, ymax, ystep, xmin, xmax, xstep

            for iteration in self.diStatResults['Molgroups']:  # cycle through all individual stat. results
                i = 0;
                dareaold = 0;
                dzold = 0
                for molgroup in iteration:  # cycle through all molecular groups
                    for l in range(len(iteration[molgroup]['zaxis'])):  # extract area profile data point by point
                        dz = round(iteration[molgroup]['zaxis'][l] / dZGrid) * dZGrid  # round to grid precision
                        darea = round(
                            iteration[molgroup]['areaaxis'][l] / dAreaGrid) * dAreaGrid  # round nSLD to grid precision
                        if l != 0:
                            fnMap2Array(liContourArray[i],
                                        liContourArrayDimensions[i],
                                        dareaold, dzold, darea, dz)

                        dzold = dz
                        dareaold = darea
                    i += 1

            i = 0
            for molgroup in self.diStatResults['Molgroups'][0]:  # loop through all models

                print('%s %i: %i x %i' % (molgroup, i, len(liContourArray[i][0]),
                                          len(liContourArray[i])))
                sFileName = 'Cont_' + molgroup + '_Array' + '.dat'  # write out array
                file = open(sFileName, "w")
                for line in liContourArray[i]:
                    sLine = ''
                    for element in line:
                        sLine = sLine + str(element) + ' '
                    sLine = sLine + '\n'
                    file.write(sLine)
                file.close()

                dZMin = liContourArrayDimensions[i][3]
                dZMax = liContourArrayDimensions[i][4]
                dZStep = liContourArrayDimensions[i][5]
                dAreaMin = liContourArrayDimensions[i][0]
                dAreaMax = liContourArrayDimensions[i][1]
                dAreaStep = liContourArrayDimensions[i][2]

                dZ = dZMin  # write out x-dimension wave
                sFileName = 'Cont_' + molgroup + '_DimZ' + '.dat'  # dimension wave has one point extra for Igor
                file = open(sFileName, "w")
                while (dZ <= dZMax + dZStep):
                    sLine = str(round(dZ / dZGrid) * dZGrid) + '\n'
                    file.write(sLine)
                    dZ = dZ + dZStep
                file.close()

                dArea = dAreaMin  # write out y-dimension wave
                sFileName = 'Cont_' + molgroup + '_DimArea' + '.dat'  # dimension wave has one point extra for Igor
                file = open(sFileName, "w")
                while (dArea <= dAreaMax + dAreaStep):
                    sLine = str(round(dArea / dAreaGrid) * dAreaGrid) + '\n'
                    file.write(sLine)
                    dArea = dArea + dAreaStep
                file.close()

                i = i + 1

            # -------------------------------------------------------------------------------

    def fnGetChiSq(self):  # export chi squared
        return self.chisq

    # -------------------------------------------------------------------------------

    def fnCreateBilayerPlotData(self):

        # integrate over 1D array
        def fnIntegrate(axis, array, start, stop):

            startaxis = float(axis[0])
            incaxis = float(axis[1] - axis[0])
            startindex = int((start - startaxis) / incaxis + 0.5)
            stopindex = int((stop - startaxis) / incaxis + 0.5)

            if stopindex > startindex:
                sum = 0.5 * (array[startindex] + array[stopindex])
            else:
                sum = array[startindex]

            for i in range(startindex + 1, stopindex):
                sum += array[i]

            return sum

        # find maximum values and indizees of half-height points assuming unique solution and steady functions
        def fnMaximumHalfPoint(data):

            max = numpy.amax(data)

            point1 = False
            point2 = False
            hm1 = 0
            hm2 = 0

            for i in range(len(data)):
                if data[i] > (max / 2) and not point1:
                    point1 = True
                    hm1 = i
                if data[i] < (max / 2) and not point2:
                    point2 = True
                    hm2 = i - 1
                    i = len(data)

            return max, hm1, hm2

        def fnStat(diarea, name, diStat):
            for i in range(len(diarea[list(diarea.keys())[0]])):
                liOnePosition = [iteration[i] for key, iteration in diarea.items() if key != 'zaxis']
                stat = self.fnCalcConfidenceLimits(liOnePosition, method=1)
                diStat[name + '_msigma'].append(stat[1])
                diStat[name].append(stat[2])
                diStat[name + '_psigma'].append(stat[3])

        # initialize Statistical Dictionary
        print('Initializing ...')
        lGroupList = ['substrate', 'siox', 'tether', 'innerhg', 'innerhc', 'outerhc', 'outerhg', 'protein',
                      'sum', 'water']
        diStat = {}
        for element in lGroupList:
            diStat[element] = []
            diStat[element + '_corr'] = []
            diStat[element + '_cvo'] = []
            diStat[element + '_corr_cvo'] = []

        keylist = list(diStat)
        for element in keylist:
            diStat[element + '_msigma'] = []
            diStat[element + '_psigma'] = []

        diIterations = {}
        for element in lGroupList:
            diIterations[element] = {}
            diIterations[element + '_corr'] = {}
            diIterations[element + '_cvo'] = {}
            diIterations[element + '_corr_cvo'] = {}

        # pull all relevant molgroups
        # headgroups are allready corrected for protein penetration (_corr) and will be copied over to the _corr
        # entries next
        print('Pulling all molgroups ...')
        print('  substrate ...')
        diIterations['substrate'], __, __ = self.fnPullMolgroupLoader(['substrate'])
        print('  siox ...')
        diIterations['siox'], __, __ = self.fnPullMolgroupLoader(['siox'])
        print('  tether ...')
        diIterations['tether'], __, __ = self.fnPullMolgroupLoader(['bME', 'tetherg', 'tether'])
        print('  innerhg ...')
        diIterations['innerhg'], __, __ = self.fnPullMolgroupLoader(['headgroup1', 'headgroup1_2', 'headgroup1_3'])
        print('  innerhc ...')
        diIterations['innerhc'], __, __ = self.fnPullMolgroupLoader(['lipid1', 'methyl1'])
        print('  outerhc ...')
        diIterations['outerhc'], __, __ = self.fnPullMolgroupLoader(['lipid2', 'methyl2'])
        print('  outerhg ...')
        diIterations['outerhg'], __, __ = self.fnPullMolgroupLoader(['headgroup2', 'headgroup2_2', 'headgroup2_3'])
        print('  protein ...')
        diIterations['protein'], __, __ = self.fnPullMolgroupLoader(['protein'])

        # save z-axis
        diStat['zaxis'] = numpy.copy(diIterations['substrate']['zaxis'])

        # shallow copies of the uncorrected data into the corrected dictionaries
        # and the values will be replaced by their modifications step by step
        for element in lGroupList:
            diIterations[element + '_corr'] = diIterations[element].copy()

        # loop over all iterations and apply the corrections / calculations
        print('Applying corrections ...\n')
        for i in range(len(list(diIterations['substrate'].keys())) - 1):

            key = 'iter' + str(i)
            substrate = numpy.array(diIterations['substrate'][key])
            siox = numpy.array(diIterations['siox'][key])
            tether = numpy.array(diIterations['tether'][key])
            innerhg_corr = numpy.array(diIterations['innerhg_corr'][key])
            innerhc = numpy.array(diIterations['innerhc'][key])
            outerhc = numpy.array(diIterations['outerhc'][key])
            outerhg_corr = numpy.array(diIterations['outerhg_corr'][key])
            protein = numpy.array(diIterations['protein'][key])
            axis = numpy.array(diIterations['substrate']['zaxis'])

            hc = innerhc + outerhc
            # this is the sum as used for the joining procedure
            sum = substrate + siox + tether + innerhg_corr + hc + outerhg_corr

            areaperlipid, _, _ = fnMaximumHalfPoint(substrate)
            maxbilayerarea, _, _ = fnMaximumHalfPoint(hc)
            # vf_bilayer = maxbilayerarea/areaperlipid

            # recuperate the non-corrected headgroup distributions that were not saved to file by the fit
            # by reversing the multiplication based on the amount of replaced hc material
            __, hc1_hm1, hc1_hm2 = fnMaximumHalfPoint(innerhc)
            __, hc2_hm1, hc2_hm2 = fnMaximumHalfPoint(outerhc)
            hg1ratio = fnIntegrate(axis, protein, axis[hc1_hm1], axis[hc1_hm2]) \
                       / fnIntegrate(axis, innerhc, axis[hc1_hm1], axis[hc1_hm2])
            hg2ratio = fnIntegrate(axis, protein, axis[hc2_hm1], axis[hc2_hm2]) \
                       / fnIntegrate(axis, outerhc, axis[hc2_hm1], axis[hc2_hm2])
            innerhg = innerhg_corr / (1 - hg1ratio)
            diIterations['innerhg'][key] = numpy.copy(innerhg)
            outerhg = outerhg_corr / (1 - hg2ratio)
            diIterations['outerhg'][key] = numpy.copy(outerhg)

            # prepare arrays for correction
            innerhc_corr = numpy.copy(innerhc)
            outerhc_corr = numpy.copy(outerhc)

            # correct the hc profiles due to protein penetration
            for i in range(len(protein)):
                if sum[i] + protein[i] > maxbilayerarea:
                    if (innerhc[i] + outerhc[i]) > 0:
                        excess = sum[i] + protein[i] - maxbilayerarea
                        if excess > (innerhc[i] + outerhc[i]):
                            excess = (innerhc[i] + outerhc[i])
                        # print (innerhc[i]+outerhc[i]) > 0, i
                        # print 'first' , innerhc_corr[i], excess, innerhc[i], outerhc[i],
                        # excess*innerhc[i]/(innerhc[i]+outerhc[i])
                        innerhc_corr[i] -= excess * innerhc[i] / (innerhc[i] + outerhc[i])
                        # print 'second' , outerhc_corr[i], excess, innerhc[i], outerhc[i],
                        # excess*outerhc[i]/(innerhc[i]+outerhc[i])
                        outerhc_corr[i] -= excess * outerhc[i] / (innerhc[i] + outerhc[i])

            # update dictionary entries for later statistics
            diIterations['innerhc_corr'][key] = numpy.copy(innerhc_corr)
            diIterations['outerhc_corr'][key] = numpy.copy(outerhc_corr)
            sum_corr = substrate + siox + tether + innerhg_corr + innerhc_corr + outerhc_corr + outerhg_corr
            diIterations['sum_corr'][key] = numpy.copy(sum_corr)

            # this is the truly non-corrected sum, different from the previous sum used for joining
            sum = substrate + siox + tether + innerhg + innerhc + outerhc + outerhg
            diIterations['sum'][key] = numpy.copy(sum)
            water_corr = areaperlipid - sum_corr - protein
            diIterations['water_corr'][key] = numpy.copy(water_corr)
            water = areaperlipid - sum - protein
            diIterations['water'][key] = numpy.copy(water)

            # calculate volume occupancy distributions by division by the area per lipid
            for element in lGroupList:
                diIterations[element + '_cvo'][key] = diIterations[element][key] / areaperlipid
                diIterations[element + '_corr_cvo'][key] = diIterations[element + '_corr'][key] / areaperlipid

        # calculate the statisics
        print('Calculating statistics ...\n')
        for element in lGroupList:
            if element != 'zaxis':
                fnStat(diIterations[element], element, diStat)
                fnStat(diIterations[element + '_corr'], element + '_corr', diStat)
                fnStat(diIterations[element + '_cvo'], element + '_cvo', diStat)
                fnStat(diIterations[element + '_corr_cvo'], element + '_corr_cvo', diStat)

        print('Saving data to bilayerplotdata.dat ...\n')
        self.Interactor.fnSaveSingleColumns(self.mcmcpath + '/bilayerplotdata.dat', diStat)

    def fnGetNumberOfModelsFromSetupC(self):

        file = open(self.setupfilename, "r")  # open setup.c
        data = file.readlines()
        file.close()
        smatch = compile(r'define\s+MODELS\s+(.+?)\n', IGNORECASE | VERBOSE)
        for line in data:  # search through setup.c
            if smatch.search(line):  # searching for MODELS constant
                i = smatch.search(line).group(1)
                return (int(i))
        return (0)

    def fnGetParameterValue(self, sname):  # export absolute parameter value
        return self.diParameters[sname]['value']  # for given name

    # -------------------------------------------------------------------------------

    def fnGetSortedParNames(self):  # return a list of sorted parameter
        litest = list(self.diParameters.keys())
        litest = sorted(litest, key=lambda keyitem: self.diParameters[keyitem]['number'])
        return litest

    # -------------------------------------------------------------------------------

    def fnGetSortedParValues(self):  # the same as above but it returns
        litest = list(self.diParameters.keys())  # a list of sorted parameter values
        litest = sorted(litest, key=lambda keyitem: self.diParameters[keyitem]['number'])
        lvalue = []
        for parameter in litest:
            lvalue.append(str(self.diParameters[parameter]['value']))
        return lvalue

    # -------------------------------------------------------------------------------

    def fnGetTaggedParameters(self):  # returns a list of the name and the
        file = open(self.setupfilename)  # range + stepsize information of parameters
        data = file.readlines()  # which are tagged for displacement error
        file.close()  # analysis
        output = []
        for line in data:
            if '!rstag' in line:
                smatch = compile(r'pars_add\(pars.*?\"(.+?)\"\s*,.+?,(.+?),(.+?)\).+!rstag\s+(.+?)\s+(.+?)\s+(.+?)\s+!',
                                 IGNORECASE | VERBOSE)
                output.append(smatch.search(line).groups())
        return output

    def fnImportMCMCBestFit(self, sPath):  # imports best-fit from MCMC

        call(['cp', 'setup.c', 'setup.back'])
        call(['cp', 'setup.cc', 'setup.backcc'])

        self.fnLoadParameters(sPath)  # Load Parameters and modify setup.cc

        # get list of parameters from setup.c/par.dat and sort by number of appearance in setup.cc
        li_parameters = list(self.diParameters.keys())
        li_parameters = sorted(li_parameters, key=lambda keyitem: self.diParameters[keyitem]['number'])

        # change setup.c to quasi fix all parameters
        li_addition = []
        for parameter in li_parameters:
            li_addition.append(('%s = %s;\n' %
                                (self.diParameters[parameter]['variable'], self.diParameters[parameter]['value'])))
        self.fnWriteConstraint2SetupC(li_addition)
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

    def fnLoadAndPrintPar(self, sPath='./'):
        self.fnLoadParameters(sPath)
        self.fnLoadCovar(sPath)
        self.fnPrintPar()

    def fnLoadCovar(self, sPath):
        # loads a covariance matrix,
        # calculates the errors and stores
        # it into self.diParameters
        # the function fnLoadParameters
        # must have been already carried out

        if path.isfile(sPath + 'covar.dat'):
            File = open(sPath + 'covar.dat')
            data = File.readlines()
            File.close()
            for i in range(len((data[1]).split())):  # number of columns in second row
                for parameter in list(self.diParameters.keys()):  # search for parameter number and
                    if self.diParameters[parameter]['number'] == i:  # retrieve value
                        fvalue = self.diParameters[parameter]['value']
                        ferror = float((data[i + 1]).split()[i])  # the first row contains a comment
                        if ferror < 0:
                            ferror = 0
                        ferror = sqrt(ferror) * fvalue  # calculate error
                        self.diParameters[parameter]['error'] = ferror
                        break
        else:
            for parameter in list(self.diParameters.keys()):
                self.diParameters[parameter]['error'] = float(0)  # fill all errors with zero

    def fnLoadFileListAndChangeToLocal(self):
        # scans the setup.c file and creates a list of the filenames of the reflectivity data files, it also
        # modifies setup.c in this way that it loads copies of reflectivity dat files located in the working directory
        # with the file ending .mce, those files are the modified files by the statistical error analysis

        File = open(self.setupfilename, "r")
        data = File.readlines()
        File.close()
        newdata = []
        filelist = []
        smatch1 = compile(r'fit_data.+?\"(.+?)\"',
                          IGNORECASE | VERBOSE)
        smatch2 = compile(r'(fit_data.+?\").+?(\")',
                          IGNORECASE | VERBOSE)
        smatch3 = compile(r'\s*//', IGNORECASE | VERBOSE)
        for line in data:
            if ('fit_data' in line) and (not smatch3.match(line)):  # scan for file loading but not comment lines
                filelist.append(smatch1.search(line).group(1))  # append path+filename to filelist
                newdata.append(smatch2.sub(r'\1' +
                                           path.split(filelist[-1])[-1] + r'.mce\2',
                                           line))  # modifiy setup.c line with new filename
            else:
                newdata.append(line)  # no file loading -> do not modify setup.c line
        file = open(self.setupfilename, "w")  # write setup.c
        file.writelines(newdata)
        file.close()
        return filelist

    @staticmethod
    def fnLoadObject(sFileName):
        import pickle
        File = open(sFileName, "rb")
        Object = pickle.load(File)
        File.close()
        return Object

    def fnLoadParameters(self):
        # loads the last row's parameters, and
        # ranges from par.dat and stores them into
        # self.diParameters; parameter names are
        # read from setup.c, since par.dat truncates
        # parameter names over 17 characters

        # after reading in the parameters check for
        # definitions of GAUSSIAN,
        # and check vs. active fit parameters
        # load them into self.molgroups dictionary
        self.diParameters = {}
        self.diParameters, self.chisq = self.Interactor.fnLoadParameters()

    def fnLoadStatData(self, sparse=0):
        self.diStatResults = self.Interactor.fnLoadStatData(sparse)
        # cycle through all parameters
        # determine length of longest parameter name for displaying
        iMaxParameterNameLength = 0
        for parname in list(self.diStatResults['Parameters'].keys()):
            if len(parname) > iMaxParameterNameLength:
                iMaxParameterNameLength = len(parname)
        self.diStatResults['MaxParameterLength'] = iMaxParameterNameLength

        self.diStatResults['NumberOfStatValues'] = \
            len(self.diStatResults['Parameters'][list(self.diStatResults['Parameters'].keys())[0]]['Values'])
        #print('StatDataLength: ')
        #print(len(self.diStatResults['Parameters'][list(self.diStatResults['Parameters'].keys())[0]]['Values']))

    def fnnSLDEnvelopes(self, fGrid, fSigma, sname, shortflag=False, iContrast=-1):

        def fnInterpolateData(xdata, ydata, fMin, fMax, fGrid):

            f = fMin  # target start
            ix = 0  # data start
            liInterpolated = [[], []]  # interpolation result

            while f <= fMax:  # main loop
                # print f,fMin,fMax,fGrid,ix,xdata[ix],len(xdata)
                if f < xdata[0]:  # fill area where no data available with first value
                    liInterpolated[0].append(f)
                    liInterpolated[1].append(ydata[0])
                    f = f + fGrid
                elif f > xdata[-1]:  # fill remaining cells with last value
                    liInterpolated[0].append(f)
                    liInterpolated[1].append(ydata[-1])
                    f = f + fGrid
                else:  # at least one data point surpassed by f
                    while (f > xdata[ix]) and (ix < (len(xdata) - 1)):  # searching first data point past f
                        ix = ix + 1
                    if f < xdata[ix]:  # was there a data point past f?
                        LowerX = ix - 1  # calculate data surrounding f
                        UpperX = ix
                        fDataGrid = xdata[UpperX] - xdata[LowerX]
                        fInterpolate = ydata[LowerX] * ((f - xdata[LowerX]) / fDataGrid  # do the interpolation
                                                        ) + ydata[UpperX] * ((xdata[UpperX] - f) / fDataGrid)
                        liInterpolated[0].append(f)
                        liInterpolated[1].append(fInterpolate)
                        f = f + fGrid
                    elif f == xdata[ix]:  # no interpolation needed
                        liInterpolated[0].append(xdata[ix])
                        liInterpolated[1].append(ydata[ix])
                        f = f + fGrid

            return liInterpolated

        def fnStoreEnvelope(liStoredEnvelopes, liStoredEnvelopeHeaders, envelope, percentile):
            iposition = len(liStoredEnvelopes) / 2
            liStoredEnvelopes.insert(iposition, envelope[1])
            liStoredEnvelopes.insert(iposition, envelope[0])
            liStoredEnvelopeHeaders.insert(iposition, str(1 - percentile / 2))
            liStoredEnvelopeHeaders.insert(iposition, str(percentile / 2))

        def fnSaveEnvelopes(liStoredEnvelopes, liStoredEnvelopeHeaders, iModel, fMin, fMax, fGrid, fSigma):
            file = open('Envelopes' + str(iModel) + '.dat', "w")

            if fSigma == 0:  # save all computed envelopes

                liSaveEnvelopeHeaders = liStoredEnvelopeHeaders
                liSaveEnvelopes = liStoredEnvelopes


            else:  # save only multiples of fSigma in Sigma

                liSaveEnvelopeHeaders = []
                liSaveEnvelopes = []

                fmult = 0.
                while 1:

                    fConfidence = special.erf(fmult / sqrt(2))  # upper and lower Percentiles for fmult*sigma
                    fLowerPerc = (1 - fConfidence) / 2
                    fUpperPerc = 1 - (1 - fConfidence) / 2

                    fStepSize = 1 / float(len(liStoredEnvelopeHeaders))  # positions in envelopes list
                    iLowerPerc = int(floor(fLowerPerc / fStepSize))
                    iUpperPerc = int(ceil(fUpperPerc / fStepSize))

                    if (iLowerPerc == 0) or (iUpperPerc >= len(liStoredEnvelopeHeaders) - 1):
                        break

                    liSaveEnvelopeHeaders.insert(0, 'minus' + str(fmult) + 'sigma')
                    liSaveEnvelopes.insert(0, liStoredEnvelopes[iLowerPerc])

                    if iUpperPerc != iLowerPerc:
                        liSaveEnvelopeHeaders.append('plus' + str(fmult) + 'sigma')
                        liSaveEnvelopes.append(liStoredEnvelopes[iUpperPerc])

                    fmult = fmult + fSigma

                liSaveEnvelopeHeaders.insert(0, 'LowerEnvelope')
                liSaveEnvelopes.insert(0, liStoredEnvelopes[0])
                liSaveEnvelopeHeaders.append('UpperEnvelope')
                liSaveEnvelopes.append(liStoredEnvelopes[len(liStoredEnvelopeHeaders) - 1])

            file.write("z ")
            for element in liSaveEnvelopeHeaders:
                if fSigma == 0:
                    file.write("p" + element + " ")
                else:
                    file.write(element + " ")
            file.write("\n")

            f = fMin
            for i in range(len(liSaveEnvelopes[0])):
                file.write(str(f) + " ")
                for element in liSaveEnvelopes:
                    file.write(str(element[i]) + " ")
                file.write("\n")
                f = f + fGrid

            file.close()

        def fnCalculateEnvelope(profilelist):

            LowerEnvelope = array(profilelist[0][1])
            UpperEnvelope = array(profilelist[0][1])

            for iteration in profilelist:
                LowerEnvelope = minimum(LowerEnvelope, array(iteration[1]))
                UpperEnvelope = maximum(UpperEnvelope, array(iteration[1]))

            fArea = sum(subtract(UpperEnvelope, LowerEnvelope))

            envelope = [LowerEnvelope.tolist(), UpperEnvelope.tolist()]

            return fArea, envelope

        self.fnRecreateStatistical()  # Load Parameters and profiles from stat data
        iNumberOfModels = self.fnGetNumberOfModelsFromSetupC()  # how many profiles to analyze?

        if iContrast == -1:
            iModelStart = 0
            iModelEnd = len(self.diStatResults['nSLDProfiles'][0])
        else:
            iModelStart = iContrast
            iModelEnd = iContrast + 1

        for iModel in range(iModelStart, iModelEnd):  # do the analysis separately for each model
            print('Analyzing model %i ...' % (iModel))

            liStoredEnvelopes = []  # List of stored envelopes
            liStoredEnvelopeHeaders = []  # and their headers

            fMin = 0
            fMax = 0  # initializing boundaries
            profilelist = []  # extracting all profiles related to the actual
            for iteration in self.diStatResults['nSLDProfiles']:  # model
                profilelist.append(iteration[iModel][:])
                if fMin > profilelist[-1][0][0]:
                    fMin = profilelist[-1][0][0]
                if fMax < profilelist[-1][0][-1]:
                    fMax = profilelist[-1][0][-1]
            fMax = floor((fMax - fMin) / fGrid) * fGrid + fMin  # make fMax compatible with fGrid and fMin

            print('Rebinning data...')
            for i in range(len(profilelist)):
                profilelist[i] = fnInterpolateData(profilelist[i][0], profilelist[i][1], fMin, fMax, fGrid)

            iNumberOfProfiles = len(profilelist)

            if not shortflag:
                for iPercentile in range(iNumberOfProfiles):

                    print('Calculating %f percentile...' % (1 - float(iPercentile) / float(iNumberOfProfiles)))

                    fArea, liEnvelope = fnCalculateEnvelope(profilelist)  # calculate envelope
                    fnStoreEnvelope(liStoredEnvelopes, liStoredEnvelopeHeaders,
                                    liEnvelope, float(iPercentile) / float(iNumberOfProfiles))  # store envlope in list

                    if iPercentile != iNumberOfProfiles - 1:
                        iMaxReduction = 0
                        for i in range(len(profilelist)):  # eliminate the profile with the largest reduction
                            # in area
                            profiletestlist = profilelist[:]
                            profiletestlist.pop(i)
                            fTestArea, fTestEnvelope = fnCalculateEnvelope(profiletestlist)
                            if fArea > fTestArea:
                                iMaxReduction = i
                                fArea = fTestArea

                        profilelist.pop(iMaxReduction)  # eliminate profile
            else:
                fSigmaStep = 0.1
                iPercentile = 0
                fArea, liEnvelope = fnCalculateEnvelope(profilelist)  # calculate envelope
                fnStoreEnvelope(liStoredEnvelopes, liStoredEnvelopeHeaders,
                                liEnvelope, float(iPercentile) / float(iNumberOfProfiles))  # store envlope in list
                iPercentile += 1

                while iPercentile < iNumberOfProfiles:

                    print('Short: Calculating %f percentile...' % (1 - float(iPercentile) / float(iNumberOfProfiles)))

                    lScoring = []  # list of profile scores
                    for i in range(len(profilelist)):  # eliminate the profile with the largest reduction
                        # in area
                        profiletestlist = profilelist[:]
                        profiletestlist.pop(i)
                        fTestArea, fTestEnvelope = fnCalculateEnvelope(profiletestlist)
                        lScoring.append([i, fTestArea])
                        # print lScoring

                    lScoring = sorted(lScoring, key=itemgetter(1))  # sort by lowest area

                    fConf = (1 - float(iPercentile) / float(iNumberOfProfiles))  # momentary confidence level
                    fSig = special.erfinv(fConf) * sqrt(2)  # related sigma value
                    fSigT = floor(fSig / fSigmaStep) * fSigmaStep  # get to the lower sigma step
                    fConfT = special.erf(fSigT / sqrt(2))  # target confidence level
                    # number of profiles to eliminate
                    iElimination = iNumberOfProfiles - iPercentile - int(fConfT * iNumberOfProfiles)

                    print("iPercentile %i iNumberOfProfiles %i iElimination %i fConf %e fSig %e fSigT %e fConfT %e" % (
                        iPercentile, iNumberOfProfiles, iElimination, fConf, fSig, fSigT, fConfT))

                    lScoring = sorted(lScoring[0:iElimination], key=itemgetter(0), reverse=True)

                    for i in range(iElimination):
                        profilelist.pop(lScoring[i][0])  # delete profiles starting with highest indices
                        fArea, liEnvelope = fnCalculateEnvelope(profilelist)  # calculate envelope
                        fnStoreEnvelope(liStoredEnvelopes, liStoredEnvelopeHeaders,
                                        liEnvelope,
                                        float(iPercentile) / float(iNumberOfProfiles))  # store envlope in list
                        iPercentile += 1

                        # raw_input("Press Enter to continue...")

            fnSaveEnvelopes(liStoredEnvelopes, liStoredEnvelopeHeaders, iModel,
                            fMin, fMax, fGrid, fSigma)

    def fnPlotMolgroups(self):

        from matplotlib.font_manager import fontManager, FontProperties

        font = FontProperties(size='x-small')
        plotname = 'Molgroups'

        self.fnLoadParameters()  # Load Parameters and modify setup.cc
        self.fnBackup()  # Backup setup.c, and other files
        try:
            liParameters = list(self.diParameters.keys())  # get list of parameters from setup.c/par.dat
            liParameters = sorted(liParameters, key=lambda keyitem: self.diParameters[keyitem][
                'number'])  # sort by number of appereance in setup.c

            liAddition = []
            for parameter in liParameters:  # cycle through all parameters
                liAddition.append(('%s = %s;\n' %  # change setup.c to quasi fix all parameters
                                   (self.diParameters[parameter]['variable'], self.diParameters[parameter]['value'])))
            self.fnWriteConstraint2SetupC(liAddition)  # write out
            call(["rm", "-f", "mol.dat"])
            self.fnMake()  # compile changed setup.c
            call(["./fit", "-o"])  # write out profile.dat and fit.dat
            call(["sync"])  # synchronize file system
            sleep(1)  # wait for system to clean up

        finally:
            self.fnRemoveBackup()

        self.fnLoadMolgroups()  # Load Molgroups into self.diMolgroups
        normarea = 100  # default, if no normarea is provided

        File = open('areatab.dat', 'w')  # save all molgroupdata in table for loading
        # into Igor

        File.write('z ')  # write header line
        for element in self.diMolgroups:
            File.write(self.diMolgroups[element]['headerdata']['ID'] + ' ')
        File.write("summol water waterperc")
        File.write('\n')

        element = list(self.diMolgroups.keys())[0]
        datalength = len(self.diMolgroups[element]['zaxis'])

        for i in range(datalength):
            File.write(str(self.diMolgroups[element]['zaxis'][i]) + ' ')
            summe = 0
            normarea = 0
            for element in self.diMolgroups:
                if element != 'normarea':
                    summe = summe + self.diMolgroups[element]['areaaxis'][i]
                else:
                    normarea = self.diMolgroups[element]['areaaxis'][i]
                File.write(str(self.diMolgroups[element]['areaaxis'][i]) + ' ')
            File.write(str(summe) + ' ' + str(normarea - summe) + ' ' + str((normarea - summe) / normarea))
            File.write('\n')

        File.close()

        plt.figure(1, figsize=(14, 10))

        plt.subplot(221)

        areasum = []  # add up area
        nSLsum = []  # calculate nSL
        zax = []
        for element in self.diMolgroups:
            area = []
            nSL = []
            zax = self.diMolgroups[element]['zaxis']
            stepsize = float(zax[1] - zax[0])
            length = len(zax)
            for i in range(length):
                area.append(self.diMolgroups[element]['areaaxis'][i])
                nSL.append(self.diMolgroups[element]['nslaxis'][i] * 1e4)
            plt.subplot(221)
            plt.plot(zax, area, label=element)
            plt.subplot(222)
            plt.plot(zax, nSL, label=element)

            if element != 'normarea':  # do not sum up normarea indicator
                if not areasum:
                    for i in range(length):
                        areasum.append(0.)
                        nSLsum.append(0.)
                for i in range(length):
                    areasum[i] = areasum[i] + area[i]
                    nSLsum[i] = nSLsum[i] + nSL[i]

            if element == 'normarea':  # Get normarea for nSLD calculations
                normarea = self.diMolgroups[element]['areaaxis'][0]

        plt.subplot(221)
        plt.plot(zax, areasum, label='Sum')
        plt.subplot(222)
        plt.plot(zax, nSLsum, label='Sum')

        nSLDSum = []  # calculate nSLD sum
        nSLDtVolFracSum = []  # calculate nSLD times vol frac sum
        for i in range(len(areasum)):
            nSLDSum.append(0)
            nSLDtVolFracSum.append(0)

        for element in self.diMolgroups:
            nSLDtVolFrac = []  # calculate nSLD times volfrac, at the moment times volfrac is not used
            if element != 'normarea':  # Get normarea for nSLD calculations
                for i in range(len(self.diMolgroups[element]['areaaxis'])):
                    area = self.diMolgroups[element]['areaaxis'][i]
                    if area:
                        nSLD = self.diMolgroups[element]['nslaxis'][i] / (stepsize * area) * 1e6
                    else:
                        nSLD = 0
                    nSLDtVolFrac.append(nSLD)  # *(area/normarea))
                    nSLDtVolFracSum[i] = nSLDtVolFracSum[i] + nSLDtVolFrac[i]
                plt.subplot(223)
                plt.plot(zax, nSLDtVolFrac, label=element)
        plt.subplot(223)
        # plt.plot(zax,nSLDtVolFracSum,label='nSLDtVolFracSum')

        for j in range(4):  # loop over contrast mixtures
            if j == 0:
                fBulknSLD = 6.34
            if j == 1:
                fBulknSLD = 4.00
            if j == 2:
                fBulknSLD = 0.00
            if j == 3:
                fBulknSLD = -0.56

            for i in range(length):  # calculate nSLD for several cases
                nSLDSum[i] = nSLsum[i] * 1E2 / (stepsize * normarea) + fBulknSLD * (1 - (areasum[i] / normarea))
            plt.subplot(224)
            plt.plot(zax, nSLDSum, label='nSLDsum CM' + str(fBulknSLD))

        plt.subplot(221)
        plt.ylabel('area / Ang+2')
        plt.xlabel('z / Ang')
        plt.legend(loc=0, prop=font)

        plt.subplot(222)
        plt.ylabel('nSL / 1e-4 Ang')
        plt.xlabel('z / Ang')
        plt.legend(loc=0, prop=font)

        plt.subplot(223)
        plt.ylabel('nSLD')  # * volfrac / 1e-6 Ang-2')
        plt.xlabel('z / Ang')
        plt.legend(loc=0, prop=font)

        plt.subplot(224)
        plt.ylabel('nSLD / 1e-6 Ang-2')
        plt.xlabel('z / Ang')
        plt.legend(loc=0, prop=font)

        plt.suptitle('%s \n\n Area and nSLD Profile' % (plotname))
        plt.savefig('%s_mol.png' % (plotname), format='png')
        plt.show()
        # plt.close()

    # -------------------------------------------------------------------------------

    def fnPlotFit(self, plotname):

        from matplotlib.font_manager import fontManager, FontProperties

        font = FontProperties(size='x-small')

        plt.figure(1, figsize=(14, 10))

        iCounter = 0
        while 1:
            sfilename = 'fit' + str(iCounter) + '.dat'
            if path.isfile(sfilename):
                file = open(sfilename, 'r')
                data = file.readlines()
                file.close()
                data = data[1:]

                k = 0
                l = 0
                qlist = []
                dqlist = []
                Rlist = []
                dRlist = []
                fitlist = []
                fitRFlist = []
                RFlist = []
                dRFlist = []
                reslist = []
                resplus = []
                resminus = []
                for line in data:
                    splitline = line.split()
                    qlist.append(float(splitline[0]))
                    dqlist.append(float(splitline[1]))
                    Rlist.append(float(splitline[2]))
                    dRlist.append(float(splitline[3]))
                    fitlist.append(float(splitline[4]))
                    RFlist.append(float(splitline[2]) * pow(float(splitline[0]), 4))
                    dRFlist.append(float(splitline[3]) * pow(float(splitline[0]), 4))
                    fitRFlist.append(float(splitline[4]) * pow(float(splitline[0]), 4))
                    reslist.append((float(splitline[2]) - float(splitline[4])) * pow(float(splitline[3]), -1))
                    resplus.append(1)
                    resminus.append(-1)

                plt.subplot(221)
                plt.errorbar(qlist, Rlist, yerr=dRlist, xerr=dqlist, fmt='.')
                plt.semilogy(qlist, fitlist, label='fit' + str(iCounter))
                plt.xlim(xmin=-0.01)

                plt.subplot(222)
                plt.errorbar(qlist, RFlist, yerr=dRFlist, xerr=dqlist, fmt='.')
                plt.semilogy(qlist, fitRFlist, label='fit' + str(iCounter))
                plt.xlim(xmin=-0.01)

                plt.subplot(223)
                plt.plot(qlist, reslist, label='fit' + str(iCounter))
                plt.plot(qlist, resplus, 'r')
                plt.plot(qlist, resminus, 'r')

                iCounter = iCounter + 1

            else:
                break

        iCounter = 0
        while 1:
            sfilename = 'profile' + str(iCounter) + '.dat'
            if path.isfile(sfilename):
                file = open(sfilename, 'r')
                data = file.readlines()
                file.close()
                data = data[1:]

                k = 0;
                l = 0
                zlist = [];
                rholist = []
                for line in data:
                    splitline = line.split()
                    zlist.append(float(splitline[0]))
                    rholist.append(float(splitline[1]) * 1e6)

                plt.subplot(224)
                plt.plot(zlist, rholist, label='profile' + str(iCounter))

                iCounter = iCounter + 1

            else:
                break

        plt.subplot(221)
        plt.ylabel('Reflectivity / R')
        plt.xlabel('q / Ang-1')
        plt.legend(loc=0, prop=font)

        plt.subplot(222)
        plt.ylabel('R*Q^4')
        plt.xlabel('q / Ang-1')
        plt.legend(loc=0, prop=font)

        plt.subplot(223)
        plt.ylabel('Normalized Residuals/ (R -F(q))/dR')
        plt.xlabel('q / Ang-1')
        plt.legend(loc=0, prop=font)

        plt.subplot(224)
        plt.ylabel('nSLD / 1e-6 Ang-2')
        plt.xlabel('z / Ang')
        plt.legend(loc=0, prop=font)

        plt.suptitle('%s \n\n Reflectivity data and fit, residuals and profile' % (plotname))
        plt.savefig('%s_fit.png' % (plotname), format='png')
        plt.show()
        # plt.close()

    # -------------------------------------------------------------------------------

    def fnPrintPar(self):  # prints parameters and their errors
        # from the covariance matrix on the screen

        litest = list(self.diParameters.keys())
        litest = sorted(litest, key=lambda keyitem: self.diParameters[keyitem]['number'])
        for parameter in litest:
            fRange = (self.diParameters[parameter]['upperlimit']
                      - self.diParameters[parameter]['lowerlimit'])
            fLowLim = self.diParameters[parameter]['lowerlimit']
            fValue = self.diParameters[parameter]['value']
            sRangeIndicator = ''
            for i in range(10):  # creates the par range ascii overview
                if ((fValue >= float(i) / 10 * fRange + fLowLim) and
                        (fValue < float(i + 1) / 10 * fRange + fLowLim)):
                    sRangeIndicator += '|'
                else:
                    if (fValue == float(i + 1) / 10 * fRange + fLowLim) and (i == 9):
                        sRangeIndicator += '|'
                    else:
                        sRangeIndicator += '.'
            print('%2i %25s  %s %15g +/- %g in [%g,%g]' % (self.diParameters[parameter]['number'],
                                                           parameter, sRangeIndicator, fValue,
                                                           self.diParameters[parameter]['error'],
                                                           fLowLim, self.diParameters[parameter]['upperlimit']))
        print('Chi squared: %g' % self.chisq)

    def fnPullMolgroup(self, liMolgroupNames, sparse=0):
        """
        Calls Function that recreates statistical data and extracts only area and nSL profiles for
        submolecular groups whose names are given in liMolgroupNames. Those groups
        are added for each iteration and a file pulledmolgroups.dat is created.
        A statistical analysis area profile containing the median, sigma, and
        2 sigma intervals are put out in pulledmolgroupsstat.dat.
        Saves results to file
        """
        diarea, dinsl, dinsld = self.fnPullMolgroupLoader(liMolgroupNames, sparse)
        diStat = self.fnPullMolgroupWorker(diarea, dinsl, dinsld)

        self.Interactor.fnSaveSingleColumns(self.mcmcpath+'/pulledmolgroups_area.dat', diarea)
        self.Interactor.fnSaveSingleColumns(self.mcmcpath+'/pulledmolgroups_nsl.dat', dinsl)
        self.Interactor.fnSaveSingleColumns(self.mcmcpath+'/pulledmolgroups_nsld.dat', dinsld)
        self.Interactor.fnSaveSingleColumns(self.mcmcpath+'/pulledmolgroupsstat.dat', diStat)

    def fnPullMolgroupLoader(self, liMolgroupNames, sparse=0):
        """
        Function recreates statistical data and extracts only area and nSL profiles for
        submolecular groups whose names are given in liMolgroupNames. Those groups
        are added for each iteration and a file pulledmolgroups.dat is created.
        A statistical analysis area profile containing the median, sigma, and
        2 sigma intervals are put out in pulledmolgroupsstat.dat.
        """

        try:
            self.diStatResults = self.fnLoadObject(self.mcmcpath+'/StatDataPython.dat')
            print('Loaded statistical data from StatDataPython.dat')

        except IOError:
            print('Failure to load StatDataPython.dat.')
            print('Recreate statistical data from sErr.dat.')
            self.fnRecreateStatistical(sparse=sparse)

        diarea = {
            'zaxis': self.diStatResults['Molgroups'][0][list(self.diStatResults['Molgroups'][0].keys())[0]]['zaxis']}
        dinsl = {
            'zaxis': self.diStatResults['Molgroups'][0][list(self.diStatResults['Molgroups'][0].keys())[0]]['zaxis']}
        dinsld = {
            'zaxis': self.diStatResults['Molgroups'][0][list(self.diStatResults['Molgroups'][0].keys())[0]]['zaxis']}
        for i, iteration in enumerate(self.diStatResults['Molgroups']):

            # add together all the molgroups that have to be analyzed
            sumareaprofile = []
            sumnslprofile = []
            for molgroup in liMolgroupNames:
                if molgroup in iteration:
                    sumareaprofile = [ii + jj for ii, jj in zip(sumareaprofile, iteration[molgroup]['areaaxis'])] if \
                        sumareaprofile else iteration[molgroup]['areaaxis']
                    sumnslprofile = [ii + jj for ii, jj in zip(sumnslprofile, iteration[molgroup]['nslaxis'])] if \
                        sumnslprofile else iteration[molgroup]['nslaxis']
                else:
                    if i == 0:
                        print('Molecular group %s does not exist.' % molgroup)
                    if not sumareaprofile:
                        sumareaprofile = [0 for _ in range(len(diarea['zaxis']))]
                        sumnslprofile = [0 for _ in range(len(diarea['zaxis']))]

                diarea['iter%i' % i] = sumareaprofile
                dinsl['iter%i' % i] = sumnslprofile
                stepsize = diarea['zaxis'][1] - diarea['zaxis'][0]
                dinsld['iter%i' % i] = []
                for j in range(len(diarea['iter%i' % i])):
                    if diarea['iter%i' % i][j] != 0:
                        dinsld['iter%i' % i].append(dinsl['iter%i' % i][j] / diarea['iter%i' % i][j] / stepsize)
                    else:
                        dinsld['iter%i' % i].append(0.0)

        return diarea, dinsl, dinsld

    def fnPullMolgroupWorker(self, diarea, dinsl, dinsld):
        """
        Function recreates statistical data and extracts only area and nSL profiles for
        submolecular groups whose names are given in liMolgroupNames. Those groups
        are added for each iteration and a file pulledmolgroups.dat is created.
        A statistical analysis area profile containing the median, sigma, and
        2 sigma intervals are put out in pulledmolgroupsstat.dat.
        """

        # do a statistics over every z-position
        diStat = dict(zaxis=[], m2sigma_area=[], msigma_area=[], median_area=[], psigma_area=[], p2sigma_area=[],
                      m2sigma_nsl=[], msigma_nsl=[], median_nsl=[], psigma_nsl=[], p2sigma_nsl=[],
                      m2sigma_nsld=[], msigma_nsld=[], median_nsld=[], psigma_nsld=[], p2sigma_nsld=[])
        for i in range(len(diarea[list(diarea.keys())[0]])):
            liOnePosition = [iteration[i] for key, iteration in diarea.items() if key != 'zaxis']
            stat = self.fnCalcConfidenceLimits(liOnePosition, method=1)
            diStat['zaxis'].append(str(diarea['zaxis'][i]))
            diStat['m2sigma_area'].append(stat[0])
            diStat['msigma_area'].append(stat[1])
            diStat['median_area'].append(stat[2])
            diStat['psigma_area'].append(stat[3])
            diStat['p2sigma_area'].append(stat[4])
            liOnePosition = [iteration[i] for key, iteration in dinsl.items() if key != 'zaxis']
            stat = self.fnCalcConfidenceLimits(liOnePosition, method=1)
            diStat['m2sigma_nsl'].append(stat[0])
            diStat['msigma_nsl'].append(stat[1])
            diStat['median_nsl'].append(stat[2])
            diStat['psigma_nsl'].append(stat[3])
            diStat['p2sigma_nsl'].append(stat[4])
            liOnePosition = [iteration[i] for key, iteration in dinsld.items() if key != 'zaxis']
            stat = self.fnCalcConfidenceLimits(liOnePosition, method=1)
            diStat['m2sigma_nsld'].append(stat[0])
            diStat['msigma_nsld'].append(stat[1])
            diStat['median_nsld'].append(stat[2])
            diStat['psigma_nsld'].append(stat[3])
            diStat['p2sigma_nsld'].append(stat[4])

        return diStat

    def fnRecreateStatistical(self, bRecreateMolgroups=True, sparse=0):

        # Recreates profile and fit data associated with stat file
        from sys import stdout

        self.fnLoadParameters()
        self.fnLoadStatData(sparse)

        problem = self.Interactor.fnRestoreFitProblem()

        j = 0
        self.diStatResults['nSLDProfiles'] = []  # delete list of all nSLD profiles
        self.diStatResults['Molgroups'] = []  # delete list of all molecular groups

        for iteration in range(self.diStatResults['NumberOfStatValues']):  # cycle through all individual stat. results
            try:
                # appends a new list for profiles for the current MC iteration
                self.diStatResults['nSLDProfiles'].append([])
                liParameters = list(self.diParameters.keys())
                # sort by number of appereance in setup.c
                liParameters = sorted(liParameters, key=lambda keyitem: self.diParameters[keyitem]['number'])
                bConsistency = True
                for element in liParameters:
                    if element not in list(self.diStatResults['Parameters'].keys()):
                        bConsistency = False
                if bConsistency:  # check for consistency
                    print('Processing parameter set %i.\n' % (j))
                    p = []
                    for parameter in liParameters:
                        val = self.diStatResults['Parameters'][parameter]['Values'][iteration]
                        # Rescaling is currently disabled
                        # if 'rho_' in parameter or 'background' in parameter:
                        #    val *= 1E6
                        p.append(val)
                    problem.setp(p)
                    # distinguish between FitProblem and MultiFitProblem
                    if 'models' in dir(problem):
                        for M in problem.models:
                            z, rho, irho = self.Interactor.fnRestoreSmoothProfile(M)
                            self.diStatResults['nSLDProfiles'][-1].append((z, rho, irho))
                            print(M.chisq())
                    else:
                        z, rho, irho = self.Interactor.fnRestoreSmoothProfile(problem)
                        self.diStatResults['nSLDProfiles'][-1].append((z, rho, irho))
                        print(problem.chisq())
                    stdout.flush()

                    if bRecreateMolgroups:
                        self.diMolgroups = self.Interactor.fnRestoreMolgroups(problem)
                        # append molgroup information to self.diStatResults
                        self.diStatResults['Molgroups'].append(self.diMolgroups)
                else:
                    raise RuntimeError('Statistical error data and setup file do not match')
            finally:
                j += 1

        self.fnSaveObject(self.diStatResults, self.mcmcpath + '/StatDataPython.dat')  # save stat data to disk

    def fnReplaceParameterLimitsInSetup(self, sname, flowerlimit, fupperlimit):  # scans setup.c file for parameter with
        file = open(self.setupfilename, 'r+')  # name sname and replaces the lower and
        data = file.readlines()  # upper fit limits by the given values
        file.close()
        smatch = compile(r'(pars_add\(pars.*?\"' + sname +
                         '.+?,.+?,).+?(,).+?(\))',
                         IGNORECASE | VERBOSE)
        newdata = []
        for line in data:
            newdata.append(smatch.sub(r'\1 ' + str(flowerlimit) + r'\2 '
                                      + str(fupperlimit) + r'\3', line))

        file = open(self.setupfilename, 'w')
        file.writelines(newdata)
        file.close()

    def fnRemoveBackup(self):  # deletes the backup directory
        self.fnRestoreBackup()
        call(['rm', '-rf', 'rsbackup'])

    def fnRestoreBackup(self):  # copies all files from the backup directory
        if path.isfile('rsbackup/pop.dat'):
            call(['cp', 'rsbackup/pop.dat', '.'])  # back to the working directory
        if path.isfile('rsbackup/par.dat'):
            call(['cp', 'rsbackup/par.dat', '.'])
        if path.isfile('rsbackup/covar.dat'):
            call(['cp', 'rsbackup/covar.dat', '.'])
        if path.isfile('rsbackup/fit*.dat'):
            call(['cp', 'rsbackup/fit*.dat', '.'])
        if path.isfile('rsbackup/fit'):
            call(['cp', 'rsbackup/fit', '.'])
        if path.isfile('rsbackup/fit*.dat'):
            call(['cp', 'rsbackup/model*.dat', '.'])
        if path.isfile('rsbackup/pop_bak.dat'):
            call(['cp', 'rsbackup/pop_bak.dat', '.'])
        if path.isfile('rsbackup/profile*.dat'):
            call(['cp', 'rsbackup/profile*.dat', '.'])
        if path.isfile('rsbackup/' + self.setupfilename):
            call(['cp', 'rsbackup/' + self.setupfilename, '.'])
        if path.isfile('rsbackup/setup.o'):
            call(['cp', 'rsbackup/setup.o', '.'])
        # -------------------------------------------------------------------------------

    def fnRestoreFileList(self, filelist):  # not used
        file = open(self.setupfilename, "r")
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
        file = open(self.setupfilename, "w")
        file.writelines(newdata)
        file.close()

    # -------------------------------------------------------------------------------

    def fnSaveObject(self, object, sFileName):

        import pickle

        File = open(sFileName, "wb")
        pickle.dump(object, File)
        File.close()

    # -------------------------------------------------------------------------------
    def fnSimulateReflectivity(self, qmin=0.008, qmax=0.325, s1min=0.108, s1max=4.397, s2min=0.108, s2max=4.397,
                               tmin=18, tmax=208, nmin=11809, rhomin=-0.56e-6, rhomax=6.34e-6, cbmatmin=1.1e-5,
                               cbmatmax=1.25e-6, mode='water', pre=1, qrange=0):
        # simulates reflectivity based on a parameter file called simpar.dat
        # requires a compiled and ready to go fit whose fit parameters are modified and fixed

        def convert_fit_to_sim_data(i):
            simdata = pandas.read_csv('fit' + str(i) + '.dat', sep=' ', header=None,
                                      names=['Q', 'dQ', 'R', 'dR', 'fit'],
                                      skip_blank_lines=True, comment='#')
            del simdata['dQ']
            del simdata['R']
            simdata = simdata.rename(columns={'fit': 'R'})
            simdata = simdata[['Q', 'R', 'dR']]
            return simdata

        def simulate_data():
            self.fnMake()  # compile changed setup.c
            call(["rm", "-f", "mol.dat"])
            call(["./fit", "-g"])  # write out profile.dat and fit.dat
            call(["sync"])  # synchronize file system
            sleep(1)  # wait for system to clean up

        def backup_simdat(i):
            if not path.isfile('simbackup' + str(i) + '.dat'):
                pr = Popen(["cp", 'sim' + str(i) + '.dat', 'simbackup' + str(i) + '.dat'])
                pr.wait()
            else:
                pr = Popen(["cp", 'simbackup' + str(i) + '.dat', 'sim' + str(i) + '.dat'])
                pr.wait()

        self.fnLoadParameters(self.sPath)  # Load Parameters and modify setup.cc
        self.fnBackup()  # Backup setup.c, and other files
        try:
            liParameters = list(self.diParameters.keys())  # get list of parameters from setup.c/par.dat
            liParameters = sorted(liParameters, key=lambda keyitem: self.diParameters[keyitem][
                'number'])  # sort by number of appereance in setup.c
            simpar = pandas.read_csv('simpar.dat', sep=' ', header=None, names=['par', 'value'], skip_blank_lines=True,
                                     comment='#')

            liAddition = []
            for parameter in liParameters:  # cycle through all parameters
                parvalue = simpar[simpar.par == parameter].iloc[0][1]
                strchange = ''
                parvaluefinal = parvalue
                if ('rho' in parameter) and fabs(parvalue) > 1E-4:
                    parvaluefinal = parvalue * 1E-6
                    strchange = ' => changed to ' + str(parvaluefinal)
                elif (
                        'background' in parameter) and parvalue > 1E-4:  # The issue with the background is two different uses
                    parvaluefinal = parvalue * 1E-6  # some setup.cc use negative numbers as parameter values and then
                    strchange = ' => changed to ' + str(
                        parvaluefinal)  # compute background=10^value, others use positive values and then
                    # compute background=value *1E-6
                    # this should catch it all
                print(str(parameter) + ' ' + str(parvalue) + strchange)

                liAddition.append(('%s = %s;\n' %  # change setup.c to quasi fix all parameters
                                   (self.diParameters[parameter]['variable'], parvaluefinal)))
            self.fnWriteConstraint2SetupC(liAddition)  # write out

            # in case q-range is changing, make sure to back up original sim dat or reload this to preserve any
            # changing q-steps in the original
            i = 0
            while path.isfile('fit' + str(i) + '.dat') and qrange != 0:
                backup_simdat(i)
                i += 1
            # create first set of simulated data
            simulate_data()

            # extend simulated data to q-range as defined
            i = 0
            while path.isfile('fit' + str(i) + '.dat') and qrange != 0:
                simdata = convert_fit_to_sim_data(i)
                # first cut data at qrange
                simdata = simdata[(simdata.Q <= qrange)]
                # now add data points in case the q-range is too short
                while simdata['Q'].iloc[-1] < qrange:
                    newframe = pandas.DataFrame([[2 * simdata['Q'].iloc[-1] - simdata['Q'].iloc[-2], 1, 1]],
                                                columns=['Q', 'R', 'dR'])
                    simdata = simdata.append(newframe)
                simdata.to_csv('sim' + str(i) + '.dat', sep=' ', index=None)
                i += 1

            # redo simulation, now with correct length of data
            simulate_data()

            i = 0
            last_defined_rho_solv = 0
            while path.isfile('fit' + str(i) + '.dat'):
                simdata = convert_fit_to_sim_data(i)

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

                simdata['dR'] = 0.0

                for index in simdata.index:
                    ns = I * simdata.iloc[index, 1] * c1 * c2 * (simdata.iloc[index, 0]) ** 2 * (
                                c3 + c4 * (simdata.iloc[index, 0]) ** 2)
                    nb = I * cbmat * c1 * c2 * (simdata.iloc[index, 0]) ** 2 * (c3 + c4 * (simdata.iloc[index, 0]) ** 2)
                    dRoR = sqrt(ns + 2 * nb) / (ns)
                    dR = simdata.iloc[index, 1] * dRoR  # see manuscript for details on calculation
                    simdata.iat[index, 2] = dR
                    simdata.iat[index, 1] = simdata.iloc[index, 1] + normalvariate(0,
                                                                                   1) * dR  # modify reflectivity within error bars
                    # if ns<0:
                    #    print index,ns,I,simdata.iloc[index,1],c1,c2,simdata.iloc[index,0],c3,c4,nb,dRoR,dR,simdata.iloc[index,2]

                simdata.to_csv('sim' + str(i) + '.dat', sep=' ', index=None)
                i += 1
        finally:
            self.fnRemoveBackup()

    # -------------------------------------------------------------------------------

    def fnStatTable(self, sTableName, fConfidence):

        def fnTexFormatf(fLow, fMed, fHigh):

            def fnDetPrec(fF):
                if fabs(fF) < 1e-4:  # applies to nSLDs
                    fF *= 1E6
                if fF > 0:  # determine precision
                    fPrec = ceil(log10(fF) * (-1))
                else:
                    fPrec = 0
                if fPrec > 0:  # takes care of numbers like fF=0.0095
                    iPrec = int(fPrec)
                    if round(fF, iPrec) == round(fF, iPrec - 1):  # which should be rounded to 0.01 and
                        fPrec -= 1  # not to 0.010
                return fF, fPrec

            fLowDiff, fLowPrec = fnDetPrec(fMed - fLow)
            fHighDiff, fHighPrec = fnDetPrec(fHigh - fMed)
            fMed, fMedPrec = fnDetPrec(fMed)
            fPrec = (max(fLowPrec, fHighPrec)) + 1.0
            iPrec = int(fPrec)

            fLowDiff = round(fLowDiff + 0.5 * pow(10, (-1) * (fPrec + 1.0)), iPrec)  # conservative rounding
            fHighDiff = round(fHighDiff + 0.5 * pow(10, (-1) * (fPrec + 1.0)), iPrec)
            fMed = round(fMed, iPrec)
            return fLowDiff, fMed, fHighDiff, iPrec

        self.fnAnalyzeStatFile(fConfidence)  # analyze stat data and
        # populate self.diStatresults

        file = open(sTableName, 'r')  # load in template
        template = file.readlines()
        file.close()

        table = []  # table to be created

        for line in template:
            splitline = line.split()
            for i, phrasetex in enumerate(splitline):  # look at each string in template
                phrase = phrasetex.replace('\\', '')  # remove all '\'
                if phrase in self.diStatResults['Parameters']:  # if it resembles a paramter name -> replace
                    fMedian = self.diStatResults['Parameters'][phrase]['Median']
                    fLowPerc = self.diStatResults['Parameters'][phrase]['LowPerc']
                    fHighPerc = self.diStatResults['Parameters'][phrase]['HighPerc']
                    fLowPercDiff, fMedianDiff, fHighPercDiff, iPrec = fnTexFormatf(fLowPerc, fMedian, fHighPerc)
                    sPrec = '%(#).' + str(iPrec) + 'f'
                    temp = '$'  # Latex format
                    temp += (sPrec % {'#': fMedianDiff})
                    temp += '_{'
                    temp = temp + '-' + (sPrec % {'#': fLowPercDiff})
                    temp += '}^{'
                    temp = temp + '+' + (sPrec % {'#': fHighPercDiff})
                    temp += '}$'
                    splitline[i] = temp
            table.append(' '.join(splitline) + '\n')

        file = open("StatTable.tex", 'w')  # save table to file
        file.writelines(table)
        file.close()

    # -------------------------------------------------------------------------------
    def fnTruncateParDat(self):
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

            # -------------------------------------------------------------------------------

    def fnWriteConstraint2SetupC(self, liExpression):
        # Writes a list of expressions at the beginning of
        # the constraint section in setup.c
        File = open(self.setupfilename, "r")  # open setup.c
        data = File.readlines()
        File.close()
        newdata = []

        for i in range(len(data)):  # need integer index for slicing
            if 'void constr_models(' in data[i]:  # searching for function declaration
                newdata = data[:i + 2] + liExpression + data[i + 2:]  # inserting Expression two lines later,
                # after the initial '{' of the function

                break  # break loop, only one insertion
        File = open(self.setupfilename, "w")  # write changed setup.c
        File.writelines(newdata)
        File.close()

    # -------------------------------------------------------------------------------

    def fnWriteOutGareflModel(self):
        pr = Popen(["./fit", "-g"])  # write out profile.dat and fit.dat
        pr.wait()


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

class CSAXS:
    def __init__(self):
        self.datainteractor = CDataInteractor()

    def fnMonteCarloPofR(self, cfilename='saxs.dat', iiterations=100):
        for i in range(iiterations):
            self.datainteractor.fnMCModifyFile([cfilename])
            # to be continued ...


# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------


def Auto(convergence=0.001):  # automatic fit
    ReflPar.fnMake()
    call(['rm', '-f', 'par.dat'], stdout=open(devnull, "w"))
    call(['rm', '-f', 'pop.dat'], stdout=open(devnull, "w"))
    call(['rm', '-f', 'covar.dat'], stdout=open(devnull, "w"))
    call(['nice', './fit', '-eS', '-n', '21'], stdout=open(devnull, "w"))
    call(['cp', 'pop_bak.dat', 'pop.dat'])  # copy newest population into pop.dat
    print('Genetic run, approximate roughness')
    ReflPar.fnLoadAndPrintPar()
    fOldChiSq = ReflPar.fnGetChiSq()

    while 1:  # genetic runs until chisq<20 or not improving
        call(['nice', './fit', '-peS', '-n', '51'], stdout=open(devnull, "w"))
        call(['cp', 'pop_bak.dat', 'pop.dat'])
        print('Genetic run, approximate roughness')
        ReflPar.fnLoadAndPrintPar()
        if (ReflPar.fnGetChiSq() < 10 or (fOldChiSq - ReflPar.fnGetChiSq()) < convergence):
            call(['cp', 'pop_bak.dat', 'pop.dat'])
            break
        fOldChiSq = ReflPar.fnGetChiSq()
    AutoFinish(convergence)  # AutoFinish takes over
    ReflPar.fnTruncateParDat()


# -------------------------------------------------------------------------------

def AutoFinish(convergence=0.001, sPath='./'):
    if sPath != './':
        ReflPar.fnImportMCMCBestFit(sPath)
        # automatic fit starting with Amoeba, approximate roughness
    ReflPar.fnLoadParameters()
    fOldChiSq = ReflPar.fnGetChiSq()

    while 1:  # Amoeba, approximate roughness
        call(['nice', './fit', '-peaS'], stdout=open(devnull, "w"))  # until chisq<10 or no improvement
        print('Amoeba, approximate roughness')
        ReflPar.fnLoadAndPrintPar()
        call(['cp', 'pop_bak.dat', 'pop.dat'])  # copy newest population into pop.dat
        if ReflPar.chisq < 5 or (fOldChiSq - ReflPar.fnGetChiSq() < convergence):
            break
        fOldChiSq = ReflPar.fnGetChiSq()
    AutoFinish2(convergence)  # AutoFinish2 takes over


# -------------------------------------------------------------------------------


def AutoFinish2(convergence=0.001, sPath='./'):  # automatic fit, only local minimum refinement

    if sPath != './':
        ReflPar.fnImportMCMCBestFit(sPath)

    ReflPar.fnLoadParameters()

    fOldChiSq = 1E24
    fBlockChiSq = 5E23
    fTempChiSq = 1E23
    iGeneticIsUseless = False
    iGeneticSkipCount = 0
    iFinishCounter = 0

    while 1:

        if iGeneticIsUseless:  # skip genetic algorithm after it proved to be
            iGeneticSkipCount += 1  # useless and give it a chance every 5 iterations
            print(' ')
            print('Genetic algorithm skipped')

        call(['cp', 'pop.dat', 'pop_rspy.dat'])
        while (not iGeneticIsUseless) and iGeneticSkipCount < 6:
            iGeneticIsUseless = False
            iGeneticSkipCount = 0
            call(['nice', './fit', '-pe', '-n', '21'], stdout=open(devnull, "w"))
            print('Genetic run, correct roughness')
            ReflPar.fnLoadAndPrintPar()
            fTempChiSq2 = ReflPar.fnGetChiSq()
            if fTempChiSq - fTempChiSq2 < convergence:
                if fTempChiSq2 < fTempChiSq:
                    fTempChiSq = fTempChiSq2
                    call(['cp', 'pop_bak.dat', 'pop.dat'])
                else:
                    pass
                    # iGeneticIsUseless=True
                break
            else:
                fTempChiSq = fTempChiSq2
                call(['cp', 'pop_bak.dat', 'pop.dat'])
        if fTempChiSq > fBlockChiSq:
            call(['cp', 'pop_rspy.dat', 'pop.dat'])
            fTempChiSq = fBlockChiSq

        call(['cp', 'pop.dat', 'pop_rspy.dat'])
        while 1:
            call(['nice', './fit', '-pea'], stdout=open(devnull, "w"))
            print('Amoeba, correct roughness')
            ReflPar.fnLoadAndPrintPar()
            fTempChiSq2 = ReflPar.fnGetChiSq()
            if fTempChiSq - fTempChiSq2 < convergence:
                if fTempChiSq2 < fTempChiSq:
                    fTempChiSq = fTempChiSq2
                    call(['cp', 'pop_bak.dat', 'pop.dat'])
                break
            else:
                call(['cp', 'pop_bak.dat', 'pop.dat'])
                fTempChiSq = fTempChiSq2
        if fTempChiSq > fBlockChiSq:
            call(['cp', 'pop_rspy.dat', 'pop.dat'])
            fTempChiSq = fBlockChiSq

        call(['cp', 'pop.dat', 'pop_rspy.dat'])
        while 1:
            call(['nice', './fit', '-pel'], stdout=open(devnull, "w"))
            print('Levenberg-Marquardt, correct roughness')
            ReflPar.fnLoadAndPrintPar()
            fTempChiSq2 = ReflPar.fnGetChiSq()
            if fTempChiSq - fTempChiSq2 < convergence:
                if fTempChiSq2 < fTempChiSq:
                    fTempChiSq = fTempChiSq2
                    call(['cp', 'pop_bak.dat', 'pop.dat'])
                break
            else:
                call(['cp', 'pop_bak.dat', 'pop.dat'])
                fTempChiSq = fTempChiSq2
        if fTempChiSq > fBlockChiSq:
            call(['cp', 'pop_rspy.dat', 'pop.dat'])
            fTempChiSq = fBlockChiSq

        print('old ChiSq: %g new ChiSq: %g' % (fOldChiSq, fTempChiSq))
        if (fOldChiSq - fTempChiSq) < convergence:
            iFinishCounter += 1
        else:
            iFinishCounter = 0

        if iFinishCounter == 2:
            break

        fOldChiSq = fTempChiSq


# -------------------------------------------------------------------------------
def AvgProfile():
    iCounter = 0
    while 1:
        sfilename = 'ContDimRho' + str(iCounter) + '.dat'
        if path.isfile(sfilename):
            file = open(sfilename, 'r')
            data = file.readlines()
            file.close()
            rho = []
            for i in range(len(data)):
                rho.append(float(data[i]))

            sfilename = 'ContDimZ' + str(iCounter) + '.dat'
            file = open(sfilename, 'r')
            data = file.readlines()
            file.close()
            z = []
            for i in range(len(data)):
                z.append(float(data[i]))

            ztemp = []  # correct for Igor-style array axis notation
            for i in range(len(z) - 1):
                ztemp.append((z[i] + z[i + 1]) / 2)
            z = ztemp

            sfilename = 'ContArray' + str(iCounter) + '.dat'
            file = open(sfilename, 'r')
            data = file.readlines()
            file.close()
            arr = []
            for i in range(len(data)):
                rhoarr = []
                tdata = (data[i]).split()
                for j in range(len(tdata)):
                    rhoarr.append(float(tdata[j]))
                arr.append(rhoarr)
            # print len(tdata), len(rho), len(z), len(data)

            median = []
            maxlikely = []
            lowperc = []
            highperc = []
            for i in range(len(z)):
                cumulative = numpy.cumsum(arr[i])
                mediancum = cumulative[-1] * 0.5
                lowperccum = cumulative[-1] * 0.17
                highperccum = cumulative[-1] * 0.83
                for j in range(len(arr[i]) - 1):
                    if cumulative[j] <= mediancum <= cumulative[j + 1]:
                        frac = (mediancum - cumulative[j]) / (cumulative[j + 1] - cumulative[j])
                        median.append(rho[j] * frac + rho[j - 1] * (1 - frac))
                    if cumulative[j] <= lowperccum <= cumulative[j + 1]:
                        frac = (lowperccum - cumulative[j]) / (cumulative[j + 1] - cumulative[j])
                        lowperc.append(rho[j] * frac + rho[j - 1] * (1 - frac))
                    if cumulative[j] <= highperccum <= cumulative[j + 1]:
                        frac = (highperccum - cumulative[j]) / (cumulative[j + 1] - cumulative[j])
                        highperc.append(rho[j] * frac + rho[j - 1] * (1 - frac))
                    if max(arr[i]) == arr[i][j]:
                        maxlikely.append(rho[j])

            sfilename = 'AvgProfile' + str(iCounter) + '.dat'
            file = open(sfilename, "w")
            file.write('z    maxlikely    median    lowperc     highperc \n')
            for i in range(len(z)):
                file.write(
                    str(z[i]) + ' ' + str(maxlikely[i]) + ' ' + str(median[i]) + ' ' + str(lowperc[i]) + ' ' + str(
                        highperc[i]) + '\n')
            file.close()
            iCounter += 1
        else:
            break


# -------------------------------------------------------------------------------

def fContour(dZGrid=0.5, dRhoGrid=1e-8, sStatMode='is'):  # Calculate contour plot data from SErr.dat
    ReflPar.fnContourData(sStatMode + 'Err.dat', dZGrid, dRhoGrid)


# -------------------------------------------------------------------------------
def fCalculateMolgroups(fConfidence):
    ReflPar.fnCalculateMolgroupProperty(fConfidence)


# -------------------------------------------------------------------------------
def fnDetermineFitSoftware():
    if path.isfile('run.py'):
        print('Refl1D setup identified.')
        return 'refl1d'
    else:
        print('Garefl setup identified.')
        return 'garefl'


# -------------------------------------------------------------------------------
def fnDetermineStatFile():
    if path.isfile('isErr.dat'):
        return 'isErr.dat'
    elif path.isfile('sErr.dat'):
        return 'sErr.dat'
    else:
        return ''


# -------------------------------------------------------------------------------

def DErr(convergence=0.001):  # function has not yet been thoroughly tested
    # Error analysis by stepwise parameter displacement,
    # fixing, and fitting within a defined range
    def fnWriteToFile(sParameterName, fTestValue, fChiSq):  # append new data set to file
        try:
            file = open('Error_' + sParameterName + '.dat', 'r')  # if file exists, open and read data
            data = file.readlines()
            file.close()
        except IOError:
            data = []  # otherwise start with an empty data file
        newdata = data[:]  # copy dat into new object
        newdata.append(str(fTestValue) + "  " + str(fChiSq) + '\n')  # append line with the data parameter to store
        file = open('Error_' + sParameterName + '.dat', 'w')  # create new file and write out the new data
        file.writelines(newdata)
        file.close()

    ReflPar.fnLoadAndPrintPar()  # print latest fit parameters
    tTaggedParameters = ReflPar.fnGetTaggedParameters()  # get a list of all tagged parameters and test ranges
    ReflPar.fnBackup()  # backup working directory
    try:
        for tParameter in tTaggedParameters:  # iterate through all parameters
            fValue = ReflPar.fnGetParameterValue(tParameter[0])  # get fitted par value, fit constraints,
            sParameterName = tParameter[0]
            fLowerLimit = float(tParameter[1])
            fUpperLimit = float(tParameter[2])
            fTestLLimit = float(tParameter[3])  # test ranges
            fTestULimit = float(tParameter[4])
            fTestStep = float(tParameter[5])
            for iTestPoint in range(int(fTestLLimit / fTestStep),  # iterate through test points, which are centered
                                    int(fTestULimit / fTestStep)):  # around the fit value
                fTestValue = fValue + float(iTestPoint) * fTestStep
                if fTestValue < fLowerLimit or fTestValue > fUpperLimit:
                    continue
                print(sParameterName, fTestValue)
                ReflPar.fnReplaceParameterLimitsInSetup(sParameterName,
                                                        # replace fit constraint in setup.c with a range, which
                                                        0.9999 * fTestValue, 1.0001 * fTestValue)  # is
                pr = Popen(["make"])
                pr.wait()
                AutoFinish(convergence)  # fit using the stored minimum
                pr = Popen(["cp", "pop_bak.dat", "pop.dat"])
                pr.wait()
                AutoFinish2(convergence)  # refit, because sometimes the fit gets trapped
                pr = Popen(["cp", "pop_bak.dat", "pop.dat"])
                pr.wait()
                ReflPar.fnLoadParameters()
                fnWriteToFile(sParameterName, fTestValue, ReflPar.fnGetChiSq())
            ReflPar.fnRestoreBackup()
    finally:
        ReflPar.fnRemoveBackup()  # in case of crash restore working dir


# -------------------------------------------------------------------------------

def fEnvelope(fZGrid=0.1, fSigma=0, sStatMode='is', bShortFlag=False,
              iContrast=-1):  # Calculate nSLD envelopes plot data from SErr.dat
    ReflPar.fnnSLDEnvelopes(fZGrid, fSigma, sStatMode + 'Err.dat', bShortFlag, iContrast)


# -------------------------------------------------------------------------------


def fMCMultiCore(iIterations=1000, fMCConvergence=0.01, iConvergence=0.01,
                 fConfidence=0.9546, sMode='is', iNodes=1, iIterationsPerCall='10',
                 sJobID=''):  # Monte Carlo for Multicore systems

    def fCleanUpDirectory(sDir):
        call(['rm', '-rf', sDir])  # remove rs directory

    def fCreateRsDir(sMode, lSubProcesses):
        while 1:  # create random name directory
            sRandomName = str(int(random() * 1000000000000))  # and check if the name already
            sDirName = '../rspy' + sRandomName  # exists as a Multicore DirName
            iNameExists = 0
            for element in lSubProcesses:
                if sDirName == element[1]:
                    iNameExists = 1
            if not iNameExists == 1:
                break

        call('mkdir ' + sDirName, shell=True)

        call('cp * ' + sDirName + '/', shell=True)  # copy whole main directory over
        call('rm -f ' + sDirName + '/' + sMode + 'Err.dat', shell=True)  # not including the stat data calculated
        # so far

        if path.isfile(sDirName + '/' + 'StatFinished'):
            call(["rm", "-f", sDirName + '/' + "StatFinished"])
        if path.isfile(sDirName + '/' + 'StatAborted'):
            call(["rm", "-f", sDirName + '/' + "StatAborted"])
        return sDirName

    def fKillAllProcesses():  # kills all running processes
        for item in lSubProcesses:
            try:
                kill(int(item[0]), 15)
            except:
                pass
            print('Delete directories ...')
            fCleanUpDirectory(item[1])

    lSubProcesses = []
    iDoneIterations = 0
    iChange = 0
    fTimeAverage = 0.
    iExitFlag = 0
    # ReflPar.fnCheckFit()

    try:

        while (iDoneIterations < iIterations) and (iExitFlag == 0):
            # start new subprocesses
            while (iDoneIterations < (iIterations - len(lSubProcesses) * iIterationsPerCall)) and (
                    len(lSubProcesses) < iNodes):
                sDirName = fCreateRsDir(sMode,
                                        lSubProcesses)  # sDirName is name of directory for Multicore architecture
                # and a Filepointer for the PBS architecture

                pid = str(Popen([sDirName + '/rs.py',
                                 '-' + sMode, str(iIterationsPerCall), str(iConvergence)],
                                cwd=sDirName, stdout=open(devnull, "w"), stderr=open(devnull, "w")).pid)
                lSubProcesses.append((pid, sDirName, time()))
                iChange = 1
                sleep(2)  # wait for ga_refl random number generator

            if iChange == 1:  # is there any change to report?
                iChange = 0
                fActualTime = time()
                fProjectedTime = 0
                fInProcessTime = 0
                for element in lSubProcesses:
                    fInProcessTime = fInProcessTime + (
                            fActualTime - element[2])  # add up time already spent in not finished proc.

                fProjectedTime = (
                                         iIterations - iDoneIterations) * fTimeAverage  # total cpu time for remaining iterations
                fProjectedTime -= fInProcessTime  # credit for total cpu time alrady done
                fProjectedTime /= len(lSubProcesses)  # get real time

                if fProjectedTime < 0:
                    fProjectedTime = 0
                lTD = gmtime(fProjectedTime)
                lTA = gmtime(fTimeAverage)
                print('')
                print(ctime(fActualTime))
                print('-%s Monte Carlo Error analysis using %i processes.' % (sMode, iNodes))
                print('Computing %i iterations per call' % iIterationsPerCall)
                print('%i of %i iterations done.' % (iDoneIterations, iIterations))
                print('')
                print('Process list:')
                for i in range(len(lSubProcesses)):
                    lTS = gmtime(fActualTime - lSubProcesses[i][2])
                    print('%i: PID %s in directory %s running for %id %ih %imin %is' % (i + 1, lSubProcesses[i][0],
                                                                                        lSubProcesses[i][1],
                                                                                        lTS[2] - 1,
                                                                                        lTS[3], lTS[4], lTS[5]))
                if fTimeAverage > 0:
                    print('')
                    print('Average time per iteration: %id %ih %imin %is' % (lTA[2] - 1, lTA[3], lTA[4], lTA[5]))
                    print('Projected finish in %id %ih %imin %is' % (lTD[2] - 1, lTD[3], lTD[4], lTD[5]))
                    print('on %s' % (ctime(fActualTime + fProjectedTime)))
                print('')
                print('')

            sleep(30)  # wait before checking for finished subprocesses

            while 1:
                iFinishedProcessFound = 0

                for i in range(len(lSubProcesses)):  # check for finished sub processes
                    sDirName = lSubProcesses[i][1]
                    pid = lSubProcesses[i][0]

                    if path.isfile(sDirName + '/' + 'StatFinished'):  # look up if process is finished
                        if path.isfile(sMode + 'Err.dat'):
                            file = open(sDirName + '/' + sMode + 'Err.dat', 'r')  # get statistical data from rs subdir
                            data = file.readlines()
                            file.close()
                            data = data[1:]  # delete headerline
                            iFailureCounter = 0
                            while 1:
                                try:
                                    file = open(sMode + 'Err.dat', 'a')  # append new data to file in root dir
                                    file.writelines(data)
                                    file.close()
                                    break
                                except:
                                    print('(i)sErr.dat is in use. Wait for 2s')
                                    iFailureCounter += 1
                                    sleep(2)
                                    if iFailureCounter == 15:
                                        print('Cannot append to (i)sErr.dat -> Abort.')
                                        break
                        else:
                            call(['cp', sDirName + '/' + sMode + 'Err.dat', '.'])

                        fCleanUpDirectory(sDirName)

                        iDoneIterations = iDoneIterations + iIterationsPerCall
                        iChange = 1
                        iFinishedProcessFound = 1
                        fDeltaTime = time() - lSubProcesses[i][2]
                        fTimeAverage = fTimeAverage * (float(iDoneIterations -
                                                             iIterationsPerCall)) / float(
                            iDoneIterations) + fDeltaTime / float(iDoneIterations)
                        del lSubProcesses[i]  # remove entry from list

                        try:
                            ReflPar.fnAnalyzeStatFile(
                                fConfidence)  # see if convergence criterium for whole MC had been reached
                            if ReflPar.diStatResults['Convergence'] <= fMCConvergence:
                                print('MC has converged ..')
                                iExitFlag = 1
                        except:
                            print('Analysis failed...')

                        break  # because we changed len(lSubProcesses)

                    if path.isfile(sDirName + '/' + 'StatAborted'):  # look up if process is finished
                        fCleanUpDirectory(sDirName)
                        if (time() - lSubProcesses[i][2]) < 180:
                            iExitFlag = 1
                            print('=========Multicore Error=========')
                            print('Process termination within 3 min.')
                            print('=================================')
                        del lSubProcesses[i]  # remove entry from list
                        iChange = 1
                        iFinishedProcessFound = 1
                        break  # because we changed len(lSubProcesses)

                if (iFinishedProcessFound == 0) or (len(lSubProcesses) == 0) or (iExitFlag == 1):
                    break

        if not sJobID == '':
            file = open(sJobID, 'w')  # indicate Multicore is finished if required
            file.write('Multicore finished \n')  # requirement comes from presence of a
            file.close()  # job id as given by a PBS caller


    finally:
        print('Exiting MC')
        print(iExitFlag, iDoneIterations, iIterations)
        fKillAllProcesses()


# -------------------------------------------------------------------------------
def fMCMC(iMaxIterations=1024000, liMolgroups=['protein'], fSparse=0):
    while True:

        if not path.isfile('run.py'):  # make sure there is a run.py
            file = open('run.py', 'w')
            file.write('from bumps.fitproblem import FitProblem\n')
            file.write('from refl1d import garefl\n')
            file.write('from refl1d.names import *\n')
            file.write('\n')
            file.write("problem = garefl.load('model.so')\n")
            file.write('\n')
            file.write('problem.penalty_limit = 50\n')
            file.close()

        iBurn = 4000  # check wether an existing MCMC exists
        bMCMCexists = False
        for i in range(1, 9):
            iBurn = iBurn * 2
            if path.isdir('MCMC_' + str(iBurn) + '_500'):
                print('Found ' + 'MCMC_' + str(iBurn) + '_500 \n')
                bMCMCexists = True
                break

                # if not bMCMCexists:                                                         #run best fit
                # Auto()

        lCommand = ['refl1d_cli.py', 'run.py', '--fit=dream', '--parallel=0']
        if bMCMCexists:
            if iBurn >= iMaxIterations:
                print('Maximum number of MCMC iterations reached\n')
                break  # end
            lCommand.append('--resume=MCMC_' + str(iBurn) + '_500')
            lCommand.append('--store=MCMC_' + str(iBurn * 2) + '_500')
            lCommand.append('--burn=' + str(iBurn))
        else:
            lCommand.append('--init=lhs')
            lCommand.append('--store=MCMC_8000_500')
            lCommand.append('--burn=8000')
            iBurn = 4000

        lCommand.append('--steps=500')

        call(lCommand)  # run MCMC

        rename('MCMC_' + str(iBurn * 2) + '_500', 'MCMC')  # create sErr.dat
        if path.isfile('isErr.dat'):
            remove('isErr.dat')
        if path.isfile('sErr.dat'):
            remove('sErr.dat')
        StatAnalysis(-1, 0.005)  # sErr.dat contains about 1000 iterations
        rename('MCMC', 'MCMC_' + str(iBurn * 2) + '_500')

        if liMolgroups:
            if path.isfile('StatDataPython.dat'):
                remove('StatDataPython.dat')
            ReflPar.fnPullMolgroup(liMolgroups, 0)

        if bMCMCexists:
            rmtree('MCMC_' + str(iBurn) + '_500')

    return


# -------------------------------------------------------------------------------

def fMCPBS(iIterations, iConvergence=0.01, fConfidence=0.9546, sMode='is', iNodes=1,
           iIterationsPerCall='10', sQueue='default',
           iNumberOfPBSJobsInQueue=1):  # Monte Carlo for PBS batch system submission
    # fConfidence not yet used

    def fCleanUpDirectory(sJobID, iNotOutputFiles=0):
        call(['rm', '-f', 'run' + sJobID + ".sh"])  # delete run.sh files
        call(['rm', '-f', sJobID])  # delete file pointer for sucessfull PBS run

        if iNotOutputFiles == 0:  # and output files
            sOutName = sJobID + '.out'
            sErrName = sJobID + '.err'
            iSleepCounter = 0
            while 1:
                if path.isfile(sOutName) and path.isfile(sErrName):  # wait for output files
                    call(['rm', '-f', sOutName])
                    call(['rm', '-f', sErrName])
                    break
                retcode = call('ls ' + path.expanduser('~/') + '*.OU',
                               shell=True, stdout=open(devnull, "w"))
                if retcode == 0:  # on some systems only
                    break  # standardfiles are written
                sleep(1)
                iSleepCounter = iSleepCounter + 1
                if iSleepCounter > 20:
                    print('Waited 20s for output files to be written ... giving up.')
                    break

    def fCreateMultiCoreJobID(lSubProcesses):
        while 1:  # create random name directory
            sRandomName = str(int(random() * 1000000000000))  # and check if the name already
            iNameExists = 0  # as a PBS Job ID
            for element in lSubProcesses:
                if sRandomName == element[1]:
                    iNameExists = 1
            if not iNameExists == 1:
                break
        return sRandomName  # randomname will be job id

    def fKillAllProcesses():  # kills all running processes
        for item in lSubProcesses:
            call(['qdel', item[0]])  # delete submitted job
            sleep(2)  # give system time
            print('Delete directories ...')
            fCleanUpDirectory(item[1])

    lSubProcesses = []
    iDoneIterations = 0
    iChange = 0
    fTimeAverage = 0.
    iExitFlag = 0

    try:

        while (iDoneIterations < iIterations) and (iExitFlag == 0):
            # start new subprocesses
            while (iDoneIterations < (iIterations - len(lSubProcesses) * iIterationsPerCall)) and (
                    len(lSubProcesses) < iNumberOfPBSJobsInQueue):
                sJobID = fCreateMultiCoreJobID(
                    lSubProcesses)  # sDirName is name of directory for Multicore architecture
                # and a Filepointer for the PBS architecture
                data = []  # create run batchfile
                # data.append('#PBS -l nodes=1:ppn=1\n')
                data.append('#PBS -l ncpus=1\n')
                data.append('#PBS -l cput=48:00:00\n')
                data.append('#PBS -l walltime=72:00:00\n')
                data.append('#PBS -e ' + getcwd() + '/' + sJobID + '.err\n')
                data.append('#PBS -o ' + getcwd() + '/' + sJobID + '.out\n')
                data.append('cd $PBS_O_WORKDIR\n')
                data.append('./rs.py -' + sMode + ' ' + str(iIterationsPerCall) + ' '
                            + '-c' + ' ' + str(iConvergence) + ' ' + '-m' + ' ' + str(iNodes) + ' '
                            + '-ipc' + ' ' + str(iIterationsPerCall / iNodes) + ' ' + '-id' + ' ' + sJobID + '\n')
                data.append('#end')
                runfile = 'run.sh'
                file = open(runfile, "w")
                file.writelines(data)
                file.close()
                pid = Popen(['qsub', '-q', sQueue, runfile], open(devnull, "w")).communicate()[
                    0]  # ged pid from queue systm
                pid = pid.split()[0]  # remove newline at end of output

                lSubProcesses.append((pid, sJobID, time()))
                iChange = 1
                sleep(2)  # wait for ga_refl random number generator

            if iChange == 1:  # is there any change to report?
                iChange = 0
                fActualTime = time()
                fProjectedTime = 0
                fInProcessTime = 0
                for element in lSubProcesses:
                    fInProcessTime = fInProcessTime + (
                            fActualTime - element[2])  # add up time already spent in not finished proc.

                fProjectedTime = (
                                         iIterations - iDoneIterations) * fTimeAverage  # total cpu time for remaining iterations
                fProjectedTime = fProjectedTime - fInProcessTime  # credit for total cpu time alrady done
                fProjectedTime = fProjectedTime / len(lSubProcesses)  # get real time

                if fProjectedTime < 0:
                    fProjectedTime = 0
                lTD = gmtime(fProjectedTime)
                lTA = gmtime(fTimeAverage)
                print('')
                print(ctime(fActualTime))
                print('-%s Monte Carlo Error analysis using %i PBS jobs in queue.' % (sMode, iNumberOfPBSJobsInQueue))
                print('Computing %i iterations per call' % iIterationsPerCall)
                print('%i of %i iterations done.' % (iDoneIterations, iIterations))
                print('')
                print('Process list:')
                for i in range(len(lSubProcesses)):
                    lTS = gmtime(fActualTime - lSubProcesses[i][2])
                    print('%i: PID %s rs.py ID %s running for %id %ih %imin %is' % (i + 1, lSubProcesses[i][0],
                                                                                    lSubProcesses[i][1], lTS[2] - 1,
                                                                                    lTS[3], lTS[4], lTS[5]))
                if fTimeAverage > 0:
                    print('')
                    print('Average time per iteration: %id %ih %imin %is' % (lTA[2] - 1, lTA[3], lTA[4], lTA[5]))
                    print('Projected finish in %id %ih %imin %is' % (lTD[2] - 1, lTD[3], lTD[4], lTD[5]))
                    print('on %s' % (ctime(fActualTime + fProjectedTime)))
                print('')
                print('')

            sleep(30)  # wait before checking for finished subprocesses

            while 1:

                iFinishedProcessFound = 0

                for i in range(len(lSubProcesses)):  # check for finished sub processes
                    sJobID = lSubProcesses[i][1]
                    pid = lSubProcesses[i][0]

                    if path.isfile(sJobID):  # look up if process is finished
                        fCleanUpDirectory(sJobID)

                        iDoneIterations = iDoneIterations + iIterationsPerCall
                        iChange = 1
                        iFinishedProcessFound = 1
                        fDeltaTime = time() - lSubProcesses[i][2]
                        fTimeAverage = fTimeAverage * (float(iDoneIterations -
                                                             iIterationsPerCall)) / float(
                            iDoneIterations) + fDeltaTime / float(iDoneIterations)
                        del lSubProcesses[i]  # remove entry from list
                        break  # because we changed len(lSubProcesses)

                        # check for incorrectly finished subprocesses
                    retcode = call(['qstat', pid],
                                   stdout=open(devnull, "w"))
                    if retcode != 0:  # check for orphaned processes
                        fCleanUpDirectory(sJobID, 1)
                        if (time() - lSubProcesses[i][2]) < 180:
                            iExitFlag = 1
                            print('============PBS Error============')
                            print('Process termination within 3 min.')
                            print('=================================')
                        del lSubProcesses[i]  # remove entry from list
                        iChange = 1
                        iFinishedProcessFound = 1
                        break  # because we changed len(lSubProcesses)

                if (iFinishedProcessFound == 0) or (len(lSubProcesses) == 0) or (iExitFlag == 1):
                    break

    finally:
        fKillAllProcesses()
        call('rm -f ' + path.expanduser('~/') + '*.OU',
             shell=True, stdout=open(devnull, "w"))  # delete standard files
        call('rm -f ' + path.expanduser('~/') + '*.ER',
             shell=True, stdout=open(devnull, "w"))  # (on some systems)
        # this has to go through shell


# -------------------------------------------------------------------------------
def DisplayMolgroups(sPath='./'):  # loads parameters and covariance matrix and prints all
    ReflPar.fnPlotMolgroups(sPath)


# -------------------------------------------------------------------------------
def DisplayFit(plotname):
    ReflPar.fnPlotFit(plotname)


# -------------------------------------------------------------------------------
def Result(sPath='./'):  # loads parameters and covariance matrix and prints all
    ReflPar.fnLoadAndPrintPar(sPath)


# -------------------------------------------------------------------------------
def SErr(iiterations, convergence=0.01, sStatMode='is', fConfidence=0.9546):
    try:
        ReflPar.fnLoadAndPrintPar()
        ReflPar.fnBackup()
        filelist = ReflPar.fnLoadFileListAndChangeToLocal()
        ReflPar.fnMake()
    except:
        file = open('StatAborted', 'w')  # indicate that all iterations are done
        file.write('statistical analysis aborted \n')
        file.close()
        return ()
    try:
        if not path.isfile(sStatMode + 'Err.dat'):
            file = open(sStatMode + 'Err.dat', 'w')
            file.write("Chisq " + " ".join(ReflPar.fnGetSortedParNames()) + '\n')
            file.close()
        for i in range(iiterations):
            print('Iteration #', i)
            ReflPar.cGaReflInteractor.fnMCModifyFile(filelist)
            if sStatMode == "is":
                Auto(convergence)  # automatic fitting, independent start
            elif sStatMode == 's':
                AutoFinish2(convergence)  # automatic fitting dependent start
            if path.isfile('pop_bak.dat'):
                call(['cp', 'pop_bak.dat', 'pop.dat'])  # copy newest population into
                # pop.dat, this is done here because
                # sometimes AutoFinish does not find
                # a better fit with the Levenberg-M.
                # and then by default, the generated
                # population is not used, but we want
                # to have it here for the statistical
                # analysis
            ReflPar.fnLoadParameters()
            chisq = ReflPar.fnGetChiSq()
            file = open(sStatMode + 'Err.dat', 'a')
            file.write(str(chisq) + " " + " ".join(ReflPar.fnGetSortedParValues()) + '\n')
            file.close()

        file = open('StatFinished', 'w')  # indicate that all iterations are done
        file.write('statistical analysis finished \n')
        file.close()
    except:
        file = open('StatAborted', 'w')  # indicate that all iterations are done
        file.write('statistical analysis aborted \n')
        file.close()
    finally:
        ReflPar.fnRemoveBackup()  # always restore working dir,
        call(['rm', '-f', '*.mce'])  # also in case of crash
        ReflPar.fnTruncateParDat()


# -------------------------------------------------------------------------------
def StatAnalysis(fConfidence=0.9546, sparse=0):  # Analyzes MC output
    ReflPar.fnAnalyzeStatFile(fConfidence, sparse)


# -------------------------------------------------------------------------------
def StatTable(sTableName, fConfidence=0.9546):  # Produces Latex formatted table
    ReflPar.fnStatTable(sTableName, fConfidence)  # and the actual statistical data


# -------------------------------------------------------------------------------


# main programm from command line
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    if len(argv) == 1:
        print('')
        print('Reflscript usage:')
        print('-----------------------------------------------------------------')
        print('-a [c]               Auto fit')
        print('-aprof               Creates average profiles from ContArray')
        print('                     ContDimZ and ContDimRho files')
        print('-bilayerplot         Creates data for bilayer plots in ')
        print('                     bilayerplotdata.dat')
        print('-conf cf             Confidence for statistical calculatioins')
        print('                     0<=cf<=1, if cf < 0 than cf is interpreted')
        print('                     in units of (negative) sigma')
        print('-cont [z rho]        Create a contour/image plot from sErr.dat')
        print('                     or isErr.dat')
        print('-d [c]               Displacement error analysis')
        print('-env [zGrid sigma]   Calculates nSLD envelopes from stat file')
        print('-f [c]               Finish fit')
        print('-fit filename        Save fit, rho, res in a png file')
        print('                     called filename_fit.png')
        print('-ipc [ipc]           iterations per call/job for -pbs and -m')
        print('-is n [c]            Independent statistical error analysis')
        print('-l [c]               Auto fit finish with correct roughness')
        print('-m [m]               use m parallel processes on one node')
        print('-mol filename        Save molgroups plot in a png file')
        print('                     called filename_mol.png')
        print('-mon                 Monitors a running fit')
        print('-pbs [j workq]       Use pbs to submit j jobs to workq queue')
        print('-pull molgroups      Creates stat profile for a number of molgroups')
        print('-r                   Print results')
        print('-s n [c]             Statistical error analysis with n iterations')
        print('-stat                Analysis of previously calculated stat. data')
        print('                     For n<1 a maximum of 1000 iterations is')
        print('                     carried out until the MC converges and all')
        print('                     parameters do not show a relative change of')
        print('                     more than n compared to the fit interval')
        print('-sparse n            Uses only a subset of n iterations from stat')
        print('                     file. If n<1 than a random chance of f is ')
        print('                     applied to each line that it is used.')
        print('-t                   Truncate par.dat')
        print('')
        print('Displacement Error Analysis:')
        print('For each parameter that should be analyzed, the script expects')
        print('a tag with the syntax"!rstag min max s !" in the line where the')
        print('parameter is initialized with the pars_add command. The tag')
        print('should be embedded in a comment section starting with //. The ')
        print('tag parameters min, max, and s give the relative(!) range and')
        print('stepsize of the displacement of the parameters. For each ')
        print('parameter, an output file is created that contains the fixed ')
        print('parameter value and the achieved chi squared')
        print('')
        print('Statistical Error Analysis')
        print('The script creates n data sets by varying the measured')
        print('reflectivity using random-generated normal deviates')
        print('applied to the uncertainty of the measured data points')
        print('An automated fit determines the fit parameters and they')
        print('are stored in a file called "sErr.dat". A histogram analysis')
        print('should be carried out afterwards with a software like Igor.')
        print('')
        print('Independent Statistical Error Analysis')
        print('This option is equal to "-s" but does not load the stored')
        print('population of a previous fit at start-up. The fit parameters')
        print('are stored in isErr.dat.')
        print('')
        print('Multiprocessor support')
        print('Designed for workstations with multicore architecture. The')
        print('attribute m defines the number of parallel processes used.')
        print('')
        print('PBS batch support')
        print('The j attribute defines the number of jobs submitted once')
        print('at a time. The workq argument is the queue name. PBS and the')
        print('Multi-process support -m can be combined. In this case the')
        print('-ipc option defines the number of MC iterations per process')
        print('The number of iterations per PBS job is then ips * m with')
        print('m being the number of parallel processes.')
        print('')
        print('Contour/Image Plots')
        print('Using SErr.dat or iSErr.dat an array with 3D-data usable')
        print('for contour or image plots is created. The output file is')
        print('ContArr#model.dat. 3D data is generated for all models')
        print('specified in setup.c. For each array two files with axis')
        print('label information are written intended to be used with Igor.')
        print('The attributes z and rho determine the bin size of the')
        print('grid used for mapping of the nSLD profiles.')
        print('')
        print('All fit and error analysis calls may have an additional')
        print('attribute c that sets the convergence condition of the fit')
        print('The fit is viewed at as converged if the change of chi2 after')
        print('a complete alternation over genetic algorigthm, amoeba and')
        print('Levenberg Marquardt is less than c. A c of 0.001 is the')
        print('default')
        print('Example: ./rs.py -iS 1000 0.01 decreases the convergence')
        print('condition by one order of magnitude')
        print('Envelopes')
        print('The option -short startes a less precise but significantly')
        print('faster algorithm. The precision of the envelopes is')
        print('+/- 0.1 sigma. -env takes two arguments, the bin size')
        print('for z and a sigma parameter. Presets are 0.5 in both cases')
        print('The sigma parameter defines the spacing in units of ')
        print('sigma for which envelopes are calculated. A sigma parameter')
        print('of 0 saves all calculates envelopes')
    else:

        bShortFlag = False
        bDisplayMolgroups = False
        bResult = False
        sStatMode = 'none'
        iSummarizeStatistics = 0
        sArchitecture = ''
        sJobSubmission = ''
        sWorkqueue = 'default'
        sJobID = ''
        sPath = './'
        sTableTemplate = 'tabletemplate.tex'
        sMolgroup = ''
        iCm = 0
        iContour = 0
        iContrast = -1
        iEnv = 0
        fZGrid = 0.5
        fSigma = 0.5
        fRhoGrid = 1e-8
        iAutoFinish = 0
        iAutoFinish2 = 0
        iIterationsPerCall = 0
        iNumberOfMCMCIterations = 64000
        iNumberOfParallelThreads = 1
        iNumberOfPBSJobsInQueue = 1
        iNumberOfPbsJobs = 10
        iNumberOfMonteCarloIterations = 1000
        iStatTable = 0
        fConvergence = 0.01
        fConfidence = 0.9546
        fMCConvergence = 0
        fSparse = 0
        iSummarizeStatistic = 0
        iPullMolgroups = 0
        liMolgroups = []
        prefactor = 1

        ReflPar = CMolStat()

        i = 1
        while i < len(argv):
            arg = argv[i]
            if argv[i][0] != '-':
                sPath = argv[i]
            elif argv[i] == '-a':
                if len(argv) > i + 1:
                    fConvergence = float(argv[i + 1])
                Auto(fConvergence)
            elif argv[i] == '-aprof':
                AvgProfile()
            elif argv[i] == '-bilayerplot':
                ReflPar.fnCreateBilayerPlotData()
            elif argv[i] == '-cm':
                iCm = 1
                sMolgroup = argv[i + 1]
                i += 1
            elif argv[i] == '-conf':
                fConfidence = float(argv[i + 1])
                i += 1
            elif argv[i] == '-contrast':
                iContrast = int(argv[i + 1])
                i += 1
            elif argv[i] == '-f':
                if len(argv) > i + 1:
                    fConvergence = float(argv[i + 1])
                AutoFinish(fConvergence)
            elif argv[i] == '-fit':
                import numpy as np
                import matplotlib.pyplot as plt

                if len(argv) > i + 1:
                    plotname = argv[i + 1]
                else:
                    plotname = ''
                DisplayFit(plotname)
            elif argv[i] == '-l':
                if len(argv) > i + 1:
                    fConvergence = float(argv[i + 1])
                iAutoFinish2 = 1
            elif argv[i] == '-mol':
                import numpy as np
                import matplotlib.pyplot as plt

                bDisplayMolgroups = True
            elif argv[i] == '-t':
                ReflPar.fnTruncateParDat()
            elif argv[i] == '-cont':
                if len(argv) > i + 2:
                    fZGrid = float(argv[i + 1])
                    fRhoGrid = float(argv[i + 2])
                    i += 2
                iContour = 1
            elif argv[i] == '-env':
                iEnv = 1
                try:
                    if len(argv) > i + 1:
                        fZGrid = float(argv[i + 1])
                        i += 1
                    if len(argv) > i + 1:
                        fSigma = float(argv[i + 1])
                        i += 1
                except:
                    pass
            elif argv[i] == '-short':
                bShortFlag = True
            elif argv[i] == '-sparse':
                fSparse = float(argv[i + 1])
                i += 1
            elif argv[i] == '-d':
                if len(argv) > i + 1:
                    fConvergence = float(argv[i + 1])
                DErr(fConvergence)
            elif argv[i] == '-is':
                sStatMode = 'is'
                if len(argv) > i + 1:
                    temp = float(argv[i + 1])
                    if temp >= 1:
                        iNumberOfMonteCarloIterations = int(argv[i + 1])
                        fMCConvergence = 0
                    else:
                        iNumberOfMonteCarloIterations = 1000
                        fMCConvergence = temp
                        iIterationsPerCall = 1
                    i += 1
                if len(argv) > i + 1:
                    try:
                        fConvergence = float(argv[i + 1])
                        i += 1
                    except:
                        i = i
            elif argv[i] == '-s':
                sStatMode = 's'
                if len(argv) > i + 1:
                    temp = float(argv[i + 1])
                    if temp >= 1:
                        iNumberOfMonteCarloIterations = int(argv[i + 1])
                        fMCConvergence = 0
                    else:
                        iNumberOfMonteCarloIterations = 1000
                        fMCConvergence = temp
                        iIterationsPerCall = 1
                    i += 1
                if len(argv) > i + 1:
                    try:
                        fConvergence = float(argv[i + 1])
                        i += 1
                    except:
                        i = i
            elif argv[i] == '-m':
                sArchitecture = 'Multicore'
                if len(argv) > i + 1:
                    iNumberOfParallelThreads = int(argv[i + 1])
                    i += 1
            elif argv[i] == '-MCMC':
                if len(argv) > i + 1:
                    iNumberOfMCMCIterations = int(argv[i + 1])
                    i += 1
                fMCMC(iNumberOfMCMCIterations)
            elif argv[i] == '-ipc':
                if len(argv) > i + 1:
                    iIterationsPerCall = int(argv[i + 1])
                    i += 1
            elif argv[i] == '-id':
                if len(argv) > i + 1:
                    sJobID = argv[i + 1]
                    i += 1
            elif argv[i] == '-c':
                if len(argv) > i + 1:
                    fConvergence = float(argv[i + 1])
                    i += 1
            elif argv[i] == '-path':
                sPath = argv[i + 1]
                i += 1
            elif argv[i] == '-pbs':
                sJobSubmission = 'PBS'
                if len(argv) > i + 1:
                    iNumberOfPBSJobsInQueue = int(argv[i + 1])
                    i += 1
                if len(argv) > i + 1:
                    sWorkqueue = argv[i + 1]
                    i += 1
            elif argv[i] == '-pull':
                iPullMolgroups = 1
                liMolgroups = argv[i + 1:]
            elif argv[i] == '-r':
                bResult = True
            elif argv[i] == '-sim':
                i += 1
                mode = 'water'
                if len(argv) > i + 1:
                    if argv[i] == '-mode':
                        mode = argv[i + 1]
                        i += 2
                    if argv[i] == '-pre':
                        prefactor = float(argv[i + 1])
                        i += 2
                ReflPar.fnSimulateReflectivity(mode=mode, pre=prefactor)
                break
            elif argv[i] == '-stat':
                iSummarizeStatistic = 1
            elif argv[i] == '-stattable':
                iStatTable = 1
                sTableTemplate = argv[i + 1]

            i += 1

        if iAutoFinish == 1:
            AutoFinish(fConvergence, sPath)
        elif iAutoFinish2 == 1:
            AutoFinish2(fConvergence, sPath)
        elif iContour == 1:
            if sStatMode == 'none':  # Contour mode needs is based on a MC statistics
                sStatMode = 'is'  # choose default
            fContour(fZGrid, fRhoGrid, sStatMode)
        elif bDisplayMolgroups:
            DisplayMolgroups(sPath)
        elif bResult:
            Result(sPath)
        elif iCm == 1:
            fCalculateMolgroups(fConfidence)
        elif iEnv == 1:
            if sStatMode == 'none':  # Contour mode needs is based on a MC statistics
                sStatMode = 'is'  # choose default
            fEnvelope(fZGrid, fSigma, sStatMode, bShortFlag, iContrast)
        elif sStatMode == 'is' or sStatMode == 's':
            if sJobSubmission == '':
                if sArchitecture == '':
                    SErr(iNumberOfMonteCarloIterations, fConvergence, sStatMode, fConfidence)
                if sArchitecture == 'Multicore':
                    if not iIterationsPerCall:
                        iIterationsPerCall = 10
                    fMCMultiCore(iNumberOfMonteCarloIterations, fMCConvergence, fConvergence, fConfidence, sStatMode,
                                 iNumberOfParallelThreads, iIterationsPerCall, sJobID)
            elif sJobSubmission == 'PBS':
                if not iIterationsPerCall:
                    iIterationsPerCall = 1
                if sArchitecture == '':
                    fMCPBS(iNumberOfMonteCarloIterations, fConvergence, fConfidence, sStatMode,
                           1, iIterationsPerCall, sWorkqueue, iNumberOfPBSJobsInQueue)
                if sArchitecture == 'Multicore':
                    fMCPBS(iNumberOfMonteCarloIterations, fConvergence, fConfidence, sStatMode,
                           iNumberOfParallelThreads, iIterationsPerCall * iNumberOfParallelThreads,
                           sWorkqueue, iNumberOfPBSJobsInQueue)
        elif iSummarizeStatistic == 1:
            StatAnalysis(fConfidence, fSparse)
        elif iStatTable == 1:
            StatTable(sTableTemplate, fConfidence)
        elif iPullMolgroups:
            ReflPar.fnPullMolgroup(liMolgroups, fSparse)
