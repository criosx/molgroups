from __future__ import print_function
from math import fabs
from os import path
from re import VERBOSE, IGNORECASE, compile
from subprocess import call, Popen
from time import sleep
import pandas
import shutil
import glob
import os

from molgroups.support import api_refl1d


class CGaReflAPI(api_refl1d.CRefl1DAPI):
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
            super(api_refl1d.CRefl1DAPI, self).fnLoadStatData(dSparse)
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
        api_refl1d.CRefl1DAPI.fnRunMCMC(self, burn=burn, steps=steps, batch=batch)

    def fnSaveMolgroups(self, problem):
        # should call the setup.cc function that saves mol.dat
        # TODO: Needs testing
        problem.active_model.fitness.output_model()

    def fnSimulateData(self, diAddition, liData, data_column='R'):
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
            liData[i][1][data_column] = simdata['fit']
            i += 1

        return liData

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
