from __future__ import print_function
from os import getcwd, remove, rename, path, kill, devnull
from random import seed, normalvariate, random
from shutil import rmtree
from sys import argv, exit, stdout
from subprocess import call, Popen
from time import time, gmtime, sleep, ctime
import numpy
import os

from molgroups.support import molstat

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

    while True:  # genetic runs until chisq<20 or not improving
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


def AutoFinish(convergence=0.001, sPath='./'):
    if sPath != './':
        ReflPar.fnImportMCMCBestFit(sPath)
        # automatic fit starting with Amoeba, approximate roughness
    ReflPar.fnLoadParameters()
    fOldChiSq = ReflPar.fnGetChiSq()

    while True:  # Amoeba, approximate roughness
        call(['nice', './fit', '-peaS'], stdout=open(devnull, "w"))  # until chisq<10 or no improvement
        print('Amoeba, approximate roughness')
        ReflPar.fnLoadAndPrintPar()
        call(['cp', 'pop_bak.dat', 'pop.dat'])  # copy newest population into pop.dat
        if ReflPar.chisq < 5 or (fOldChiSq - ReflPar.fnGetChiSq() < convergence):
            break
        fOldChiSq = ReflPar.fnGetChiSq()
    AutoFinish2(convergence)  # AutoFinish2 takes over


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

    while True:
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
        while True:
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
        while True:
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


def AvgProfile():
    iCounter = 0
    while path.isfile(f"ContDimRho{iCounter}.dat"):
        with open(f"ContDimRho{iCounter}.dat", 'r') as file:
            data = file.readlines()
        rho = list(map(float, data))
        
        with open(f"ContDimZ{iCounter}.dat", 'r') as file:
            data = file.readlines()
        z = list(map(float, data))
        z = [(z[i] + z[i + 1]) / 2 for i in range(len(z) - 1)] #correct for Igor-style array axis notation

        with open(f"ContArray{iCounter}.dat", 'r') as file:
            data = file.readlines()
        arr = [list(map(float, rhoarr.split())) for rhoarr in data]

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

        with open(f"AvgProfile{iCounter}.dat", "w") as file:
            file.write('z    maxlikely    median    lowperc     highperc \n')
            for i in range(len(z)):
                file.write(f"{z[i]} {maxlikely[i]} {median[i]} {lowperc[i]} {highperc[i]}\n")
        iCounter += 1

def fContour(dZGrid=0.5, dRhoGrid=1e-8, sStatMode='is'):  # Calculate contour plot data from SErr.dat
    ReflPar.fnContourData(sStatMode + 'Err.dat', dZGrid, dRhoGrid)


def fCalculateMolgroups(fConfidence):
    ReflPar.fnCalculateMolgroupProperty(fConfidence)


def fnDetermineFitSoftware():
    if path.isfile('run.py'):
        print('Refl1D setup identified.')
        return 'refl1d'
    else:
        print('Garefl setup identified.')
        return 'garefl'


def fnDetermineStatFile():
    if path.isfile('isErr.dat'):
        return 'isErr.dat'
    elif path.isfile('sErr.dat'):
        return 'sErr.dat'
    else:
        return ''


def DErr(convergence=0.001):  # function has not yet been thoroughly tested
    # Error analysis by stepwise parameter displacement,
    # fixing, and fitting within a defined range
    def fnWriteToFile(sParameterName, fTestValue, fChiSq):  # append new data set to file
        try:
            with open(f"Error_{sParameterName}.dat", "r") as file:  # if file exists, open and read data
                data = file.readlines()
        except IOError:
            data = []  # otherwise start with an empty data file
        newdata = data[:]  # copy dat into new object
        newdata.append(f"{fTestValue}  {fChiSq}\n")  # append line with the data parameter to store
        with open(f"Error_{sParameterName}.dat", "w") as file:  # create new file and write out the new data
            file.writelines(newdata)

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


def fEnvelope(fZGrid=0.1, fSigma=0, sStatMode='is', bShortFlag=False,
              iContrast=-1):  # Calculate nSLD envelopes plot data from SErr.dat
    ReflPar.fnnSLDEnvelopes(fZGrid, fSigma, sStatMode + 'Err.dat', bShortFlag, iContrast)


def fMCMultiCore(iIterations=1000, fMCConvergence=0.01, iConvergence=0.01,
                 fConfidence=0.9546, sMode='is', iNodes=1, iIterationsPerCall='10',
                 sJobID=''):  # Monte Carlo for Multicore systems

    def fCleanUpDirectory(sDir):
        call(['rm', '-rf', sDir])  # remove rs directory

    def fCreateRsDir(sMode, lSubProcesses):
        while True:  # create random name directory
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

            while True:
                iFinishedProcessFound = 0

                for i in range(len(lSubProcesses)):  # check for finished sub processes
                    sDirName = lSubProcesses[i][1]
                    pid = lSubProcesses[i][0]

                    if path.isfile(os.path.join(sDirName, 'StatFinished')):  # look up if process is finished
                        if path.isfile(sMode + 'Err.dat'):
                            with open(os.path.join(sDirName, sMode) + 'Err.dat', 'r') as file:  # get statistical data from rs subdir
                                data = file.readlines()
                            data = data[1:]  # delete headerline
                            iFailureCounter = 0
                            while True:
                                try:
                                    with open(sMode + 'Err.dat', 'a') as file:  # append new data to file in root dir
                                        file.writelines(data)
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

                    if path.isfile(os.path.join(sDirName, 'StatAborted')):  # look up if process is finished
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
            with open(sJobID, 'w') as file: # indicate Multicore is finished if required
                file.write('Multicore finished \n')  # requirement comes from presence of a
                                                     # job id as given by a PBS caller

    finally:
        print('Exiting MC')
        print(iExitFlag, iDoneIterations, iIterations)
        fKillAllProcesses()


def fMCMC(iMaxIterations=1024000, liMolgroups=['protein'], fSparse=0):
    while True:
        if not path.isfile('run.py'):  # make sure there is a run.py
            with open('run.py', 'w') as file:
                file.write('from bumps.fitproblem import FitProblem\n')
                file.write('from refl1d import garefl\n')
                file.write('from refl1d.names import *\n')
                file.write('\n')
                file.write("problem = garefl.load('model.so')\n")
                file.write('\n')
                file.write('problem.penalty_limit = 50\n')

        iBurn = 4000  # check whether an existing MCMC exists
        bMCMCexists = False
        for i in range(1, 9):
            iBurn = iBurn * 2
            if path.isdir(f'MCMC_{iBurn}_500'):
                print(f'Found MCMC_{iBurn}_500 \n')
                bMCMCexists = True
                break

        lCommand = ['refl1d_cli.py', 'run.py', '--fit=dream', '--parallel=0']
        if bMCMCexists:
            if iBurn >= iMaxIterations:
                print('Maximum number of MCMC iterations reached\n')
                break  # end
            lCommand.append(f'--resume=MCMC_{iBurn}_500')
            lCommand.append(f'--store=MCMC_{iBurn*2}_500')
            lCommand.append(f'--burn={iBurn}')
        else:
            lCommand.append('--init=lhs')
            lCommand.append('--store=MCMC_8000_500')
            lCommand.append('--burn=8000')
            iBurn = 4000

        lCommand.append('--steps=500')

        call(lCommand)  # run MCMC

        rename(f'MCMC_{iBurn*2}_500', 'MCMC')  # create sErr.dat
        if path.isfile('isErr.dat'):
            remove('isErr.dat')
        if path.isfile('sErr.dat'):
            remove('sErr.dat')
        StatAnalysis(-1, 0.005)  # sErr.dat contains about 1000 iterations
        rename('MCMC', f'MCMC_{iBurn*2}_500')

        if liMolgroups:
            if path.isfile('StatDataPython.dat'):
                remove('StatDataPython.dat')
            ReflPar.fnPullMolgroup(liMolgroups, 0)

        if bMCMCexists:
            rmtree(f'MCMC_{iBurn}_500')

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
            while True:
                if path.isfile(sOutName) and path.isfile(sErrName):  # wait for output files
                    call(['rm', '-f', sOutName])
                    call(['rm', '-f', sErrName])
                    break
                retcode = call('ls ' + path.expanduser('~/') + '*.OU',
                               shell=True, stdout=open(devnull, "w"))
                if retcode == 0:  # on some systems only
                    break  # standardfiles are written
                sleep(1)
                iSleepCounter += 1
                if iSleepCounter > 20:
                    print('Waited 20s for output files to be written ... giving up.')
                    break

    def fCreateMultiCoreJobID(lSubProcesses):
        while True:  # create random name directory
            sRandomName = str(int(random() * 1000000000000))  # and check if the name already
            iNameExists = False  # as a PBS Job ID
            for element in lSubProcesses:
                if sRandomName == element[1]:
                    iNameExists = True
            if not iNameExists:
                return sRandomName # randomname will be job id

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
                data = ['#PBS -l ncpus=1\n', '#PBS -l cput=48:00:00\n', '#PBS -l walltime=72:00:00\n',
                        f'#PBS -e {getcwd()}/{sJobID}.err\n',
                        f'#PBS -o {getcwd()}/{sJobID}.out\n', 'cd $PBS_O_WORKDIR\n',
                        f'./rs.py -{sMode} {iIterationsPerCall} -c {iConvergence} -m {iNodes} -ipc ' \
                            f'{iIterationsPerCall / iNodes} -id {sJobID}\n',
                        '#end']  # create run batchfile
                # data.append('#PBS -l nodes=1:ppn=1\n')
                runfile = 'run.sh'
                with open(runfile, "w") as file:
                    file.writelines(data)
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

            while True:
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
        call(f"rm -f {path.expanduser('~/')}*.OU",
             shell=True, stdout=open(devnull, "w"))  # delete standard files
        call(f"rm -f {path.expanduser('~/')}*.ER",
             shell=True, stdout=open(devnull, "w"))  # (on some systems)
        # this has to go through shell


def DisplayMolgroups(sPath='./'):  # loads parameters and covariance matrix and prints all
    ReflPar.fnPlotMolgroups(sPath)


def DisplayFit(plotname):
    ReflPar.fnPlotFit(plotname)


def Result(sPath='./'):  # loads parameters and covariance matrix and prints all
    ReflPar.fnLoadAndPrintPar(sPath)


def SErr(iiterations, convergence=0.01, sStatMode='is', fConfidence=0.9546):
    try:
        ReflPar.fnLoadAndPrintPar()
        ReflPar.fnBackup()
        filelist = ReflPar.fnLoadFileListAndChangeToLocal()
        ReflPar.fnMake()
    except:
        with open('StatAborted', 'w') as file: # indicate that all iterations are done
            file.write('statistical analysis aborted \n')
            return ()
    try:
        if not path.isfile(sStatMode + 'Err.dat'):
            with open(sStatMode + 'Err.dat', 'w') as file:
                file.write("Chisq " + " ".join(ReflPar.fnGetSortedParNames()) + '\n')
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
            with  open(sStatMode + 'Err.dat', 'a') as file:
                file.write(str(chisq) + " " + " ".join(ReflPar.fnGetSortedParValues()) + '\n')

        with open('StatFinished', 'w') as file:  # indicate that all iterations are done
            file.write('statistical analysis finished \n')
    except:
        with open('StatAborted', 'w') as file: # indicate that all iterations are done
            file.write('statistical analysis aborted \n')
    finally:
        ReflPar.fnRemoveBackup()  # always restore working dir,
        call(['rm', '-f', '*.mce'])  # also in case of crash
        ReflPar.fnTruncateParDat()


def StatAnalysis(fConfidence=0.9546, sparse=0):  # Analyzes MC output
    ReflPar.fnAnalyzeStatFile(fConfidence, sparse)


def StatTable(sTableName, fConfidence=0.9546):  # Produces Latex formatted table
    ReflPar.fnStatTable(sTableName, fConfidence)  # and the actual statistical data


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

        ReflPar = molstat.CMolStat()

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
