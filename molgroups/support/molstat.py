from __future__ import print_function
from math import sqrt
from os import path
from scipy import stats, special
import numpy
import pandas
import os
import bumps.curve


class CMolStat:
    def __init__(self, fitsource="refl1d", spath=".", mcmcpath=".",
                 runfile="run", state=None, problem=None,
                 load_state=True, save_stat_data=False):
        """
        self.diParameters is a dictionary containing all parameters. Structure:

        dictionary: sparameter : dictionary
                                 'number'    : int       # order of initialization in ga_refl
                                 'lowerlimit'  : float   # constraint range lower limit
                                 'upperlimit'  : float   # constraint range upper limit
                                 'value'  : float        # absolute value
                                 'error'  : float        # error derived from covariance matrix
                                 'relval' : float        # relative value between 0 and 1 in terms of the constraints
                                 'variable': string       # associated ga_refl variable

        self.diStatRestuls is a dictionary of statistical results with entries from various routines:
            (from fnAnalyzeStatFile)
            'Parameters' is itself a dictionary with the following keys:
                      'Values' key contains ordered list of all MCMC parameters
                      'LowPerc' contains lower percentile
                      'Median' contains median
                      'HighPerc' contains higher percentile
                      'MaxParameterLength' contains length of the longest parameter name

            (from fnLoadStatData)
            'NumberOfStatValues' number of MC points
            'nSLDProfiles' contains all profiles 'nSLDProfiles'[MCiteration][model][z,rho]
            'Molgroups' contains a list of dictionaries storing all molgroups
                'Molgroups'[MCiteration] {'molgroupname' {'zaxis' [z]
                                                          'are' [area]
                                                          'sl' [sl]
                                                          'sld' [sld]
                                                          'headerdata' {'property' value}
                                                          'property1' value
                                                          'property2' value}}
           'Results' contains a dictionary of derived values from fit paramters for post-analysis
                'Results' is itself a dictionary with the following keys:
                        'Values' key contains ordered list of all MCMC derived parameters
        """

        self.diParameters = {}
        self.diMolgroups = {}
        self.diStatResults = {}
        self.diResults = {}

        self.fitsource = fitsource  # define fitting software
        self.spath = spath
        self.mcmcpath = mcmcpath
        self.runfile = runfile

        self.liStatResult = []  # a list that contains the isErr.dat or sErr.dat file line by line
        self.sStatResultHeader = ''  # headerline from isErr.dat or sErr.dat
        self.sStatFileName = ''  # Name of the statistical File
        self.chisq = 0.
        self.fMolgroupsStepsize = 0.
        self.iMolgroupsDimension = 0
        self.fMolgroupsNormArea = 0
        # check for system and type of setup file

        self.Interactor = None
        if self.fitsource == "bumps":
            from molgroups.support import api_bumps
            self.Interactor = api_bumps.CBumpsAPI(spath, mcmcpath, runfile, state, problem, load_state=load_state)
        elif self.fitsource == 'refl1d':
            from molgroups.support import api_refl1d
            self.Interactor = api_refl1d.CRefl1DAPI(spath, mcmcpath, runfile, load_state=load_state)
        elif self.fitsource == 'garefl':
            from molgroups.support import api_garefl
            self.Interactor = api_garefl.CGaReflAPI(spath, mcmcpath, runfile, load_state=load_state)
        elif self.fitsource == 'SASView':
            from molgroups.support import api_sasview
            self.Interactor = api_sasview.CSASViewAPI(spath, mcmcpath, runfile, load_state=load_state)

        self.save_stat_data = save_stat_data

    def fnAnalyzeStatFile(self, fConfidence=-1, sparse=0):

        def data_append(data, origin, name, vis, lower_limit, upper_limit, lower_percentile, median_percentile,
                        upper_percentile, interval_lower, interval_upper, confidence):
            data['origin'].append(origin)
            data['name'].append(name)
            data['vis'].append(vis)
            data['lower limit'].append(lower_limit)
            data['upper limit'].append(upper_limit)
            data['lower percentile'].append(lower_percentile)
            data['median percentile'].append(median_percentile)
            data['upper percentile'].append(upper_percentile)
            data['interval lower'].append(interval_lower)
            data['interval upper'].append(interval_upper)
            data['confidence'].append(confidence)
            return data

        # self.fnLoadParameters()
        self.fnLoadStatData(sparse)

        data = {'origin': [],
                'name': [],
                'vis': [],
                'lower limit': [],
                'upper limit': [],
                'lower percentile': [],
                'median percentile': [],
                'upper percentile': [],
                'interval lower': [],
                'interval upper': [],
                'confidence': []
                }

        fConfidence = min(fConfidence, 1)
        if fConfidence < 0:
            fConfidence = special.erf(-1 * fConfidence / sqrt(2))
        percentiles = (100.0 * (1 - fConfidence) / 2, 50., 100.0 - 100.0 * (1 - fConfidence) / 2)
        iNumberOfMCIterations = self.diStatResults['NumberOfStatValues']

        print('Analysis of MCMC fit ...')
        print('Number of iterations: %(ni)d' % {'ni': iNumberOfMCIterations})
        print('')
        print('Fit Parameters:')

        # Fit Parameters
        for element in sorted(list(self.diParameters.keys()),
                              key=lambda sParameter: self.diParameters[sParameter]['number']):
            vals = self.diStatResults['Parameters'][element]['Values']
            perc = stats.scoreatpercentile(vals, percentiles)
            self.diStatResults['Parameters'][element]['LowPerc'] = perc[0]
            self.diStatResults['Parameters'][element]['Median'] = perc[1]
            self.diStatResults['Parameters'][element]['HighPerc'] = perc[2]

            flowerlimit = self.diParameters[element]['lowerlimit']
            fupperlimit = self.diParameters[element]['upperlimit']
            temp = abs(fupperlimit - flowerlimit)

            sGraphOutput = '['
            itemp1 = int((perc[0] - flowerlimit) / temp * 10 + 0.5)
            itemp2 = int((perc[1] - flowerlimit) / temp * 10 + 0.5)
            itemp3 = int((perc[2] - flowerlimit) / temp * 10 + 0.5)
            for i in range(11):
                s1 = ' '
                if itemp1 == i or itemp3 == i:
                    s1 = '|'
                if itemp2 == i:
                    if s1 == '|':
                        s1 = '+'
                    else:
                        s1 = '-'
                sGraphOutput += s1
            sGraphOutput += ']'
            if (perc[0] - flowerlimit) < temp * 0.01:
                self.diStatResults['Parameters'][element]['LowerLimitCollision'] = True
                sGraphOutput = '#' + sGraphOutput[1:]
            else:
                self.diStatResults['Parameters'][element]['LowerLimitCollision'] = False
            if (fupperlimit - perc[2]) < temp * 0.01:
                self.diStatResults['Parameters'][element]['UpperLimitCollision'] = True
                sGraphOutput = sGraphOutput[:-1] + '#'
            else:
                self.diStatResults['Parameters'][element]['UpperLimitCollision'] = False

            data = data_append(data, 'fit', element, sGraphOutput, flowerlimit, fupperlimit, perc[0], perc[1],
                               perc[2], perc[0]-perc[1], perc[2]-perc[1], fConfidence)

        # Derived parameters – results
        for origin in self.diStatResults['Results']:
            for name in self.diStatResults['Results'][origin]:
                vals = self.diStatResults['Results'][origin][name]
                perc = stats.scoreatpercentile(vals, percentiles)
                data = data_append(data, origin, name, '', None, None, perc[0], perc[1], perc[2], perc[0]-perc[1],
                                   perc[2]-perc[1], fConfidence)

        # Pandas Dataframe
        return pandas.DataFrame(data)

    @staticmethod
    def fnCalcConfidenceLimits(data, method=1):
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
            if numpy.isscalar(ci):
                ci = [ci]

            # Simple solution: ci*N is the number of points in the interval, so
            # find the width of every interval of that size and return the smallest.
            result = [_find_interval(x, i) for i in ci]

            if len(ci) == 1:
                result = result[0]
            return result

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
            a = c = maxindex
            if 0 < maxindex < len(histo) - 1:
                a = maxindex - 1
                c = maxindex + 1
            maxindexsmooth = a * histo[a] + maxindex * histo[maxindex] + c * histo[c]
            maxindexsmooth = maxindexsmooth / (histo[a] + histo[maxindex] + histo[c])
            maxvaluesmooth = low_range + (maxindexsmooth + 0.5) * bin_size

            a = c = maxindex
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

            return [twosigmam, onesigmam, maxvaluesmooth, onesigmap, twosigmap]

        # shortest confidence interval method, NIST recommended
        else:
            twosigmam, twosigmap = credible_interval(data, 0.95)
            onesigmam, onesigmap = credible_interval(data, 0.68)
            reported = 0.5 * (onesigmam + onesigmap)
            return [twosigmam, onesigmam, reported, onesigmap, twosigmap]

    def fnGetChiSq(self):  # export chi squared
        return self.chisq

    def fnGetParameterValue(self, sname):  # export absolute parameter value
        return self.diParameters[sname]['value']  # for given name

    def fnGetSortedParNames(self):  # return a list of sorted parameter
        litest = list(self.diParameters.keys())
        litest = sorted(litest, key=lambda keyitem: self.diParameters[keyitem]['number'])
        return litest

    def fnGetSortedParValues(self):  # the same as above but it returns
        litest = list(self.diParameters.keys())  # a list of sorted parameter values
        litest = sorted(litest, key=lambda keyitem: self.diParameters[keyitem]['number'])
        lvalue = []
        for parameter in litest:
            lvalue.append(str(self.diParameters[parameter]['value']))
        return lvalue

    def fnLoadAndPrintPar(self, sPath='./'):
        self.fnLoadParameters()
        self.fnLoadCovar(sPath)
        self.fnPrintPar()

    def fnLoadCovar(self, sPath):
        """
        loads a covariance matrix, calculates the errors and stores it into self.diParameters
        The function fnLoadParameters must have been already carried out
        """
        if path.isfile(f'{sPath}covar.dat'):
            with open(f'{sPath}covar.dat') as file:
                data = file.readlines()
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

    @staticmethod
    def fnLoadObject(sFileName):
        import pickle
        with open(sFileName, "rb") as file:
            Object = pickle.load(file)
        return Object

    def fnLoadParameters(self):
        if self.diParameters == {}:
            self.diParameters, self.chisq = self.Interactor.fnLoadParameters()

    def fnLoadStatData(self, sparse=0):
        self.fnLoadParameters()
        if self.diStatResults != {}:
            return

        try:
            self.diStatResults = self.fnLoadObject(os.path.join(self.spath, self.mcmcpath, 'StatDataPython.dat'))
            print('Loaded statistical data from StatDataPython.dat')
            return
        except IOError:
            print('No StatDataPython.dat.')
            print('Recreate statistical data from sErr.dat.')

        self.diStatResults = self.Interactor.fnLoadStatData(sparse)
        # cycle through all parameters, determine length of the longest parameter name for displaying
        iMaxParameterNameLength = 0
        for parname in list(self.diStatResults['Parameters'].keys()):
            if len(parname) > iMaxParameterNameLength:
                iMaxParameterNameLength = len(parname)
        self.diStatResults['MaxParameterLength'] = iMaxParameterNameLength
        self.diStatResults['NumberOfStatValues'] = \
            len(self.diStatResults['Parameters'][list(self.diStatResults['Parameters'].keys())[0]]['Values'])

        # Recreates profile and fit data associated with parameter stats
        if self.Interactor.problem is not None:
            problem = self.Interactor.problem
        else:
            problem = self.Interactor.fnRestoreFitProblem()

        j = 0
        self.diStatResults['nSLDProfiles'] = []
        self.diStatResults['Molgroups'] = {}
        self.diStatResults['Results'] = {}

        for iteration in range(self.diStatResults['NumberOfStatValues']):
            try:
                # appends a new list for profiles for the current MC iteration
                self.diStatResults['nSLDProfiles'].append([])
                liParameters = list(self.diParameters.keys())

                # sort by number of appereance in setup file
                liParameters = sorted(liParameters, key=lambda keyitem: self.diParameters[keyitem]['number'])
                bConsistency = set(liParameters).issubset(self.diStatResults['Parameters'].keys())
                if not bConsistency:
                    raise RuntimeError('Statistical error data and setup file do not match')

                p = []
                for parameter in liParameters:
                    val = self.diStatResults['Parameters'][parameter]['Values'][iteration]
                    p.append(val)
                problem.setp(p)
                problem.model_update()

                # TODO: By calling .chisq() I currently force an update of the BLM function. There must be a better
                #   way, also implement which contrast to use for pulling groups garefl based code should save a
                #   mol.dat automatically on updating the model
                if 'models' in dir(problem):
                    for M in problem.models:
                        M.chisq()
                        break
                else:
                    problem.chisq()

                # distinguish between FitProblem and MultiFitProblem
                if 'models' in dir(problem):
                    for M in problem.models:
                        if not isinstance(M.fitness, bumps.curve.Curve):
                            z, rho, irho = self.Interactor.fnRestoreSmoothProfile(M)
                            self.diStatResults['nSLDProfiles'][-1].append((z, rho, irho))
                else:
                    z, rho, irho = self.Interactor.fnRestoreSmoothProfile(problem)
                    self.diStatResults['nSLDProfiles'][-1].append((z, rho, irho))

                # Recreate Molgroups and Derived Results
                self.diMolgroups, self.diResults = self.Interactor.fnLoadMolgroups(problem)

                # Store Molgroups
                # self.diStatResults['Molgroups'].append(self.diMolgroups)
                for name in self.diMolgroups:
                    if name not in self.diStatResults['Molgroups']:
                        self.diStatResults['Molgroups'][name] = {}
                    for entry in self.diMolgroups[name]:
                        if entry == 'zaxis':
                            # store z-axis only once
                            if entry not in self.diStatResults['Molgroups'][name]:
                                self.diStatResults['Molgroups'][name][entry] = self.diMolgroups[name][entry]
                        else:
                            if entry not in self.diStatResults['Molgroups'][name]:
                                self.diStatResults['Molgroups'][name][entry] = []
                            self.diStatResults['Molgroups'][name][entry].append(self.diMolgroups[name][entry])

                # Store Derived Results
                # origin is the name of the object that provided a result with a certain name
                for origin in self.diResults:
                    if origin not in self.diStatResults['Results']:
                        self.diStatResults['Results'][origin] = {}
                    for name in self.diResults[origin]:
                        if name not in self.diStatResults['Results'][origin]:
                            self.diStatResults['Results'][origin][name] = []
                        self.diStatResults['Results'][origin][name].append(self.diResults[origin][name])

            finally:
                j += 1

        # save stat data to disk, if flag is set
        if self.save_stat_data:
            self.fnSaveObject(self.diStatResults, f'{self.spath}/{self.mcmcpath}/StatDataPython.dat')

    def fnPrintPar(self):
        # prints parameters and their errors from the covariance matrix onto the screen

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

    def fnProfilesStat(self, sparse=0):
        self.fnLoadStatData(sparse)

        for group in self.diStatResults['Molgroups']:
            median_area = numpy.percentile(self.diStatResults['Molgroups'][group]['area'], 50., axis=0)
            psigma_area = numpy.percentile(self.diStatResults['Molgroups'][group]['area'], 82., axis=0)
            msigma_area = numpy.percentile(self.diStatResults['Molgroups'][group]['area'], 18., axis=0)
            median_sl = numpy.percentile(self.diStatResults['Molgroups'][group]['sl'], 50., axis=0)
            psigma_sl = numpy.percentile(self.diStatResults['Molgroups'][group]['sl'], 82., axis=0)
            msigma_sl = numpy.percentile(self.diStatResults['Molgroups'][group]['sl'], 18., axis=0)
            median_sld = numpy.percentile(self.diStatResults['Molgroups'][group]['sld'], 50., axis=0)
            psigma_sld = numpy.percentile(self.diStatResults['Molgroups'][group]['sld'], 82., axis=0)
            msigma_sld = numpy.percentile(self.diStatResults['Molgroups'][group]['sld'], 18., axis=0)

            self.diStatResults['Molgroups'][group]['median area'] = median_area
            self.diStatResults['Molgroups'][group]['psigma area'] = psigma_area
            self.diStatResults['Molgroups'][group]['msigma area'] = msigma_area
            self.diStatResults['Molgroups'][group]['median sl'] = median_sl
            self.diStatResults['Molgroups'][group]['psigma sl'] = psigma_sl
            self.diStatResults['Molgroups'][group]['msigma sl'] = msigma_sl
            self.diStatResults['Molgroups'][group]['median sld'] = median_sld
            self.diStatResults['Molgroups'][group]['psigma sld'] = psigma_sld
            self.diStatResults['Molgroups'][group]['msigma sld'] = msigma_sld

    def fnPullMolgroup(self, liMolgroupNames, sparse=0, verbose=True):
        """
        Calls Function that recreates statistical data and extracts only area and nSL profiles for
        submolecular groups whose names are given in liMolgroupNames. Those groups
        are added for each iteration and a file pulledmolgroups.dat is created.
        A statistical analysis area profile containing the median, sigma, and
        2 sigma intervals are put out in pulledmolgroupsstat.dat.
        Save results to file
        """
        diarea, dinsl, dinsld = self.fnPullMolgroupLoader(liMolgroupNames, sparse, verbose=verbose)
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

        self.Interactor.fnSaveSingleColumns(self.mcmcpath + '/pulledmolgroups_area.dat', diarea)
        self.Interactor.fnSaveSingleColumns(self.mcmcpath + '/pulledmolgroups_nsl.dat', dinsl)
        self.Interactor.fnSaveSingleColumns(self.mcmcpath + '/pulledmolgroups_nsld.dat', dinsld)
        self.Interactor.fnSaveSingleColumns(self.mcmcpath + '/pulledmolgroupsstat.dat', diStat)

    def fnPullMolgroupLoader(self, liMolgroupNames, sparse=0, verbose=True):
        """
        Function recreates statistical data and extracts only area and nSL profiles for
        submolecular groups whose names are given in liMolgroupNames. Those groups
        are added for each iteration and a file pulledmolgroups.dat is created.
        A statistical analysis area profile containing the median, sigma, and
        2 sigma intervals are put out in pulledmolgroupsstat.dat.
        """

        if self.diStatResults == {}:
            self.fnLoadStatData()

        zaxis = self.diStatResults['Molgroups'][0][list(self.diStatResults['Molgroups'][0].keys())[0]]['zaxis']
        diarea = {'zaxis': zaxis}
        dinsl = {'zaxis': zaxis}
        dinsld = {'zaxis': zaxis}

        for i, iteration in enumerate(self.diStatResults['Molgroups']):
            # add together all the molgroups that have to be analyzed
            sumareaprofile = numpy.zeros_like(zaxis)
            sumnslprofile = numpy.zeros_like(zaxis)

            for molgroup in liMolgroupNames:
                if molgroup in iteration:
                    sumareaprofile += iteration[molgroup]['areaaxis']
                    sumnslprofile += iteration[molgroup]['nslaxis']
                elif i == 0 and verbose:
                    print(f'Molecular group {molgroup} does not exist.')

                diarea[f'iter{i}'] = sumareaprofile
                dinsl[f'iter{i}'] = sumnslprofile
                stepsize = diarea['zaxis'][1] - diarea['zaxis'][0]
                dinsld[f'iter{i}'] = []
                for j in range(len(diarea[f'iter{i}'])):
                    if diarea[f'iter{i}'][j] != 0:
                        dinsld[f'iter{i}'].append(dinsl[f'iter{i}'][j] / diarea[f'iter{i}'][j] / stepsize)
                    else:
                        dinsld[f'iter{i}'].append(0.0)

        return diarea, dinsl, dinsld

    def fnRestoreFit(self):
        self.Interactor.fnRestoreFit()

    def fnRunFit(self, burn=2000, steps=500, batch=False):
        path1 = os.path.join(self.spath, self.mcmcpath)
        if os.path.isfile(os.path.join(path1, "sErr.dat")):
            os.remove(os.path.join(path1, "sErr.dat"))
        if os.path.isfile(os.path.join(path1,  "isErr.dat")):
            os.remove(os.path.join(path1,  "isErr.dat"))
        if os.path.isfile(os.path.join(path1,  "StatDataPython.dat")):
            os.remove(os.path.join(path1,  "StatDataPython.dat"))
        self.Interactor.fnRunMCMC(burn, steps, batch=False)

    @staticmethod
    def fnSaveObject(save_object, sFileName):
        import pickle

        with open(sFileName, "wb") as file:
            pickle.dump(save_object, file)

    def fnSimulateData(self, basefilename='sim.dat', liConfigurations=None, qmin=None, qmax=None, qrangefromfile=False,
                       t_total=None, mode='water', lambda_min=0.1, verbose=True, simpar=None, save_file=True):
        """
        Simulates scattering based on a parameter file called simpar.dat
        requires a ready-to-go fit whose fit parameters are modified and fixed
        The basename can refer to a set of data files with integer indizes before the suffix
        """

        # Load Parameters
        self.diParameters, _ = self.Interactor.fnLoadParameters()

        liParameters = list(self.diParameters.keys())
        liParameters = sorted(liParameters, key=lambda keyitem: self.diParameters[keyitem]['number'])

        if simpar is None:
            # the file simpar.dat contains the parameter values to be simulated
            # this could be done fileless
            simpar = pandas.read_csv(self.spath + '/simpar.dat', sep='\s+', header=None, names=['par', 'value'],
                                     skip_blank_lines=True, comment='#')
        else:
            # simpar is provided as a dataframe with the appropriate shape
            simpar.columns = ['par', 'value']

        if verbose:
            print(simpar)
            print(liConfigurations)

        diModelPars = {}
        for parameter in liParameters:
            diModelPars[parameter] = simpar[simpar.par == parameter].iloc[0][1]
        # load all data files into a list of Pandas dataframes
        # each element is itself a list of [comments, simdata]
        liData = self.Interactor.fnLoadData(basefilename)
        liData = self.Interactor.fnSimulateDataPlusErrorBars(liData, diModelPars, simpar=simpar,
                                                             basefilename=basefilename,
                                                             liConfigurations=liConfigurations, qmin=qmin,
                                                             qmax=qmax, qrangefromfile=qrangefromfile,
                                                             lambda_min=lambda_min, mode=mode, t_total=t_total)
        if save_file:
            self.Interactor.fnSaveData(basefilename, liData)
        else:
            return liData

