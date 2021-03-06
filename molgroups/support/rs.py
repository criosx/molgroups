from __future__ import print_function
from math import fabs, pow, floor, ceil, sqrt, log10
from collections import defaultdict
from numpy import subtract, minimum, maximum, average, array
from operator import itemgetter
from os import getcwd, remove, rename, path, kill, devnull
from random import seed, normalvariate, random
from re import VERBOSE, IGNORECASE, compile
from scipy import stats, special
from shutil import rmtree
from sys import argv, exit, stdout
from subprocess import call, Popen
from time import time, gmtime, sleep, ctime
import numpy
import pandas

from . import general
from . import rsdi


class CMolStat:
    def __init__(self, fitsource="refl1d", spath=".", mcmcpath=".",
                 runfile="run", state=None, problem=None,
                 load_state=True, save_stat_data=False):
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
        if self.fitsource == "bumps":
            self.Interactor = rsdi.CBumpsInteractor(spath, mcmcpath, runfile, state, problem, load_state=load_state)
        elif self.fitsource == 'refl1d':
            self.Interactor = rsdi.CRefl1DInteractor(spath, mcmcpath, runfile, load_state=load_state)
        elif self.fitsource == 'garefl':
            self.Interactor = rsdi.CGaReflInteractor(spath, mcmcpath, runfile, load_state=load_state)
        elif self.fitsource == 'SASView':
            self.Interactor = rsdi.CSASViewInteractor(spath, mcmcpath, runfile, load_state=load_state)

        self.save_stat_data = save_stat_data

    def fnAnalyzeStatFile(self, fConfidence=-1, sparse=0):  # summarizes stat file

        self.fnLoadParameters()                         # Load Parameters for limits
        self.fnLoadStatData(sparse)                     # Load data from file into list

        fConfidence = min(fConfidence, 1)
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
        print('Confidence level: {fConfidence:.4f}')

    def fnCalculateMolgroupProperty(self, fConfidence, verbose=True):
        def fnFindMaxFWHM(lList):
            maxvalue = max(lList)
            imax = lList.index(maxvalue)
            ifwhmplus = ifwhmminus = imax
            while lList[ifwhmplus] > (maxvalue / 2) and ifwhmplus < (len(lList) - 1):
                ifwhmplus += 1
            while lList[ifwhmminus] > (maxvalue / 2) and ifwhmminus > 0:
                ifwhmminus -= 1
            return imax, maxvalue, ifwhmminus, ifwhmplus

        self.fnLoadStatData()

        # Try to import any fractional protein profiles stored in envelopefrac1/2.dat
        # after such an analysis has been done separately
        try:
            pdf_frac1 = pandas.read_csv(f"{self.spath}/{self.mcmcpath}/envelopefrac1.dat", sep='\s+')
            pdf_frac2 = pandas.read_csv(f"{self.spath}/{self.mcmcpath}/envelopefrac2.dat", sep='\s+')

            for mcmc_iter in range(len(self.diStatResults['Molgroups'])):
                striter = f'iter{mcmc_iter}'
                self.diStatResults['Molgroups'][mcmc_iter]['frac1'] = {}
                self.diStatResults['Molgroups'][mcmc_iter]['frac2'] = {}
                self.diStatResults['Molgroups'][mcmc_iter]['frac1']['zaxis'] = pdf_frac1['zaxis'].tolist()
                self.diStatResults['Molgroups'][mcmc_iter]['frac1']['areaaxis'] = pdf_frac1[striter].tolist()
                self.diStatResults['Molgroups'][mcmc_iter]['frac2']['zaxis'] = pdf_frac2['zaxis'].tolist()
                self.diStatResults['Molgroups'][mcmc_iter]['frac2']['areaaxis'] = pdf_frac2[striter].tolist()
        except IOError:
            print('Did not find any fractional envelopes ...')

        diResults = defaultdict(list)
        l_molgroups = list(self.diStatResults['Molgroups'][0].keys())

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
                    f_int = 0  # total integral in volume per surface area (unit: ??)
                    # taking into account grid spacing of z-axis
                f_int = f_int / mgdict['bilayer.normarea']['areaaxis'][0]

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
            if {'frac1','frac2'}.issubset(l_molgroups):
                diResults['ratio_f1f2'].append(diResults['frac1_INT'][-1] / diResults['frac2_INT'][-1])
                diResults['ratio_f2f1'].append(diResults['frac2_INT'][-1] / diResults['frac1_INT'][-1])

            # percentage water in sub-membrane space and other regions for tBLM and other bilayer
            # get vf_bilayer from molecular group methylene2_x
            vf_bilayer = 0
            i = 1
            while f'bilayer.methylene2_{i}' in mgdict:
                vf_bilayer += float(mgdict[f'bilayer.methylene2_{i}']['headerdata']['nf'])
                i +=1
            # prepare arrays for summing up molgroups
            total_components = numpy.zeros_like(mgdict['bilayer.normarea']['areaaxis'])

            if 'bilayer.headgroup1_1' in l_molgroups:
                f_vol_headgroup1 = numpy.zeros_like(total_components)
                j = 1
                while f'bilayer.headgroup1_{j}' in l_molgroups:
                    total_components += numpy.array(mgdict[f'bilayer.headgroup1_{j}']['areaaxis'])
                    f_vol_headgroup1 += float(mgdict[f'bilayer.headgroup1_{j}']['headerdata']['vol']) * \
                        float(mgdict[f'bilayer.headgroup1_{j}']['headerdata']['nf'])
                    j += 1

                if {'bilayer.tether','bilayer.tetherg','bilayer.normarea','bilayer.bME'}.issubset(l_molgroups):
                    f_vol_submembrane = mgdict['bilayer.normarea']['areaaxis'][0] * \
                        (float(mgdict['bilayer.tether']['headerdata']['l']) +
                         float(mgdict['bilayer.tetherg']['headerdata']['l'])) * vf_bilayer

                    # sub membrane components
                    f_vol_components = float(mgdict['bilayer.bME']['headerdata']['vol']) * \
                        float(mgdict['bilayer.bME']['headerdata']['nf']) + \
                        float(mgdict['bilayer.tether']['headerdata']['vol']) * \
                        float(mgdict['bilayer.tether']['headerdata']['nf']) + \
                        float(mgdict['bilayer.tetherg']['headerdata']['vol']) * \
                        float(mgdict['bilayer.tetherg']['headerdata']['nf']) + f_vol_headgroup1

                    f_total_tether_length = float(mgdict['bilayer.tether']['headerdata']['l']) + \
                        float(mgdict['bilayer.tetherg']['headerdata']['l'])
                    diResults['fTotalTetherLength'].append(f_total_tether_length)

                    f_tether_density = float(mgdict['bilayer.tether']['headerdata']['nf']) / \
                        mgdict['bilayer.normarea']['areaaxis'][0]
                    diResults['fTetherDensity'].append(f_tether_density)

                    total_components += numpy.array(mgdict['bilayer.bME']['areaaxis']) + \
                        numpy.array(mgdict['bilayer.tether']['areaaxis']) + \
                        numpy.array(mgdict['bilayer.tetherg']['areaaxis'])

            if {'bilayer.methylene1_1', 'bilayer.methyl1_1'}.issubset(l_molgroups):
                f_total_lipid1_length = float(mgdict['bilayer.methylene1_1']['headerdata']['l']) + float(
                    mgdict['bilayer.methyl1_1']['headerdata']['l'])
                diResults['fTotalLipid1Length'].append(f_total_lipid1_length)
                i = 1
                while f'bilayer.methylene1_{i}' in mgdict:
                    total_components += numpy.array(mgdict[f'bilayer.methylene1_{i}']['areaaxis']) + \
                                        numpy.array(mgdict[f'bilayer.methyl1_{i}']['areaaxis'])
                    i += 1

            if {'bilayer.methylene2_1','bilayer.methyl2_1'}.issubset(l_molgroups):
                f_total_lipid2_length = float(mgdict['bilayer.methylene2_1']['headerdata']['l']) + float(
                    mgdict['bilayer.methyl2_1']['headerdata']['l'])
                diResults['fTotalLipid2Length'].append(f_total_lipid2_length)

                f_area_per_lipid2 = 0
                i = 1
                while f'bilayer.methylene2_{i}' in mgdict:
                    f_area_per_lipid2 += float(mgdict[f'bilayer.methylene2_{i}']['headerdata']['vol']) * \
                                            float(mgdict[f'bilayer.methylene2_{i}']['headerdata']['nf']) / \
                                            float(mgdict[f'bilayer.methylene2_{i}']['headerdata']['l'])
                    i += 1
                    
                diResults['fAreaPerLipid2'].append(f_area_per_lipid2)

                i = 1
                while f'bilayer.methylene2_{i}' in mgdict:
                    total_components += numpy.array(mgdict[f'bilayer.methylene2_{i}']['areaaxis']) + \
                                        numpy.array(mgdict[f'bilayer.methyl2_{i}']['areaaxis'])
                    i += 1

            if 'bilayer.headgroup2_1' in l_molgroups:
                f_vol_headgroup2 = numpy.zeros_like(total_components)
                j = 1
                while f'bilayer.headgroup2_{j}' in l_molgroups:
                    total_components += numpy.array(mgdict[f'bilayer.headgroup2_{j}']['areaaxis'])
                    f_vol_headgroup2 += float(mgdict[f'bilayer.headgroup2_{j}']['headerdata']['vol']) * \
                        float(mgdict[f'bilayer.headgroup2_{j}']['headerdata']['nf'])
                    j += 1

            if 'bilayer.defect_hc' in l_molgroups:
                total_components = total_components + mgdict['bilayer.defect_hc']['areaaxis']

            if 'bilayer.defect_hg' in l_molgroups:
                total_components = total_components + mgdict['bilayer.defect_hg']['areaaxis']

            if 'protein' in l_molgroups:
                total_components = total_components + mgdict['protein']['areaaxis']
                for i in range(len(total_components)):
                    total_components[i] = min(total_components[i], vf_bilayer * mgdict['normarea']['areaaxis'][i])

            # calculate water and protein fractions if bilayer is present
            if 'bilayer.headgroup1_1' in l_molgroups:
                f_start_hg1 = float(mgdict['bilayer.headgroup1_1']['headerdata']['z']) - 0.5 * float(
                    mgdict['bilayer.headgroup1_1']['headerdata']['l'])
                f_start_hc = float(mgdict['bilayer.methylene1_1']['headerdata']['z']) - 0.5 * float(
                    mgdict['bilayer.methylene1_1']['headerdata']['l'])
                f_start_methyl2 = float(mgdict['bilayer.methyl2_1']['headerdata']['z']) - 0.5 * float(
                    mgdict['bilayer.methyl2_1']['headerdata']['l'])
                f_start_hg2 = float(mgdict['bilayer.headgroup2_1']['headerdata']['z']) - 0.5 * float(
                    mgdict['bilayer.headgroup2_1']['headerdata']['l'])
                f_startbulk = float(mgdict['bilayer.headgroup2_1']['headerdata']['z']) + 0.5 * float(
                    mgdict['bilayer.headgroup2_1']['headerdata']['l'])

                f_step_size = mgdict['bilayer.normarea']['zaxis'][1] - mgdict['bilayer.normarea']['zaxis'][0]

                i_start_hg1 = int(floor(f_start_hg1 / f_step_size + 0.5))
                i_start_hc = int(floor(f_start_hc / f_step_size + 0.5))
                i_start_methyl2 = int(floor(f_start_methyl2 / f_step_size + 0.5))
                i_start_hg2 = int(floor(f_start_hg2 / f_step_size + 0.5))
                i_start_bulk = int(floor(f_startbulk / f_step_size + 0.5))

                ref = mgdict['bilayer.normarea']['areaaxis']
                ratio = 1 - sum(total_components[0:i_start_hg1]) / sum(ref[0:i_start_hg1])
                diResults['WaterFracSubMembrane'].append(ratio)
                ratio = 1 - sum(total_components[i_start_hg1:i_start_hc]) / sum(ref[i_start_hg1:i_start_hc])
                diResults['WaterFracHeadgroup1'].append(ratio)
                ratio = 1 - sum(total_components[i_start_hc:i_start_methyl2]) / sum(ref[i_start_hc:i_start_methyl2])
                diResults['WaterFracLipid1'].append(ratio)
                ratio = 1 - sum(total_components[i_start_methyl2:i_start_hg2]) / sum(ref[i_start_methyl2:i_start_hg2])
                diResults['WaterFracLipid2'].append(ratio)
                ratio = 1 - sum(total_components[i_start_hc:i_start_hg2]) / sum(ref[i_start_hc:i_start_hg2])
                diResults['WaterFracHydrocarbon'].append(ratio)
                ratio = 1 - sum(total_components[i_start_hg2:i_start_bulk]) / sum(ref[i_start_hg2:i_start_bulk])
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
                            ref = mgdict[group]['areaaxis']
                            f_frac_submembrane = sum(ref[0:i_start_hg1]) / sum(ref)
                            f_frac_inner_headgroup = sum(ref[i_start_hg1:i_start_hc]) / sum(ref)
                            f_frac_inner_hydrocarbon = sum(ref[i_start_hc:i_start_methyl2]) / sum(ref)
                            f_frac_outer_hydrocarbon = sum(ref[i_start_methyl2:i_start_hg2]) / sum(ref)
                            f_frac_hydrocarbon = sum(ref[i_start_hc:i_start_hg2]) / sum(ref)
                            f_frac_outer_headgroup = sum(ref[i_start_hg2:i_start_bulk]) / sum(ref)
                            f_frac_headgroups = f_frac_inner_headgroup + f_frac_outer_headgroup
                            f_frac_bulk = sum(ref[i_start_bulk:]) / sum(ref)
                            f_frac_inner_leaflet = f_frac_inner_headgroup + f_frac_inner_hydrocarbon
                            f_frac_outer_leaflet = f_frac_outer_headgroup + f_frac_outer_hydrocarbon

                        diResults[f'FracSubmembrane_{group}'].append(f_frac_submembrane)
                        diResults[f'FracHydrocarbon_{group}'].append(f_frac_hydrocarbon)
                        diResults[f'FracInnerHydrocarbon_{group}'].append(f_frac_inner_hydrocarbon)
                        diResults[f'FracOuterHydrocarbon_{group}'].append(f_frac_outer_hydrocarbon)
                        diResults[f'FracInnerHeadgroup_{group}'].append(f_frac_inner_headgroup)
                        diResults[f'FracOuterHeadgroup_{group}'].append(f_frac_outer_headgroup)
                        diResults[f'FracHeadgroups_{group}'].append(f_frac_headgroups)
                        diResults[f'FracBulk_{group}'].append(f_frac_bulk)
                        diResults[f'FracInnerLeaflet_{group}'].append(f_frac_inner_leaflet)
                        diResults[f'FracOuterLeaflet_{group}'].append(f_frac_outer_leaflet)

                        # calculate peak position and FWHM for spline profile
                        imax, __, ifwhmminus, ifwhmplus = fnFindMaxFWHM(mgdict[group]['areaaxis'])
                        diResults[f'PeakPosition_{group}'].append(mgdict[group]['zaxis'][imax] - f_startbulk)
                        diResults[f'PeakValue_{group}'].append(mgdict[group]['areaaxis'][imax])
                        diResults[f'FWHMMinusPosition_{group}'].append(mgdict[group]['zaxis'][ifwhmminus] - f_startbulk)
                        diResults[f'FWHMPlusPosition_{group}'].append(mgdict[group]['zaxis'][ifwhmplus] - f_startbulk)
                        diResults[f'FWHM_{group}'].append(
                            mgdict[group]['zaxis'][ifwhmplus] - mgdict[group]['zaxis'][ifwhmminus])

        fConfidence = min(1, fConfidence)
        if fConfidence < 0:
            fConfidence = special.erf(-1 * fConfidence / sqrt(2))

        fLowerPercentileMark = 100.0 * (1 - fConfidence) / 2
        fHigherPercentileMark = (100 - fLowerPercentileMark)
        
        results = pandas.DataFrame()
        with open(self.mcmcpath + "/CalculationResults.dat", "w") as file:
            for element, __ in sorted(diResults.items()):
                fLowPerc = stats.scoreatpercentile(diResults[element], fLowerPercentileMark)  # Calculate Percentiles
                fMedian = stats.scoreatpercentile(diResults[element], 50.)
                fHighPerc = stats.scoreatpercentile(diResults[element], fHigherPercentileMark)
                interval = {'element': element, 'lower_conf': fLowPerc, 'median': fMedian, 'upper_conf': fHighPerc}
                results = results.append(interval, ignore_index=True)
                
                sPrintString = '%(el)s  [%(lp)10.4g, %(m)10.4g, %(hp)10.4g] (-%(ld)10.4g, +%(hd)10.4g)'

                soutput = sPrintString % {'el': element, 'lp': fLowPerc, 'ld': (fMedian - fLowPerc), 'm': fMedian,
                                        'hd': (fHighPerc - fMedian), 'hp': fHighPerc}
                file.write(soutput + '\n')
                if verbose: print(soutput)
        return results

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
            diff = max(abs(diffX), abs(diffY)) # which direction more steps?

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
        # result for contour plot of profiles
        for i in range(iNumberOfModels):  # list of lists for each profile
            liContourArray = liContourArray + [[]]
            liContourArrayDimensions = liContourArrayDimensions + [[0., 0., dRhoGrid, 0., 0., dZGrid]]
            # dimension format ymin, ymax, ystep, xmin, xmax, xstep

        for iteration in self.diStatResults['nSLDProfiles']:  # cycle through all individual stat. results
            for i, model in enumerate(iteration):  # cycle through all models
                for l in range(len(model[0])):  # extract nSLD profile data point by point
                    dz = round(model[0][l] / dZGrid) * dZGrid  # round to grid precision
                    drho = round(model[1][l] / dRhoGrid) * dRhoGrid  # round nSLD to grid precision
                    if l != 0:
                        fnMap2Array(liContourArray[i],
                                    liContourArrayDimensions[i],
                                    drhoold, dzold, drho, dz)

                    dzold = dz
                    drhoold = drho

        print('Processing data for output ...')
        for i in range(iNumberOfModels):  # loop through all models
            print('Model %i: %i x %i' % (i, len(liContourArray[i][0]),
                                         len(liContourArray[i])))
            sFileName = f'Cont_nSLD_Array{i}.dat'  # write out array
            with open(sFileName, "w") as file:
                for line in liContourArray[i]:
                    sLine = ' '.join(map(str, line))
                    file.write(sLine + ' \n')

            dRhoMin = liContourArrayDimensions[i][0]
            dRhoMax = liContourArrayDimensions[i][1]
            dRhoStep = liContourArrayDimensions[i][2]
            dZMin = liContourArrayDimensions[i][3]
            dZMax = liContourArrayDimensions[i][4]
            dZStep = liContourArrayDimensions[i][5]

            dZ = dZMin  # write out x-dimension wave
            sFileName = f'Cont_nSLD_DimZ{i}.dat'  # dimension wave has one point extra for Igor
            with open(sFileName, "w") as file:
                while dZ <= dZMax + dZStep:
                    sLine = str(round(dZ / dZGrid) * dZGrid) + '\n'
                    file.write(sLine)
                    dZ = dZ + dZStep

            dRho = dRhoMin  # write out y-dimension wave
            sFileName = f'Cont_nSLD_DimRho{i}.dat'  # dimension wave has one point extra for Igor
            with open(sFileName, "w") as file:
                while dRho <= dRhoMax + dRhoStep:
                    sLine = f'{round(dRho / dRhoGrid) * dRhoGrid}\n'
                    file.write(sLine)
                    dRho += dRhoStep

        if 'Molgroups' in self.diStatResults:
            liContourArray = []  # initiate array
            liContourArrayDimensions = []  # array dimensions
            for _ in self.diStatResults['Molgroups'][0]:  # iterate through molecular group names
                liContourArray = liContourArray + [[]]
                liContourArrayDimensions = liContourArrayDimensions + [[0., 0., dAreaGrid, 0., 0., dZGrid]]
                # dimension format ymin, ymax, ystep, xmin, xmax, xstep

            for iteration in self.diStatResults['Molgroups']:  # cycle through all individual stat. results
                dareaold = 0
                dzold = 0
                for i, molgroup in enumerate(iteration):  # cycle through all molecular groups
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

            for i, molgroup in enumerate(self.diStatResults['Molgroups'][0]):  # loop through all models

                print('%s %i: %i x %i' % (molgroup, i, len(liContourArray[i][0]),
                                          len(liContourArray[i])))
                sFileName = 'Cont_' + molgroup + '_Array' + '.dat'  # write out array
                with open(sFileName, "w") as file:
                    for line in liContourArray[i]:
                        sLine = ' '.join(map(str, line)) + ' \n'
                        file.write(sLine)

                dAreaMin = liContourArrayDimensions[i][0]
                dAreaMax = liContourArrayDimensions[i][1]
                dAreaStep = liContourArrayDimensions[i][2]
                dZMin = liContourArrayDimensions[i][3]
                dZMax = liContourArrayDimensions[i][4]
                dZStep = liContourArrayDimensions[i][5]

                dZ = dZMin  # write out x-dimension wave
                sFileName = 'Cont_' + molgroup + '_DimZ' + '.dat'  # dimension wave has one point extra for Igor
                with open(sFileName, "w") as file:
                    while dZ <= dZMax + dZStep:
                        sLine = str(round(dZ / dZGrid) * dZGrid) + '\n'
                        file.write(sLine)
                        dZ = dZ + dZStep

                dArea = dAreaMin  # write out y-dimension wave
                sFileName = 'Cont_' + molgroup + '_DimArea' + '.dat'  # dimension wave has one point extra for Igor
                with open(sFileName, "w") as file:
                    while dArea <= dAreaMax + dAreaStep:
                        sLine = str(round(dArea / dAreaGrid) * dAreaGrid) + '\n'
                        file.write(sLine)
                        dArea = dArea + dAreaStep

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

            if startindex > stopindex:
                startindex, stopindex = stopindex, startindex

            fsum = 0.5 * (array[startindex] + array[stopindex])
            for i in range(startindex + 1, stopindex):
                fsum += array[i]

            return fsum

        # find maximum values and indizees of half-height points assuming unique solution and steady functions
        def fnMaximumHalfPoint(data):
            fmax = numpy.amax(data)
            point1 = False
            point2 = False
            hm1 = 0
            hm2 = 0

            for i in range(len(data)):
                if data[i] > (fmax / 2) and not point1:
                    point1 = True
                    hm1 = i
                if data[i] < (fmax / 2) and point1 and not point2:
                    point2 = True
                    hm2 = i - 1
                    break

            return fmax, hm1, hm2

        def fnStat(diarea, name, diStat):
            for i in range(len(diarea[list(diarea.keys())[0]])):
                liOnePosition = [iteration[i] for key, iteration in diarea.items() if key != 'zaxis']
                stat = self.fnCalcConfidenceLimits(liOnePosition, method=1)
                diStat[name + '_msigma'].append(stat[1])
                diStat[name].append(stat[2])
                diStat[name + '_psigma'].append(stat[3])

        # initialize Statistical Dictionary
        print('Initializing ...')
        lGroupList = ['substrate', 'siox', 'tether', 'innerhg', 'inner_cg', 'inner_phosphate', 'inner_choline',
                      'innerhc', 'innerch2', 'innerch3', 'outerhc', 'outerch2', 'outerch3', 'outerhg',
                      'outer_cg', 'outer_phosphate', 'outer_choline', 'protein', 'sum', 'water']
        diStat = {}
        for element in lGroupList:
            diStat[element] = []
            diStat[f'{element}_corr'] = []
            diStat[f'{element}_cvo'] = []
            diStat[f'{element}_corr_cvo'] = []

        keylist = list(diStat)
        for element in keylist:
            diStat[f'{element}_msigma'] = []
            diStat[f'{element}_psigma'] = []

        diIterations = {}
        for element in lGroupList:
            diIterations[element] = {}
            diIterations[f'{element}_corr'] = {}
            diIterations[f'{element}_cvo'] = {}
            diIterations[f'{element}_corr_cvo'] = {}

        # pull all relevant molgroups
        # headgroups are allready corrected for protein penetration (_corr) and will be copied over to the _corr
        # entries next
        print('Pulling all molgroups ...')
        print('  substrate ...')
        diIterations['substrate'], __, __ = self.fnPullMolgroupLoader(['bilayer.substrate'])
        print('  siox ...')
        diIterations['siox'], __, __ = self.fnPullMolgroupLoader(['bilayer.siox'])
        print('  tether ...')
        diIterations['tether'], __, __ = self.fnPullMolgroupLoader(['bilayer.bME', 'bilayer.tetherg', 'bilayer.tether',
                                                                    'bilayer.tether_bme', 'bilayer.tether_free',
                                                                    'bilayer.tether_hg'])

        print('  innerhg ...')
        grouplist = []
        i = 1
        while f'bilayer.headgroup1_{i}' in self.diStatResults['Molgroups'][0]:
            grouplist.append(f'bilayer.headgroup1_{i}')
            i += 1

        diIterations['innerhg'], __, __ = self.fnPullMolgroupLoader(grouplist)
        diIterations['inner_cg'], __, __ = self.fnPullMolgroupLoader(['bilayer.headgroup1_1.carbonyl_glycerol'])
        diIterations['inner_phosphate'], __, __ = self.fnPullMolgroupLoader(['bilayer.headgroup1_1.phosphate'])
        diIterations['inner_choline'], __, __ = self.fnPullMolgroupLoader(['bilayer.headgroup1_1.choline'])

        print('  innerhc ...')
        gl_methylene = []
        gl_methyl = []
        
        i = 1
        while f'bilayer.methylene1_{i}' in self.diStatResults['Molgroups'][0]:
            gl_methylene.append(f'bilayer.methylene1_{i}')
            i += 1
        gl_methylene.append('bilayer.tether_methylene')
        
        i = 1
        while f'bilayer.methyl1_{i}' in self.diStatResults['Molgroups'][0]:
            gl_methyl.append(f'bilayer.methyl1_{i}')
            i += 1
        gl_methyl.append('bilayer.tether_methyl')
        
        diIterations['innerhc'], __, __ = self.fnPullMolgroupLoader(gl_methylene+gl_methyl)
        diIterations['innerch2'], __, __ = self.fnPullMolgroupLoader(gl_methylene)
        diIterations['innerch3'], __, __ = self.fnPullMolgroupLoader(gl_methyl)

        print('  outerhc ...')
        gl_methylene = []
        i = 1
        while f'bilayer.methylene2_{i}' in self.diStatResults['Molgroups'][0]:
            gl_methylene.append(f'bilayer.methylene2_{i}')
            i += 1

        gl_methyl = []
        i = 1
        while f'bilayer.methyl2_{i}' in self.diStatResults['Molgroups'][0]:
            gl_methyl.append(f'bilayer.methyl2_{i}')
            i += 1

        diIterations['outerhc'], __, __ = self.fnPullMolgroupLoader(gl_methylene+gl_methyl)
        diIterations['outerch2'], __, __ = self.fnPullMolgroupLoader(gl_methylene)
        diIterations['outerch3'], __, __ = self.fnPullMolgroupLoader(gl_methyl)

        print('  outerhg ...')
        grouplist = []
        i = 1
        while f'bilayer.headgroup2_{i}' in self.diStatResults['Molgroups'][0]:
            grouplist.append(f'bilayer.headgroup2_{i}')
            i += 1
                
        diIterations['outerhg'], __, __ = self.fnPullMolgroupLoader(grouplist)
        diIterations['outer_cg'], __, __ = self.fnPullMolgroupLoader(['bilayer.headgroup2_1.carbonyl_glycerol'])
        diIterations['outer_phosphate'], __, __ = self.fnPullMolgroupLoader(['bilayer.headgroup2_1.phosphate'])
        diIterations['outer_choline'], __, __ = self.fnPullMolgroupLoader(['bilayer.headgroup2_1.choline'])

        print('  protein ...')
        diIterations['protein'], __, __ = self.fnPullMolgroupLoader(['protein'])

        # save z-axis
        diStat['zaxis'] = numpy.copy(diIterations['substrate']['zaxis'])

        # shallow copies of the uncorrected data into the corrected dictionaries
        # and the values will be replaced by their modifications step by step
        for element in lGroupList:
            diIterations[f'{element}_corr'] = diIterations[element].copy()

        # loop over all iterations and apply the corrections / calculations
        print('Applying corrections ...\n')
        for i in range(len(list(diIterations['substrate'].keys())) - 1):

            key = f'iter{i}'
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
            # if no substrate, use maximum bilayer area as area per lipid
            if substrate.min() == substrate.max() == 0:
                areaperlipid = maxbilayerarea

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
                excess = min(sum[i] + protein[i] - maxbilayerarea, (innerhc[i] + outerhc[i]))
                if excess > 0:
                    innerhc_corr[i] -= excess * innerhc[i] / (innerhc[i] + outerhc[i])
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
                diIterations[f'{element}_cvo'][key] = diIterations[element][key] / areaperlipid
                diIterations[f'{element}_corr_cvo'][key] = diIterations[f'{element}_corr'][key] / areaperlipid

        # calculate the statisics
        print('Calculating statistics ...\n')
        for element in lGroupList:
            if element != 'zaxis':
                fnStat(diIterations[element], element, diStat)
                fnStat(diIterations[f'{element}_corr'], f'{element}_corr', diStat)
                fnStat(diIterations[f'{element}_cvo'], f'{element}_cvo', diStat)
                fnStat(diIterations[f'{element}_corr_cvo'], f'{element}_corr_cvo', diStat)

        print('Saving data to bilayerplotdata.dat ...\n')
        self.Interactor.fnSaveSingleColumns(f'{self.mcmcpath}/bilayerplotdata.dat', diStat)

    def fnGetNumberOfModelsFromSetupC(self):
        with open(self.setupfilename, "r") as file:  # open setup.c
            data = file.readlines()
        smatch = compile(r'define\s+MODELS\s+(.+?)\n', IGNORECASE | VERBOSE)
        for line in data:  # search through setup.c
            if smatch.search(line):  # searching for MODELS constant
                i = smatch.search(line).group(1)
                return int(i)
        return 0

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
        if self.diStatResults == {}:
            try:
                self.diStatResults = self.fnLoadObject(f'{self.mcmcpath}StatDataPython.dat')
                print('Loaded statistical data from StatDataPython.dat')
            except IOError:
                print('No StatDataPython.dat.')
                print('Recreate statistical data from sErr.dat.')
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
                self.fnRecreateStatistical()

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
                    f += fGrid
                elif f > xdata[-1]:  # fill remaining cells with last value
                    liInterpolated[0].append(f)
                    liInterpolated[1].append(ydata[-1])
                    f += fGrid
                else:  # at least one data point surpassed by f
                    while (f > xdata[ix]) and (ix < (len(xdata) - 1)):  # searching first data point past f
                        ix += 1
                    if f < xdata[ix]:  # was there a data point past f?
                        LowerX = ix - 1  # calculate data surrounding f
                        UpperX = ix
                        fDataGrid = xdata[UpperX] - xdata[LowerX]
                        fInterpolate = ydata[LowerX] * ((f - xdata[LowerX]) / fDataGrid  # do the interpolation
                                                        ) + ydata[UpperX] * ((xdata[UpperX] - f) / fDataGrid)
                        liInterpolated[0].append(f)
                        liInterpolated[1].append(fInterpolate)
                        f += fGrid
                    elif f == xdata[ix]:  # no interpolation needed
                        liInterpolated[0].append(xdata[ix])
                        liInterpolated[1].append(ydata[ix])
                        f += fGrid

            return liInterpolated

        def fnStoreEnvelope(liStoredEnvelopes, liStoredEnvelopeHeaders, envelope, percentile):
            iposition = len(liStoredEnvelopes) / 2
            liStoredEnvelopes.insert(iposition, envelope[1])
            liStoredEnvelopes.insert(iposition, envelope[0])
            liStoredEnvelopeHeaders.insert(iposition, str(1 - percentile / 2))
            liStoredEnvelopeHeaders.insert(iposition, str(percentile / 2))

        def fnSaveEnvelopes(liStoredEnvelopes, liStoredEnvelopeHeaders, iModel, fMin, fMax, fGrid, fSigma):
            file = open(f'Envelopes{iModel}.dat', "w")

            if fSigma == 0:  # save all computed envelopes
                liSaveEnvelopeHeaders = liStoredEnvelopeHeaders
                liSaveEnvelopes = liStoredEnvelopes

            else:  # save only multiples of fSigma in Sigma
                liSaveEnvelopeHeaders = []
                liSaveEnvelopes = []

                fmult = 0.
                while True:
                    fConfidence = special.erf(fmult / sqrt(2))  # upper and lower Percentiles for fmult*sigma
                    fLowerPerc = (1 - fConfidence) / 2
                    fUpperPerc = 1 - (1 - fConfidence) / 2

                    fStepSize = 1 / float(len(liStoredEnvelopeHeaders))  # positions in envelopes list
                    iLowerPerc = int(floor(fLowerPerc / fStepSize))
                    iUpperPerc = int(ceil(fUpperPerc / fStepSize))

                    if (iLowerPerc == 0) or (iUpperPerc >= len(liStoredEnvelopeHeaders) - 1):
                        break

                    liSaveEnvelopeHeaders.insert(0, f'minus{fmult}sigma')
                    liSaveEnvelopes.insert(0, liStoredEnvelopes[iLowerPerc])

                    if iUpperPerc != iLowerPerc:
                        liSaveEnvelopeHeaders.append(f'plus{fmult}sigma')
                        liSaveEnvelopes.append(liStoredEnvelopes[iUpperPerc])

                    fmult += fSigma

                liSaveEnvelopeHeaders.insert(0, 'LowerEnvelope')
                liSaveEnvelopes.insert(0, liStoredEnvelopes[0])
                liSaveEnvelopeHeaders.append('UpperEnvelope')
                liSaveEnvelopes.append(liStoredEnvelopes[len(liStoredEnvelopeHeaders) - 1])

            file.write("z ")
            for element in liSaveEnvelopeHeaders:
                if fSigma == 0:
                    file.write(f"p{element} ")
                else:
                    file.write(f"{element} ")
            file.write("\n")

            f = fMin
            for i in range(len(liSaveEnvelopes[0])):
                file.write(f'{f} ')
                for element in liSaveEnvelopes:
                    file.write(f'{element[i]} ')
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
            print(f'Analyzing model {iModel} ...')

            liStoredEnvelopes = []  # List of stored envelopes
            liStoredEnvelopeHeaders = []  # and their headers

            fMin = fMax = 0  # initializing boundaries
            profilelist = []  # extracting all profiles related to the actual
            for iteration in self.diStatResults['nSLDProfiles']:  # model
                profilelist.append(iteration[iModel][:])
                fMin = min(fMin, profilelist[-1][0][0])
                fMax = max(fMax, profilelist[-1][0][-1])
            fMax = floor((fMax - fMin) / fGrid) * fGrid + fMin  # make fMax compatible with fGrid and fMin

            print('Rebinning data...')
            iNumberOfProfiles = len(profilelist)
            for i in range(iNumberOfProfiles):
                profilelist[i] = fnInterpolateData(profilelist[i][0], profilelist[i][1], fMin, fMax, fGrid)

            if not shortflag:
                for iPercentile in range(iNumberOfProfiles):

                    print(f'Calculating {1 - float(iPercentile) / float(iNumberOfProfiles)} percentile...')

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
                    print(f'Short: Calculating {1 - float(iPercentile) / float(iNumberOfProfiles)} percentile...')

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
            # change setup.c to quasi fix all parameters
            self.fnWriteConstraint2SetupC(liAddition)  # write out
            liAddition = [f"{self.diParameters[parameter]['variable']} = {self.diParameters[parameter]['value']};\n" for parameter in liParameters] 
            call(["rm", "-f", "mol.dat"])
            self.fnMake()  # compile changed setup.c
            call(["./fit", "-o"])  # write out profile.dat and fit.dat
            call(["sync"])  # synchronize file system
            sleep(1)  # wait for system to clean up

        finally:
            self.fnRemoveBackup()

        self.fnLoadMolgroups()  # Load Molgroups into self.diMolgroups
        normarea = 100  # default, if no normarea is provided

        with open('areatab.dat', 'w') as File: # save all molgroupdata in table for loading
        # into Igor
            File.write('z ')  # write header line
            for element in self.diMolgroups:
                File.write(self.diMolgroups[element]['headerdata']['ID'] + ' ')
            File.write('summol water waterperc\n')

            element = list(self.diMolgroups.keys())[0]
            datalength = len(self.diMolgroups[element]['zaxis'])

            for i in range(datalength):
                File.write(f"{self.diMolgroups[element]['zaxis'][i]} ")
                sum = 0
                normarea = 0
                for element in self.diMolgroups:
                    if element != 'normarea':
                        sum += self.diMolgroups[element]['areaaxis'][i]
                    else:
                        normarea = self.diMolgroups[element]['areaaxis'][i]
                    File.write(f"{self.diMolgroups[element]['areaaxis'][i]} ")
                File.write(f'{sum} {normarea-sum} {(normarea-sum) / normarea}\n')

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

            if element == 'normarea':  # Get normarea for nSLD calculations
                normarea = self.diMolgroups[element]['areaaxis'][0]
            else:  # do not sum up normarea indicator
                if not areasum:
                    for i in range(length):
                        areasum.append(0.)
                        nSLsum.append(0.)
                for i in range(length):
                    areasum[i] = areasum[i] + area[i]
                    nSLsum[i] = nSLsum[i] + nSL[i]

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
        
        fBulknSLDs = [6.34, 4.00, 0.00, -0.56]
        for fBulknSLD in fBulknSLDs:  # loop over contrast mixtures            
            for i in range(length):  # calculate nSLD for several cases
                nSLDSum[i] = nSLsum[i] * 1E2 / (stepsize * normarea) + fBulknSLD * (1 - (areasum[i] / normarea))
            plt.subplot(224)
            plt.plot(zax, nSLDSum, label=f'nSLDsum CM{fBulknSLD}')

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

        plt.suptitle(f'{plotname} \n\n Area and nSLD Profile')
        plt.savefig(f'{plotname}', format='png')
        plt.show()

    # -------------------------------------------------------------------------------

    def fnPlotFit(self, plotname):

        from matplotlib.font_manager import fontManager, FontProperties

        font = FontProperties(size='x-small')

        plt.figure(1, figsize=(14, 10))

        iCounter = 0
        sfilename = 'fit' + str(iCounter) + '.dat'
        while path.isfile(sfilename):
            with open(sfilename, 'r') as file:
                data = file.readlines()
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

            iCounter += 1
            sfilename = 'fit' + str(iCounter) + '.dat'

        iCounter = 0
        sfilename = 'profile' + str(iCounter) + '.dat'
        while path.isfile(sfilename):
            with open(sfilename, 'r') as file:
                data = file.readlines()
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
            plt.plot(zlist, rholist, label=f'profile{iCounter}')

            iCounter = iCounter + 1
            sfilename = f'profile{iCounter}.dat'

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

        plt.suptitle(f'{plotname} \n\n Reflectivity data and fit, residuals and profile')
        plt.savefig(f'{plotname}_fit.png', format='png')
        plt.show()

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

    def fnPullMolgroup(self, liMolgroupNames, sparse=0, verbose=True):
        """
        Calls Function that recreates statistical data and extracts only area and nSL profiles for
        submolecular groups whose names are given in liMolgroupNames. Those groups
        are added for each iteration and a file pulledmolgroups.dat is created.
        A statistical analysis area profile containing the median, sigma, and
        2 sigma intervals are put out in pulledmolgroupsstat.dat.
        Saves results to file
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

        self.Interactor.fnSaveSingleColumns(self.mcmcpath+'/pulledmolgroups_area.dat', diarea)
        self.Interactor.fnSaveSingleColumns(self.mcmcpath+'/pulledmolgroups_nsl.dat', dinsl)
        self.Interactor.fnSaveSingleColumns(self.mcmcpath+'/pulledmolgroups_nsld.dat', dinsld)
        self.Interactor.fnSaveSingleColumns(self.mcmcpath+'/pulledmolgroupsstat.dat', diStat)

    def fnPullMolgroupLoader(self, liMolgroupNames, sparse=0, verbose=True):
        """
        Function recreates statistical data and extracts only area and nSL profiles for
        submolecular groups whose names are given in liMolgroupNames. Those groups
        are added for each iteration and a file pulledmolgroups.dat is created.
        A statistical analysis area profile containing the median, sigma, and
        2 sigma intervals are put out in pulledmolgroupsstat.dat.
        """

        if self.diStatResults != {}:
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

    def fnRecreateStatistical(self, bRecreateMolgroups=True, verbose=False):
        """
        Recreates profile and fit data associated with parameter stats
        """
        if self.Interactor.problem is not None:
            problem = self.Interactor.problem
        else:
            problem = self.Interactor.fnRestoreFitProblem()

        j = 0
        self.diStatResults['nSLDProfiles'] = []
        self.diStatResults['Molgroups'] = []

        for iteration in range(self.diStatResults['NumberOfStatValues']):
            try:
                # appends a new list for profiles for the current MC iteration
                self.diStatResults['nSLDProfiles'].append([])
                liParameters = list(self.diParameters.keys())
                # sort by number of appereance in setup.c
                liParameters = sorted(liParameters, key=lambda keyitem: self.diParameters[keyitem]['number'])
                bConsistency = set(liParameters).issubset(self.diStatResults['Parameters'].keys())
                if bConsistency:
                    if verbose:
                        print(f'Processing parameter set {j}.\n')
                    p = []
                    for parameter in liParameters:
                        val = self.diStatResults['Parameters'][parameter]['Values'][iteration]
                        # Rescaling is currently disabled
                        # if 'rho_' in parameter or 'background' in parameter:
                        #    val *= 1E6
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

                    # explicit saving is not needed anymore
                    # garefl is automatically saved and refl1d / bumps contains a directory in problem.moldat
                    # self.Interactor.fnSaveMolgroups(problem)

                    # distinguish between FitProblem and MultiFitProblem
                    if 'models' in dir(problem):
                        for M in problem.models:
                            z, rho, irho = self.Interactor.fnRestoreSmoothProfile(M)
                            self.diStatResults['nSLDProfiles'][-1].append((z, rho, irho))
                            if verbose:
                                print(M.chisq())
                    else:
                        z, rho, irho = self.Interactor.fnRestoreSmoothProfile(problem)
                        self.diStatResults['nSLDProfiles'][-1].append((z, rho, irho))
                        if verbose:
                            print(problem.chisq())

                    if bRecreateMolgroups:
                        self.diMolgroups = self.Interactor.fnRestoreMolgroups(problem)
                        # append molgroup information to self.diStatResults
                        self.diStatResults['Molgroups'].append(self.diMolgroups)
                else:
                    raise RuntimeError('Statistical error data and setup file do not match')
            finally:
                j += 1

        # save stat data to disk, if flag is set
        if self.save_stat_data:
            self.fnSaveObject(self.diStatResults, f'{self.spath}/{self.mcmcpath}/StatDataPython.dat')

    @staticmethod
    def fnSaveObject(save_object, sFileName):
        import pickle
        
        with open(sFileName, "wb") as file:
            pickle.dump(save_object, file)

    def fnSimulateData(self, qmin=0.008, qmax=0.325, s1min=0.108, s1max=4.397, s2min=0.108, s2max=4.397, tmin=18,
                       tmax=208, nmin=11809, rhomin=-0.56e-6, rhomax=6.34e-6, cbmatmin=1.1e-5, cbmatmax=1.25e-6,
                       mode='water', pre=1, qrange=0):

        # simulates reflectivity based on a parameter file called simpar.dat
        # requires a compiled and ready to go fit whose fit parameters are modified and fixed
        # data files currently need to be named sim#.dat
        # The method is designed to be useable with reflectivity and SANS.

        # Load Parameters and modify setup.cc
        self.diParameters, _ = self.Interactor.fnLoadParameters()

        try:
            liParameters = list(self.diParameters.keys())
            liParameters = sorted(liParameters, key=lambda keyitem: self.diParameters[keyitem]['number'])

            # simpar contains the parameter values to be simulated
            # this could be done fileless
            simpar = pandas.read_csv(self.spath+'/simpar.dat', sep='\s+', header=None, names=['par', 'value'],
                                     skip_blank_lines=True, comment='#')

            diAddition = {}
            for parameter in liParameters:
                diAddition[parameter] = simpar[simpar.par == parameter].iloc[0][1]

            # if q-range is changing, back up original sim dat or reload previously backed up data
            # to always work with the same set of original data
            if qrange != 0:
                self.Interactor.fnBackupSimdat()

            # extend simulated data to q-range as defined
            i = 0
            while path.isfile(f'{self.spath}/sim{i}.dat') and qrange != 0:
                comments = general.extract_comments_from_file(f'{self.spath}/sim{i}.dat', "#")
                simdata = pandas.read_csv(f'{self.spath}/sim{i}.dat', sep='\s+', skip_blank_lines=True,
                                          comment='#')
                # first cut data at qrange
                simdata = simdata[(simdata.Q <= qrange)]
                # now add data points in case the q-range is too short
                while simdata['Q'].iloc[-1] < qrange:
                    newframe = pandas.DataFrame(simdata[-1:], columns=simdata.columns)
                    newframe['Q'].iloc[-1] = 2 * simdata['Q'].iloc[-1] - simdata['Q'].iloc[-2]
                    simdata = simdata.append(newframe)
                simdata.to_csv(f'{self.spath}/sim{i}.dat', sep=' ', index=None)
                general.add_comments_to_start_of_file(f'{self.spath}/sim{i}.dat', comments)
                i += 1

            # simulate data, works on sim.dat files
            self.Interactor.fnSimulateData(diAddition)
            # simulate error bars, works on sim.dat files
            self.Interactor.fnSimulateErrorBars(simpar, qmin=qmin, qmax=qmax, s1min=s1min, s1max=s1max, s2min=s2min,
                                                s2max=s2max, tmin=tmin, tmax=tmax, nmin=nmin, rhomin=rhomin,
                                                rhomax=rhomax, cbmatmin=cbmatmin, cbmatmax=cbmatmax, mode=mode,
                                                pre=pre)
        finally:
            pass

    def fnStatTable(self, sTableName, fConfidence):

        def fnTexFormatf(fLow, fMed, fHigh):

            def fnDetPrec(fF):
                if fabs(fF) < 1e-4:  # applies to nSLDs
                    fF *= 1E6
                fPrec = ceil(log10(fF) * (-1)) if fF > 0 else 0 #determine precision
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

        with open(sTableName, 'r') as file: # load in template
            template = file.readlines()

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
                    sPrec = f'%(#).{iPrec}f'
                    temp = '$'  # Latex format
                    temp += (sPrec % {'#': fMedianDiff})
                    temp += '_{'
                    temp = temp + '-' + (sPrec % {'#': fLowPercDiff})
                    temp += '}^{'
                    temp = temp + '+' + (sPrec % {'#': fHighPercDiff})
                    temp += '}$'
                    splitline[i] = temp
            table.append(' '.join(splitline) + '\n')

        with open("StatTable.tex", 'w') as file: # save table to file
            file.writelines(table)

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

                    if path.isfile(sDirName + '/' + 'StatFinished'):  # look up if process is finished
                        if path.isfile(sMode + 'Err.dat'):
                            with open(sDirName + '/' + sMode + 'Err.dat', 'r') as file:  # get statistical data from rs subdir
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
