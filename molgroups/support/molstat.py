from __future__ import print_function
from math import floor, sqrt
from collections import defaultdict
from numpy import average
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
                      'Values' key contains ordered list of all MC parameters
                      'LowPerc' contains lower percentile
                      'Median' contains median
                      'HighPerc' contains higher percentile
                      'MaxParameterLength' contains length of the longest parameter name

            (from fnLoadStatData)
            'NumberOfStatValues' number of MC points
            'nSLDProfiles' contains all profiles 'nSLDProfiles'[MCiteration][model][z,rho]
            'Molgroups' contains a list of dictionaries storing all molgroups
                'Molgroups'[MCiteration] {'molgroupname' {'zaxis' [z]
                                                          'areaxis' [area]
                                                          'nslaxis' [nsl]
                                                          'headerdata' {'property' value}
                                                          'property1' value
                                                          'property2' value}}
           'Results' contains a dictionary of derived values from fit paramters for post-analysis
                'Results' {'valuename' {value}}
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

    def fnAnalyzeStatFile(self, fConfidence=-1, sparse=0):  # summarizes stat file

        self.fnLoadParameters()  # Load Parameters for limits
        self.fnLoadStatData(sparse)  # Load data from file into list

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
        # print('Maximum deviation from average over last %(iHL)d iterations: %(maxc).4f' %
        #       {'iHL': iHistoryLength, 'maxc': fMaxConvergence})
        print(f'Confidence level: {fConfidence:.4f}')

    def fnCalculateMolgroupProperty(self, fConfidence, verbose=True):
        def fnFindMaxFWHM(lList):
            maxvalue = max(lList)

            if isinstance(lList, list):
                imax = lList.index(maxvalue)
                length = len(lList)
            else:
                imax = numpy.argmax(lList)
                length = lList.shape[0]

            ifwhmplus = ifwhmminus = imax
            while lList[ifwhmplus] > (maxvalue / 2) and ifwhmplus < (length - 1):
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
                    f_int = 0  # total integral in volume per surface area (unit: Å)
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
            if {'frac1', 'frac2'}.issubset(l_molgroups):
                diResults['ratio_f1f2'].append(diResults['frac1_INT'][-1] / diResults['frac2_INT'][-1])
                diResults['ratio_f2f1'].append(diResults['frac2_INT'][-1] / diResults['frac1_INT'][-1])

            # percentage water in sub-membrane space and other regions for tBLM and other bilayer
            # get vf_bilayer from molecular group methylene2_x
            vf_bilayer = 0
            i = 1
            while f'bilayer.methylene2_{i}' in mgdict:
                vf_bilayer += float(mgdict[f'bilayer.methylene2_{i}']['headerdata']['nf'])
                i += 1
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

                if {'bilayer.tether', 'bilayer.tetherg', 'bilayer.normarea', 'bilayer.bME'}.issubset(l_molgroups):
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

            if {'bilayer.methylene2_1', 'bilayer.methyl2_1'}.issubset(l_molgroups):
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
                    total_components[i] = min(total_components[i], vf_bilayer *
                                              mgdict['bilayer.normarea']['areaaxis'][i])

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
                interval = pandas.DataFrame(data={'element': element, 'lower_conf': fLowPerc, 'median': fMedian,
                                                  'upper_conf': fHighPerc}, index=[0])
                results = pandas.concat([results, interval], axis=0, ignore_index=True)

                sPrintString = '%(el)s  [%(lp)10.4g, %(m)10.4g, %(hp)10.4g] (-%(ld)10.4g, +%(hd)10.4g)'

                soutput = sPrintString % {'el': element, 'lp': fLowPerc, 'ld': (fMedian - fLowPerc), 'm': fMedian,
                                          'hd': (fHighPerc - fMedian), 'hp': fHighPerc}
                file.write(soutput + '\n')
                if verbose:
                    print(soutput)

        return results

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

    def fnCreateBilayerPlotData(self, custom_groups=None):
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
            hm1 = 0
            hm2 = 0

            for i in range(len(data)):
                if data[i] > (fmax / 2) and not point1:
                    point1 = True
                    hm1 = i
                if data[i] < (fmax / 2) and point1:
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
        if custom_groups is None:
            custom_groups = []
        lGroupList = ['substrate', 'siox', 'tether', 'innerhg', 'inner_cg', 'inner_phosphate', 'inner_choline',
                      'innerhc', 'innerch2', 'innerch3', 'outerhc', 'outerch2', 'outerch3', 'outerhg',
                      'outer_cg', 'outer_phosphate', 'outer_choline', 'protein', 'sum', 'water'] + custom_groups

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

        diIterations['innerhc'], __, __ = self.fnPullMolgroupLoader(gl_methylene + gl_methyl)
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

        diIterations['outerhc'], __, __ = self.fnPullMolgroupLoader(gl_methylene + gl_methyl)
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

        for grp in custom_groups:
            print(f'  {grp} ...')
            diIterations[grp], __, __ = self.fnPullMolgroupLoader([grp])

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
        self.diStatResults['Molgroups'] = []
        self.diStatResults['Results'] = []

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

                # Recreate Molgroups and Derived Results Dictionaries
                self.diMolgroups, self.diResults = self.Interactor.fnLoadMolgroups(problem)
                self.diStatResults['Molgroups'].append(self.diMolgroups)

            finally:
                j += 1

        # save stat data to disk, if flag is set
        if self.save_stat_data:
            self.fnSaveObject(self.diStatResults, f'{self.spath}/{self.mcmcpath}/StatDataPython.dat')

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
                       t_total=None, mode='water', lambda_min=0.1, verbose=True):
        """
        Simulates scattering based on a parameter file called simpar.dat
        requires a ready-to-go fit whose fit parameters are modified and fixed
        The basename can refer to a set of data files with integer indizes before the suffix
        """

        # Load Parameters
        self.diParameters, _ = self.Interactor.fnLoadParameters()

        try:
            liParameters = list(self.diParameters.keys())
            liParameters = sorted(liParameters, key=lambda keyitem: self.diParameters[keyitem]['number'])

            # simpar contains the parameter values to be simulated
            # this could be done fileless
            simpar = pandas.read_csv(self.spath + '/simpar.dat', sep='\s+', header=None, names=['par', 'value'],
                                     skip_blank_lines=True, comment='#')
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
            self.Interactor.fnSaveData(basefilename, liData)

        finally:
            pass
