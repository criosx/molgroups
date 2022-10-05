from __future__ import print_function
from math import fabs, pow, sqrt
import pathlib
from random import seed, normalvariate, random
from re import VERBOSE, IGNORECASE, compile
import pandas
import shutil
import os

from molgroups.support import general
from molgroups.support import api_bumps


# Refl1D methods will be used if a storage directory for a Markov Chain Monte Carlo (MCMC)
# error analysis are found.
# The MCMC directory is called 'MCMC'
# The refl1d script name has to be run.py.
class CRefl1DAPI(api_bumps.CBumpsAPI):
    def __init__(self, spath='.', mcmcpath='.', runfile='', load_state=True):
        super().__init__(spath, mcmcpath, runfile, load_state=load_state)

    def fnLoadData(self, filename):
        """
        Load all data files with the basefilenam filename into a list of Pandas dataframes.
        Each list element is itself a list of [comments, simdata]. It will load n files with the name
        basestem{i}.basesuffix, whereby 'i' is an index from 0 to n-1.
        """

        def _load(stem, suffix):
            comments = general.extract_comments_from_file(os.path.join(self.spath, stem + suffix), "#")
            data = pandas.read_csv(os.path.join(self.spath, stem + suffix), sep='\s+', skip_blank_lines=True,
                                   comment='#', header=None)
            # if there was a header move it from first row to header
            if any(data.iloc[0].apply(lambda x: isinstance(x, str))):
                data = data[1:].reset_index(drop=True).rename(columns=data.iloc[0])
                data = data.astype(float)
            return [comments, data]

        stem = pathlib.Path(filename).stem
        suffix = pathlib.Path(filename).suffix
        liData = []
        if os.path.isfile(os.path.join(self.spath, stem + suffix)):
            liData.append(_load(stem, suffix))
        else:
            i = 0
            while True:
                if os.path.isfile(os.path.join(self.spath, stem + str(i) + suffix)):
                    liData.append(_load(stem + str(i), suffix))
                    i += 1
                else:
                    break
        return liData

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

    def fnSaveData(self, basefilename, liData):
        """
        Saves all frames and comments in liData to files with the basefilenam filename.
        Each list element is itself a list of [comments, simdata]. It will save n files with the name
        basestem{i}.basesuffix, whereby 'i' is an index from 0 to n-1.
        """
        def _save(stem, suffix, frame, comment):
            frame.to_csv(os.path.join(self.spath, stem + suffix), sep=' ', index=None, header=None)
            general.add_comments_to_start_of_file(os.path.join(self.spath, stem + suffix), comment)

        stem = pathlib.Path(basefilename).stem
        suffix = pathlib.Path(basefilename).suffix
        if len(liData) == 1:
            _save(stem, suffix, liData[0][1], liData[0][0])
        else:
            for i in range(len(liData)):
                _save(stem + str(i), suffix, liData[i][1], liData[i][0])


    def fnSimulateErrorBars(self, simpar, liData, qmin=0.008, qmax=0.325, s1min=0.108, s1max=4.397, s2min=0.108,
                            s2max=4.397, tmin=18, tmax=208, nmin=11809, rhomin=-0.56e-6, rhomax=6.34e-6,
                            cbmatmin=1.1e-5, cbmatmax=1.25e-6, mode='water', pre=1):

        last_defined_rho_solv = 0
        for i in range(len(liData)):
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

            simdata = liData[i][1]
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

        return liData

    def fnSimulateData(self, diNewPars, liData, data_column='R'):
        liParameters = list(self.diParameters.keys())
        # sort by number of appereance in runfile
        liParameters = sorted(liParameters, key=lambda keyitem: self.diParameters[keyitem]['number'])
        for element in liParameters:
            if element not in list(diNewPars.keys()):
                print('Parameter '+element+' not specified.')
                # check failed -> exit method
                return
            else:
                print(element + ' ' + str(diNewPars[element]))

        p = [diNewPars[parameter] for parameter in liParameters]
        self.problem.setp(p)
        self.problem.model_update()

        # TODO: By calling .chisq() I currently force an update of the cost function. There must be a better way
        if 'models' in dir(self.problem):
            i = 0
            for M in self.problem.models:
                M.chisq()
                qvec, scatt = M.fitness.reflectivity()
                liData[i][1][data_column] = scatt
                liData[i][1]['Q'] = qvec
                i += 1
        else:
            self.problem.chisq()
            qvec, scatt = self.problem.fitness.reflectivity()
            liData[0][1][data_column] = scatt
            liData[0][1]['Q'] = qvec

        return liData


