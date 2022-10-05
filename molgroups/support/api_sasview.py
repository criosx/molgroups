from __future__ import print_function
import pandas
import os
import pathlib

import sasmodels.data

from molgroups.support import api_bumps


class CSASViewAPI(api_bumps.CBumpsAPI):
    def __init__(self, spath='.', mcmcpath='.', runfile='', load_state=True):
        super().__init__(spath, mcmcpath, runfile, load_state=load_state)

    def fnLoadData(self, filename):
        """
        Load all data files with the basefilenam filename into a list of Pandas dataframes.
        Each list element is itself a list of [comments, simdata]. It will load n files with the name
        basestem{i}.basesuffix, whereby 'i' is an index from 0 to n-1.
        """
        def _load(stem, suffix):
            ds = sasmodels.data.load_data(os.path.join(self.spath, stem + suffix))
            data = pandas.DataFrame({'Q': ds.x, 'I': ds.y, 'dI': ds.dy, 'dQ': ds.dx})
            return [[], data]
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


    def fnSimulateData(self, diNewPars, liData, data_column='I'):
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
                qvec, scatt = M.fitness.theory()
                liData[i][1][data_column] = scatt
                # liData[i][1]['Q'] = qvec
                i += 1
        else:
            self.problem.chisq()
            qvec, scatt = self.problem.fitness.theory()
            liData[0][1][data_column] = scatt
            # liData[0][1]['Q'] = qvec

        return liData

    def fnSimulateErrorBars(self, simpar, liData):
        """
        Placeholder.
        """
        return liData
