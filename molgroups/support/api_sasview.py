from __future__ import print_function
import pandas
import os
import pathlib

import sasmodels.data

from molgroups.support import general
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
            # https://github.com/sansigormacros/ncnrsansigormacros/wiki/NCNROutput1D_IQ
            ds = sasmodels.data.load_data(os.path.join(self.spath, stem + suffix))
            data = pandas.DataFrame({'Q': ds.x, 'I': ds.y, 'dI': ds.dy, 'dQ': ds.dx})
            data.columns = ['Q', 'I', 'dI', 'dQ']
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

    def fnSaveData(self, basefilename, liData):
        """
        Saves all frames and comments in liData to files with the basefilenam filename.
        Each list element is itself a list of [comments, simdata]. It will save n files with the name
        basestem{i}.basesuffix, whereby 'i' is an index from 0 to n-1.
        """
        def _save(stem, suffix, frame, comment):
            frame.to_csv(os.path.join(self.spath, stem + suffix), sep=' ', index=None)
            general.add_comments_to_start_of_file(os.path.join(self.spath, stem + suffix), comment)

        stem = pathlib.Path(basefilename).stem
        suffix = pathlib.Path(basefilename).suffix
        if len(liData) == 1:
            _save(stem, suffix, liData[0][1], liData[0][0])
        else:
            for i in range(len(liData)):
                _save(stem + str(i), suffix, liData[i][1], liData[i][0])

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
                scatt = M.fitness.theory()
                liData[i][1][data_column] = scatt
                i += 1
        else:
            self.problem.chisq()
            scatt = self.problem.fitness.theory()
            liData[0][1][data_column] = scatt

        return liData

    def fnSimulateErrorBars(self, simpar, liData):
        """
        Placeholder.
        """
        for i in range(len(liData)):
            liData[i][1]['dI'] = 0.1 * liData[i][1]['I']

        return liData
