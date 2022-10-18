from __future__ import print_function

import numpy
import pandas
import os
import pathlib

import shapely.geometry

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

    def fnSimulateErrorBars(self, simpar=None, liData=None, liConfigurations=None, percent_error=0.1):
        """
        Placeholder.
        Neutron flux can be given per area or already times area, in which the preset of beam_area=1
        ensures the correct incident intensity.
        """
        def _add_default(dictonary, key_value_list):
            for pair in key_value_list:
                if not pair[0] in dictonary:
                    dictonary[pair[0]] = pair[1]

        # if no instrument configurations are given, calculate uncertainties as percent of data value
        # using the percent_error argument. Exit function afterwards.
        if liConfigurations is None:
            for i in range(len(liData)):
                liData[i][1]['dI'] = percent_error * liData[i][1]['I']
            return liData

        # list of [n_cell(q), dq/q] tuples for every data set
        limask = []

        for dataset_n in range(len(liData)):
            limask.append([])

            for configuration in liConfigurations:
                # fill in defaults if no preset entries
                if configuration is None:
                    configuration = {}
                kl_list = [
                    ["pre", 1],
                    ["neutron_flux", 1e5],
                    ["beam_area", 1],
                    ["beam_center_x", 64.5],
                    ["beam_center_y", 64.5],
                    ["beamstop_diameter", 25.4],
                    ["detector_efficiency", 1],
                    ["detector_sample_distance", 4000],
                    ["detector_pixel_x", 128],
                    ["detector_pixel_y", 128],
                    ["detector_pixelsize_x", 0.00508],
                    ["detector_pixelsize_y", 0.00508],
                    ["time", 60],
                    ["sample_trans", 1],
                    ["cuvette_trans", 1],
                    ["dark_count_f", 0],
                    ["unit_solid_angle", 1],
                    ["cross_section_cuvette", 0],
                    ["cross_section_buffer", 0],
                    ["q_min", 0.001],
                    ["q_max", 0.4],
                    ["dq_q", 0.06],
                    ["lambda", 6],
                    ["sascalc", False]
                ]
                _add_default(configuration, kl_list)
                # provide compatibility/redundancy between the usages of q_max and qrange
                if 'qrange' in configuration:
                    configuration['q_max'] = configuration['qrange']

                if configuration["sascalc"]:
                    # calculate n_cell(q) and dq/q using Jeff's sascalc implementation
                    # li_n_cell, li_dq_q = jeff.sascalc(liData, configuration)
                    # TODO: Update once Jeff's routine is available, until then use this placeholder:
                    li_n_cell = numpy.full(liData[dataset_n][1]['Q'], 1)
                    li_dq_q = numpy.full(li_n_cell.shape, configuration['dq_q'])

                else:
                    # estimate n_cell(q) and dq/q from configuration parameters

                    q = liData[dataset_n][1]['Q']
                    theta = 2 * numpy.arcsin(q * configuration['lambda'] / 4 / numpy.pi)
                    r = configuration['detector_sample_distance'] * numpy.tan(theta)

                    # create an upper radius for the bin, which is typically the starting value of the next bin,
                    # the last element being the edge case
                    r_gradient = numpy.gradient(r)
                    r_gradient_roll = numpy.roll(r_gradient, -1)
                    r_gradient_roll[-1] = r_gradient_roll[-2]
                    r_plus = r + r_gradient + r_gradient_roll

                    omega = 2 * numpy.pi * (1 - numpy.cos(theta))
                    omega_gradient = numpy.gradient(omega)
                    omega_gradient_roll = numpy.roll(omega_gradient, -1)
                    omega_gradient_roll[-1] = omega_gradient_roll[-2]
                    delta_omega = omega_gradient + omega_gradient_roll

                    # draw detector dimensions and solid angle projections in shapely
                    dimx = configuration['detector_pixel_x'] * configuration['detector_pixelsize_x']
                    dimy = configuration['detector_pixel_y'] * configuration['detector_pixelsize_y']
                    detector = shapely.geometry.Polygon([(-dimx/2, -dimy/2), (dimx/2, -dimy/2), (dimx/2, dimy/2),
                                                         (-dimx/2, dimy/2)])
                    cx = configuration["beam_center_x"]
                    cy = configuration["beam_center_y"]
                    circle_inner = [shapely.geometry.Point(cx, cy).buffer(radius) for radius in r]
                    circle_outer = [shapely.geometry.Point(cx, cy).buffer(radius) for radius in r_plus]
                    omega_circle = [circle_outer[i].difference(circle_inner[i]) for i in range(len(circle_outer))]
                    omega_circle_area = numpy.array([omega_circle[i].area for i in range(len(omega_circle))])

                    beamstop_circle = shapely.geometry.Point(cx, cy).buffer(configuration["beamstop_diameter"])
                    detector = detector.difference(beamstop_circle)

                    area_detected = numpy.array([omega_circle[i].intersection(detector).area
                                                 for i in range(len(omega_circle))])

                    # n_cell(q) is the entire solid angle of this q-value times the fraction that is detected
                    li_n_cell = delta_omega * (area_detected / omega_circle_area)
                    li_dq_q = numpy.full(li_n_cell.shape, configuration['dq_q'])

                limask[dataset_n].append([li_n_cell, li_dq_q])

        return liData
