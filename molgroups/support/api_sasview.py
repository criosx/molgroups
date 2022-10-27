from __future__ import print_function

import copy
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

    def fnSimulateDataPlusErrorBars(self, liData, diModelPars, simpar=None, basefilename='sim.dat', qmin=None,
                                    qmax=None, liConfigurations=None, lambda_min=0.1):
        # if qmin and qmax is not given, take from first data set
        if qmin is None:
            qmin = liData[0][1]['Q'].iloc[0]
        if qmax is None:
            qmax = liData[0][1]['Q'].iloc[-1]

        # Change q-range to ridculuously high value that allows full integration over 4Pi for
        # wavelengths down to 0.1Å. Ideally, one would choose exactly the right q-range for the integration for
        # the uncertainty calculation. This is difficult, as Lambda might vary with every configuration and a
        # new data simulation would be required. Instead, we simulate over a very large q-range and use only
        # a part of it during integration on the uncertainty function.
        # This is only needed for SANS

        liData = self.fnExtendQRange(liData=liData, qmin=0, qmax=4 * numpy.pi / lambda_min)
        self.fnSaveData(basefilename=basefilename, liData=liData)
        self.fnRestoreFit()

        # simulate data, works on sim.dat files
        liData = self.fnSimulateData(diModelPars, liData)
        # simulate error bars, works on sim.dat files
        liData = self.fnSimulateErrorBars(simpar, liData, liConfigurations)

        # recover requested q-range
        liData = self.fnExtendQRange(liData, qmin=qmin, qmax=qmax)

        return liData

    def fnSimulateData(self, diNewPars=None, liData=None, data_column='I'):

        if diNewPars is not None:
            self.fnUpdateModelPars(diNewPars)

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

    @staticmethod
    def fnSimulateErrorBars(simpar=None, liData=None, liConfigurations=None, percent_error=0.1):
        """
        Calculates uncertainties on a dataset liData given a list of configurations.

        liData is a list of experimental SANS curves: liData = [SANScurve1, SANScurve2, ...]
            Each SANScurve is a list of header and data SANScurve = [header, diData]
            diData is a dictionary of numpy 1D arrays with the keys: 'Q', 'I', 'dI', 'dQ'

        liConfiguration is a list of list: liConfiguration = [ConfSANScurve1 ,ConfSANScurve2, ...]
            Each SANScurve-specific configuration is a list of sub-configurations (i.e. different detector distances)
            ConfSANScurve1 = [diConf1, diConf2, ...]
            Each subconfiguration is a dictionary with the following possible entries. If no entry is provided,
            the indicated default values are used.

                TBD

                Comments: Neutron flux can be given per area or already times area, in which the preset of beam_area=1
                ensures the correct incident intensity.

        Comments:
            - theoretical SANS curves provided in liData are supposed to cover the entire possible q-range, because
              the 4Pi integrated intensity will be calculated from it. Q-values not covered by the provided instrumental
              configurations will be eliminated afterwards using a mask. The mask is either provided by the SASCALC
              instrument calculator, or it is approximated in this routine. If no configuration is provided, the entire
              curve will be returned using the percent_error attribute passed to this function.
        """
        def _add_default(dictonary, key_value_list):
            for pair in key_value_list:
                if not pair[0] in dictonary:
                    dictonary[pair[0]] = pair[1]

        def _configuration_defaults(configuration):
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
                ["detector_sample_distance", 400],
                ["detector_pixel_x", 128],
                ["detector_pixel_y", 128],
                ["detector_pixelsize_x", 0.00508],
                ["detector_pixelsize_y", 0.00508],
                ["time", 60],
                ["cuvette_thickness", 0.2],
                ["dark_count_f", 0],
                ["unit_solid_angle", 1],
                ["differential_cross_section_cuvette", 0],
                ["differential_cross_section_buffer", 0],
                ["dq_q", 0.06],
                ["lambda", 6],
                ["sascalc", False]
            ]
            # fill in defaults if no preset entries
            _add_default(configuration, kl_list)

            return configuration

        def _calc_configuration(Q, Iq, configuration):
            theta = 2 * numpy.arcsin(Q * configuration['lambda'] / 4 / numpy.pi)
            r = configuration['detector_sample_distance'] * numpy.tan(theta)
            # create an upper radius for the bin, which is typically the starting value of the next bin,
            # the last element being the edge case, otherwise a single roll would be sufficient
            r_gradient = numpy.gradient(r)
            r_gradient_roll = numpy.roll(r_gradient, -1)
            r_gradient_roll[-1] = r_gradient_roll[-2]
            r_plus = r + r_gradient + r_gradient_roll
            # omega = 2 * numpy.pi * (1 - numpy.cos(theta))
            omega = (Q * configuration['lambda'])**2 / (4 * numpy.pi)
            delta_omega = numpy.gradient(omega)
            Sigma = numpy.dot(Iq, delta_omega)

            if configuration["sascalc"]:
                # calculate n_cell(q) and dq/q using Jeff's sascalc implementation
                # li_n_cell, li_dq_q = jeff.sascalc(liData, configuration)
                # TODO: Update once Jeff's routine is available, until then use this placeholder:
                li_n_cell = numpy.full(Q, 1)
                li_dq_q = numpy.full(li_n_cell.shape, configuration['dq_q'])
            else:
                # estimate n_cell(q) and dq/q from configuration parameters
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

            return li_n_cell, li_dq_q, delta_omega, Sigma

        def _divide(a, b):
            return numpy.divide(a, b, out=numpy.zeros_like(a), where=b != 0)

        for dataset_n in range(len(liData)):
            # prepare data set
            dataset = liData[dataset_n]
            Q = liData[dataset_n][1]['Q'].to_numpy()
            Iq = liData[dataset_n][1]['I'].to_numpy()
            dQ = liData[dataset_n][1]['dQ'].to_numpy().fill(0)
            dIq = liData[dataset_n][1]['dI'].to_numpy().fill(0)

            if liConfigurations is None:
                # if no instrument configurations are given, calculate uncertainties as percent of data value
                # using the percent_error argument. Exit function afterwards.
                dIq = percent_error * Iq
            else:
                counts_sample = numpy.zeros_like(Iq)
                counts_cuvette = numpy.zeros_like(Iq)
                counts_dark = numpy.zeros_like(Iq)

                for configuration in liConfigurations[dataset_n]:

                    # set unused configuration keys to default values
                    configuration = _configuration_defaults(configuration)

                    li_n_cell, li_dq_q, d_omeg, Sigma_s_4pi = _calc_configuration(Q, Iq, configuration)

                    D = configuration["cuvette_thickness"]
                    neutron_intensity = configuration["neutron_flux"] * configuration["beam_area"]
                    neutron_intensity *= configuration["time"] * configuration["detector_efficiency"] * D

                    # Cuvette counts
                    # Not sure, how to treat empty cuvette scattering that takes place over a shorter path length
                    # through material than a buffer-filled cuvette. At the moment it is treated as it would occur over
                    # the entire length of the cuvette.
                    Sigma_c = configuration["differential_cross_section_cuvette"]
                    T_c = numpy.exp(-1 * Sigma_c * 4 * numpy.pi * D)
                    I_cuvette = neutron_intensity * Sigma_c

                    # Buffer counts
                    Sigma_b = configuration["differential_cross_section_buffer"]
                    T_b = numpy.exp(-1 * Sigma_b * 4 * numpy.pi * D)
                    I_buffer = neutron_intensity * Sigma_b

                    # Sample counts
                    Sigma_s = Iq
                    T_s = numpy.exp(-1 * Sigma_s_4pi * D)
                    I_sample = neutron_intensity * Sigma_s

                    # Dark counts
                    Sigma_d = configuration['dark_count_f']
                    I_dark = neutron_intensity * Sigma_d

                    counts_sample += li_n_cell * (I_sample + I_buffer + I_cuvette) * T_s * T_b * T_c + I_dark
                    counts_cuvette += li_n_cell * (I_buffer + I_cuvette) * T_b * T_c + I_dark
                    counts_dark += li_n_cell * I_dark

                    # calculate new relative uncertaint
                    # Approximat Poisson by Gaussian, can be changed
                    delta_I_s = (numpy.sqrt(counts_sample)/(T_s * T_b * T_c))**2
                    delta_I_s += (numpy.sqrt(counts_cuvette)/(T_b * T_c))**2
                    delta_I_s += (numpy.sqrt(counts_dark) * (1 / T_s / T_b / T_c) - (1 / T_b / T_c))**2
                    delta_I_s = numpy.sqrt(delta_I_s)
                    delta_I_I = delta_I_s / neutron_intensity
                    delta_I_I = _divide(delta_I_I, li_n_cell)
                    delta_I_I = _divide(delta_I_I, Iq)

                    # update current dI and dQ entries in data set with data from this frame
                    # old and new relative uncertainties
                    r_1 = _divide(dIq, Iq)
                    r_2 = _divide(delta_I_I, Iq)
                    # old and new equivalent count numbers (Gaussian amplitudes) associated with those rs
                    N_1 = _divide(numpy.ones_like(r_1), r_1*r_1)
                    N_2 = _divide(numpy.ones_like(r_2), r_2*r_2)
                    # old and new dQ
                    q_1 = dQ
                    q_2 = li_dq_q * Q

                    # update data frames
                    dataset[1]['dI'] = _divide(numpy.ones_like(N_1), N_1 + N_2) * Iq
                    dataset[1]['dQ'] = numpy.sqrt((N_1 * N_1 * q_1 * q_1 + N_2 * N_2 * q_2 * q_2) /
                                                  (N_1 + N_2) * (N_1 + N_2))

        return liData
