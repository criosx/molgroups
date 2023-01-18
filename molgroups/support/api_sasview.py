from __future__ import print_function

import copy
import numpy
import pandas
import os
import pathlib
import shapely.geometry

import numpy.random

import sasmodels.data
import bumps.curve

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

    def fnRestoreSmoothProfile(self, M):
        # TODO: Decide what and if to return SLD profile for Bumps fits
        # Returns nothing for the moment, here could be a SANS profile of any kind
        z, rho, irho = [], [], []
        return z, rho, irho

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

    def fnSimulateDataPlusErrorBars(self, liData, diModelPars, simpar=None, basefilename='sim.dat', qrange=None,
                                    liConfigurations=None, qmin=None, qmax=None, qrangefromfile=True, mode=None,
                                    lambda_min=0.1, t_total=None):
        # if qmin and qmax is not given, take from first data set
        # resolve amiguity / compatibility of the qrange parameter
        if qrangefromfile:
            # specified qmin or qmax values trump qrangefromfile
            if qmin is None:
                qmin = liData[0][1]['Q'].iloc[0]
            if qmax is None:
                qmax = liData[0][1]['Q'].iloc[-1]

        # Change q-range to logspacing. This is due to real SANS data having overlap regions between configuration
        # with irregular spacing. We need a regular spacing on which we will calculate SANS from different configs
        # and either average or keep data points with the same q separately

        numpoints = int((numpy.log10(qmax) - numpy.log10(qmin)) * 60)
        qvec = numpy.logspace(numpy.log10(qmin), numpy.log10(qmax), num=numpoints, endpoint=True)
        ivec = numpy.zeros_like(qvec)
        divec = numpy.zeros_like(qvec)
        dqvec = numpy.zeros_like(qvec)
        for dataset in liData:
            dataset[1] = pandas.DataFrame([qvec, ivec, divec, dqvec]).T
            dataset[1].columns = ['Q', 'I', 'dI', 'dQ']

        # Change q-range to ridculuously high value that allows full integration over 4Pi for
        # wavelengths down to 0.1Å. Ideally, one would choose exactly the right q-range for the integration for
        # the uncertainty calculation. This is difficult, as Lambda might vary with every configuration and a
        # new data simulation would be required. Instead, we simulate over a very large q-range and use only
        # a part of it during integration on the uncertainty function.
        # This is only needed for SANS

        liData = self.fnExtendQRange(liData=liData, qmin=0, qmax=4 * numpy.pi / lambda_min, conserve_dq_q=True)
        self.fnSaveData(basefilename=basefilename, liData=liData)
        self.fnRestoreFit()

        # simulate data, works on sim.dat files
        liData = self.fnSimulateData(diModelPars, liData)
        # simulate error bars, works on sim.dat files
        liData = self.fnSimulateErrorBars(simpar, liData=liData, liConfigurations=liConfigurations, t_total=t_total)
        # length of liData might have changed when different dQ/Q were calculated for the same q-values and
        # data was not averaged
        self.fnSaveData(basefilename=basefilename, liData=liData)
        self.fnRestoreFit()
        # Simulate twice. dQ/Q is only determined when simulating error bars because it needs the configuration
        # Might make sense to decouple dQ/Q from dIq, but might not be worth it. Neglect further iterations.
        liData = self.fnSimulateData(diModelPars, liData)

        # recover requested q-range
        if qmin is not None and qmax is not None:
            liData = self.fnExtendQRange(liData, qmin=qmin, qmax=qmax)
        # add random normalvariates
        for dataset in range(len(liData)):
            liData[dataset][1]['I'] = liData[dataset][1]['I'] + numpy.random.normal(size=len(liData[dataset][1]['dI']))\
                                     * liData[dataset][1]['dI']

        # Removes datapoints for which dI/I is zero, which indicates that no data was taken there
        # Let's also drop data points with more than 200% uncertainty
        for dataset_n in range(len(liData)):
            df = liData[dataset_n][1]
            df = df[df['dI'] > 0]
            liData[dataset_n][1] = df[df['dI'] < df['I']]

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
                if not isinstance(M.fitness, bumps.curve.Curve):
                    # ignore Curve models, as they are sometimes used as auxiliary functions that do not contain
                    # scattering data
                    liData[i][1][data_column] = scatt
                    i += 1
        else:
            self.problem.chisq()
            scatt = self.problem.fitness.theory()
            liData[0][1][data_column] = scatt

        return liData

    @staticmethod
    def fnSimulateErrorBars(simpar=None, liData=None, liConfigurations=None, percent_error=0.1, average=False,
                            t_total=None):
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

        average determines whether intensities at the same q-value but of different dQ/Q are averaged according to the
            underlying counting statistics or not. If not, separate data entries are maintained.

        t_total sets a total time spent over all configurations if it is not None

        Comments:
            - Neutron flux can be given per area or already times area, in which the preset of beam_area=1
              ensures the correct incident intensity.
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
                ["neutron_flux", 1e5],
                ["beam_area", 1],
                ["beam_center_x", 0],
                ["beam_center_y", 0],
                ["beamstop_diameter", 2.54],
                ["detector_efficiency", 1],
                ["sample_detector_distance", 400],
                ["source_sample_distance", 1614],
                ["source_aperture_radius", 0.715],
                ["sample_aperture_radius", 0.635],
                ["detector_pixel_x", 128],
                ["detector_pixel_y", 128],
                ["detector_pixelsize_x", 0.508],
                ["detector_pixelsize_y", 0.508],
                ["time", 60],
                ["cuvette_thickness", 0.2],
                ["dark_count_f", 0],
                ["unit_solid_angle", 1],
                ["differential_cross_section_cuvette", 0],
                ["differential_cross_section_buffer", 0],
                ["dq_q", None],
                ["lambda", 6],
                ["dlambda_lambda", 0.125],
                ["sascalc", False]
            ]
            # fill in defaults if no preset entries
            _add_default(configuration, kl_list)

            return configuration

        def _calc_configuration(Q, Iq, configuration):
            # 4PI integration is valid only until q_max given by lambda
            q_max = 4 * numpy.pi / configuration['lambda']
            # index for closest q-value to q_max
            idx = numpy.abs(Q-q_max).argmin()
            if Q[idx] > q_max and idx > 0:
                idx -= 1
            initial_q_length = len(Q)
            # limit Q and Iq to allowed range
            Q = Q[0:idx+1]
            Iq = Iq[0:idx+1]
            omega = (Q * configuration['lambda'])**2 / (4 * numpy.pi)
            omega_roll = numpy.roll(omega, -1)
            omega_roll[-1] = omega[-1] + numpy.gradient(omega)[-1]
            delta_omega = omega_roll - omega
            # dot product only over this limited index range
            Sigma = numpy.dot(Iq, delta_omega)

            # now limit omega to 2Pi for SANS geometry reasons (no backscattering)
            q_max2 = numpy.sqrt(8) * numpy.pi / configuration['lambda']
            idx = numpy.abs(Q-q_max2).argmin()
            if Q[idx] > q_max2 and idx > 0:
                idx -= 1
            Q = Q[0:idx+1]
            # Iq = Iq[0:idx+1]

            theta = 2 * numpy.arcsin(Q * configuration['lambda'] / 4 / numpy.pi)
            r = configuration['sample_detector_distance'] * numpy.tan(theta)
            # create an upper radius for the bin, which is typically the starting value of the next bin,
            r_gradient = numpy.gradient(r)
            r_plus = numpy.roll(r, -1)
            r_plus[-1] = r[-1] + r_gradient[-1]
            # omega = 2 * numpy.pi * (1 - numpy.cos(theta))
            omega = (Q * configuration['lambda'])**2 / (4 * numpy.pi)
            delta_omega = numpy.gradient(omega)

            if configuration["sascalc"]:
                # calculate n_cell(q) and dq/q using Jeff's sascalc implementation
                # li_n_cell, li_dq_q = jeff.sascalc(liData, configuration)
                # TODO: Update once Jeff's routine is available, until then use this placeholder:
                n_cell = delta_omega * 1
                dq_q = numpy.full(Q.shape, 0.15)
            else:
                # estimate n_cell(q) and dq/q from configuration parameters
                # draw detector dimensions and solid angle projections in shapely
                dimx = configuration['detector_pixel_x'] * configuration['detector_pixelsize_x']
                dimy = configuration['detector_pixel_y'] * configuration['detector_pixelsize_y']
                detector = shapely.geometry.Polygon([(-dimx/2, -dimy/2), (dimx/2, -dimy/2), (dimx/2, dimy/2),
                                                     (-dimx/2, dimy/2)])
                cx = configuration["beam_center_x"]
                cy = configuration["beam_center_y"]
                circle_inner = [shapely.geometry.Point(cx, cy).buffer(radius) for radius in numpy.nditer(r)]
                circle_outer = [shapely.geometry.Point(cx, cy).buffer(radius) for radius in numpy.nditer(r_plus)]
                omega_circle = [circle_outer[i].difference(circle_inner[i]) for i in range(len(circle_outer))]
                omega_circle_area = numpy.array([omega_circle[i].area for i in range(len(omega_circle))])

                beamstop_circle = shapely.geometry.Point(cx, cy).buffer(configuration["beamstop_diameter"]/2)
                detector = detector.difference(beamstop_circle)

                area_detected = numpy.array([omega_circle[i].intersection(detector).area
                                             for i in range(len(omega_circle))])

                # n_cell(q) is the entire solid angle of this q-value times the fraction that is detected
                n_cell = delta_omega * (area_detected / omega_circle_area)

                if configuration['dq_q'] is not None:
                    dqq = numpy.full(Q.shape, configuration['dq_q'])
                else:
                    # calculate dq_q based on SANS toolbox page 154, equation (30), approximated for sigma_x = sigma_y
                    L1 = configuration["source_sample_distance"]
                    L2 = configuration["sample_detector_distance"]
                    R1 = configuration["source_aperture_radius"]
                    R2 = configuration["sample_aperture_radius"]
                    dx = configuration["detector_pixelsize_x"]
                    dy = configuration["detector_pixelsize_y"]
                    dL_L = configuration["dlambda_lambda"]
                    L = configuration["lambda"]

                    dqx2 = (L2*R1/L1/2)**2+((L1+L2)*R2/L1/2)**2+(1/3)*(dx/2)**2
                    dqx2 *= (numpy.pi*2/L/L2)**2
                    dqx2 += (1/6)*dL_L**2*Q*Q
                    A = L2 * (L1 + L2) * 3.073e-9
                    dqy2 = (L2*R1/L1/2)**2+((L1+L2)*R2/L1/2)**2+(1/3)*(dy/2)**2+A**2*L**4*(2/3)*dL_L**2
                    dqy2 *= (numpy.pi*2/L/L2)**2
                    dqy2 += (1/6)*dL_L**2*Q*Q

                    dqq = _divide(numpy.sqrt(dqx2 + dqy2)/numpy.sqrt(2), Q)

                # pad numpy arrays with zeros back to original length
                append_array = numpy.zeros(initial_q_length - len(n_cell))
                n_cell = numpy.append(n_cell, append_array)
                dq_q = numpy.append(dqq, append_array)
                delta_omega = numpy.append(delta_omega, append_array)

            return n_cell, dq_q, delta_omega, Sigma

        def _divide(a, b):
            return numpy.divide(a, b, out=numpy.zeros_like(a), where=b != 0)

        # Prepare configurations.
        if liConfigurations is not None:
            # check if a single set of configurations is provided for more than one dataset and
            # multiply the configurations accordingly
            if len(liConfigurations) < len(liData):
                if len(liConfigurations) == 1:
                    for _ in range(len(liData)-len(liConfigurations)):
                        liConfigurations.append(copy.deepcopy(liConfigurations[0]))
                else:
                    raise Exception('Number of configuration sets should be one or must match the number of datasets!')

            # supplement configurations with default data and adjust total time, if required
            actual_total_time = 0
            for dataset_n in range(len(liData)):
                for configuration in liConfigurations[dataset_n]:
                    # set unused configuration keys to default values
                    configuration = _configuration_defaults(configuration)
                    actual_total_time += configuration['time']

            # adjust measurement times if t_total is set
            if t_total is not None:
                if actual_total_time != 0:
                    time_factor = float(t_total) / float(actual_total_time)
                else:
                    time_factor = 0
                for dataset_n in range(len(liData)):
                    for configuration in liConfigurations[dataset_n]:
                        configuration['time'] *= time_factor

        # compute dI/I and dQ/Q
        for dataset_n in range(len(liData)):
            # prepare data set
            dataset = liData[dataset_n]
            Q = liData[dataset_n][1]['Q'].to_numpy(copy=True)
            Iq = liData[dataset_n][1]['I'].to_numpy(copy=True)
            dQ = liData[dataset_n][1]['dQ'].to_numpy(copy=True)
            dIq = liData[dataset_n][1]['dI'].to_numpy(copy=True)

            # in case no average = False
            Q_append = numpy.empty(0)
            Iq_append = numpy.empty(0)
            dQ_append = numpy.empty(0)
            dIq_append = numpy.empty(0)

            if liConfigurations is None:
                # if no instrument configurations are given, calculate uncertainties as percent of data value
                # using the percent_error argument. Exit function afterwards.
                dIq = percent_error * Iq
                dataset[1]['dI'] = dIq
            else:
                dQ.fill(0)
                dIq.fill(0)

                for configuration in liConfigurations[dataset_n]:

                    li_n_cell, li_dq_q, d_omeg, Sigma_s_4pi = _calc_configuration(Q, Iq, configuration)

                    D = configuration["cuvette_thickness"]
                    neutron_intensity = configuration["neutron_flux"] * configuration["beam_area"]
                    neutron_intensity *= configuration["time"] * configuration["detector_efficiency"]

                    # avoid zero intensity, as they produce data points that bumps cannot handle
                    if neutron_intensity == 0:
                        neutron_intensity = 1e-12

                    # Cuvette counts
                    # Not sure, how to treat empty cuvette scattering that takes place over a shorter path length
                    # through material than a buffer-filled cuvette. At the moment it is treated as it would occur over
                    # the entire length of the cuvette. See theory document.

                    Sigma_c = configuration["differential_cross_section_cuvette"]
                    Sigma_c_4pi = Sigma_c * 4 * numpy.pi
                    Sigma_b = configuration["differential_cross_section_buffer"]
                    Sigma_b_4pi = Sigma_b * 4 * numpy.pi
                    Sigma_s = Iq
                    Sigma_T_4pi = Sigma_s_4pi + Sigma_b_4pi + Sigma_c_4pi

                    T_T = numpy.exp(-1 * Sigma_T_4pi * D)
                    T_bc = numpy.exp(-1 * (Sigma_b_4pi + Sigma_b_4pi) * D)
                    T_c = numpy.exp(-1 * Sigma_c_4pi * D)
                    T_b = numpy.exp(-1 * Sigma_b_4pi * D)
                    T_s = numpy.exp(-1 * Sigma_s_4pi * D)

                    Sigma_bmT = Sigma_b_4pi - Sigma_T_4pi
                    Sigma_smT = Sigma_s_4pi - Sigma_T_4pi
                    T_1mT = 1 - T_T
                    T_1mb = 1 - T_b

                    # sample intensity
                    I_sample = neutron_intensity * (Sigma_s / Sigma_bmT) * (T_T - T_b)
                    # buffer intensity
                    I_buffer = neutron_intensity * (1 / (Sigma_bmT * Sigma_T_4pi))
                    I_buffer *= (Sigma_b_4pi**2 * T_1mT + Sigma_b_4pi * Sigma_smT * T_1mT - Sigma_s_4pi * Sigma_s_4pi \
                                * T_1mb)
                    I_buffer /= (4 * numpy.pi)
                    # self-shielding correction (Barker, Mildner, J. Appl. Cr., 2015)
                    I_buffer *= (1 + 0.625 * Sigma_T_4pi * D)
                    # cuvette intensity
                    I_cuvette = neutron_intensity * Sigma_c_4pi / Sigma_T_4pi * T_1mT / (4 * numpy.pi)
                    # Dark count intensity
                    Sigma_d = configuration['dark_count_f']
                    I_dark = neutron_intensity * Sigma_d

                    counts_sample = li_n_cell * (I_sample + I_buffer + I_cuvette + I_dark)
                    counts_cuvette = li_n_cell * (I_buffer + I_cuvette + I_dark)
                    counts_dark = li_n_cell * I_dark

                    # calculate new relative uncertainty, approximate Poisson by Gaussian
                    # make sure that the uncertainty is not larger than the count (counts below one).
                    sqrt_counts_sample = numpy.minimum(numpy.sqrt(counts_sample), counts_sample)
                    sqrt_counts_cuvette = numpy.minimum(numpy.sqrt(counts_cuvette), counts_cuvette)
                    sqrt_counts_dark = numpy.minimum(numpy.sqrt(counts_dark), counts_dark)

                    # calculate new relative uncertainty
                    # Approximate Poisson by Gaussian
                    delta_I_s = (sqrt_counts_sample / T_T)**2
                    delta_I_s += (sqrt_counts_cuvette / T_bc)**2
                    delta_I_s += (sqrt_counts_dark * ((1 / T_T) - (1 / T_bc)))**2
                    delta_I_s = numpy.sqrt(delta_I_s)
                    delta_I_I = delta_I_s / neutron_intensity
                    delta_I_I = _divide(delta_I_I, li_n_cell)
                    delta_I_I = _divide(delta_I_I, D)
                    delta_I_I = _divide(delta_I_I, Iq)

                    if average:
                        # update current dI and dQ entries in data set with data from this frame
                        # old and new relative uncertainties are averaged according to underlying counting statistics
                        r_1 = _divide(dIq, Iq)
                        r_2 = delta_I_I
                        # old and new equivalent count numbers (Gaussian amplitudes) associated with those rs
                        N_1 = _divide(numpy.ones_like(r_1), r_1*r_1)
                        N_2 = _divide(numpy.ones_like(r_2), r_2*r_2)
                        # old and new dQ
                        q_1 = dQ
                        q_2 = li_dq_q * Q

                        # update data frames
                        dIq = _divide(numpy.ones_like(N_1), numpy.sqrt(N_1 + N_2)) * Iq
                        dQ = numpy.sqrt(_divide((N_1 * N_1 * q_1 * q_1 + N_2 * N_2 * q_2 * q_2),
                                                (N_1 * N_1 + N_2 * N_2)))
                    else:
                        Q_append = numpy.append(Q_append, Q)
                        Iq_append = numpy.append(Iq_append, Iq)
                        dIq_append = numpy.append(dIq_append, delta_I_I * Iq)
                        dQ_append = numpy.append(dQ_append, li_dq_q * Q)

                if average:
                    dataset[1]['dI'] = dIq
                    dataset[1]['dQ'] = dQ
                else:
                    dataset[1] = pandas.DataFrame([Q_append, Iq_append, dIq_append, dQ_append])
                    dataset[1] = dataset[1].T
                    dataset[1].columns = ['Q', 'I', 'dI', 'dQ']
                    dataset[1] = dataset[1].sort_values(by=['Q'])

        return liData






