import itertools

import numpy
import numpy as np
import os
import pickle
import pandas
import sys
import matplotlib.pyplot as plt
import shutil

from numpy.random import permutation
from os import path, mkdir
from math import fabs, pow, floor, ceil
from subprocess import Popen
from time import sleep
from IPython.display import clear_output
from gpcam.autonomous_experimenter import AutonomousExperimenterGP

from molgroups.support import molstat
from molgroups.infotheory import MVN
from molgroups.infotheory import GMM


# static methods
def average(alist):
    while True:
        s = np.std(alist)
        result = np.mean(alist)
        maxs = 0
        val = 0
        for element in alist:
            s2 = fabs(result - element)
            if s2 > maxs:
                maxs = s2
                val = element
        if maxs > 3 * s:
            alist.remove(val)
        else:
            break
    return result, s


# multidimensional convolution of nearest horizontal, vertical, and diagonal neighbors but not the center
def convolute(a):
    conv_arr = np.zeros(a.shape)
    offsetshape = tuple([3] * len(a.shape))
    offsetparent = np.zeros(offsetshape)
    it = np.nditer(conv_arr, flags=['multi_index'])
    while not it.finished:
        itindex = tuple(it.multi_index[i] for i in range(len(a.shape)))
        it2 = np.nditer(offsetparent, flags=['multi_index'])
        result = 0.0
        counter = 0.0
        while not it2.finished:
            itindex2 = tuple(it2.multi_index[i] - 1 for i in range(len(a.shape)))
            allzero = True
            for element in itindex2:
                if element != 0:
                    allzero = False
            if not allzero:
                index = tuple(itindex[i] + itindex2[i] for i in range(len(itindex)))
                inarray = True
                for i, element in enumerate(index):
                    if element < 0 or element >= a.shape[i]:
                        inarray = False
                if inarray:
                    result += a[index]
                    counter += 1.0
            it2.iternext()

        conv_arr[itindex] = result / counter
        it.iternext()

    return conv_arr


def nice_interval(start=0, stop=1, step=None, numsteps=10):
    """
    Paul Kienzle's method to obtain nicely spaced intervals for plot axes
    """

    if step is None:
        step = (stop-start)/numsteps

    sign = 1.0
    if step < 0:
        sign = -1.0
        step = step * (-1)

    if step == 0:
        return [start, stop+1]

    exponent = floor(np.log(step)/np.log(10))
    mantisse = step / pow(10, exponent)

    new_mantisse = 1
    if fabs(mantisse-2) < fabs(mantisse-new_mantisse):
        new_mantisse = 2
    if fabs(mantisse-5) < fabs(mantisse-new_mantisse):
        new_mantisse = 5
    if fabs(mantisse-10) < fabs(mantisse-new_mantisse):
        new_mantisse = 10

    new_step = sign * new_mantisse * pow(10, exponent)
    new_start = floor(start/new_step) * new_step
    new_stop = ceil(stop/new_step) * new_step

    return np.arange(new_start, new_stop+0.05*new_step, new_step)


def rm_file(filename):
    try:
        os.remove(filename)
    except OSError:
        pass


def running_mean(current_mean, n, new_point):
    # see Tony Finch, Incremental calculation of weighted mean and variance
    return current_mean * (n - 1) / n + new_point * (1 / n)


def running_sqstd(current_sqstd, n, new_point, previous_mean, current_mean):
    # see Tony Finch, Incremental calculation of weighted mean and variance
    return (current_sqstd * (n - 1) + (new_point - previous_mean) * (new_point - current_mean)) / n


def save_plot_1d(x, y, dy=None, xlabel='', ylabel='', color='blue', filename="plot", ymin=None, ymax=None, levels=5,
                 niceticks=False, keep_plots=False):
    import matplotlib.pyplot as plt
    import matplotlib

    if ymin is None:
        ymin = np.amin(y)
    if ymax is None:
        ymax = np.amax(y)

    font = {'family': 'sans-serif', 'weight': '200', 'size': 14}
    matplotlib.rc('font', **font)

    fig, ax = plt.subplots()
    if dy is None:
        ax.plot(x, y, color=color)
    else:
        ax.errorbar(x, y, dy, color=color)
    if niceticks:
        bounds = nice_interval(start=ymin, stop=ymax, numsteps=levels)
        ax.set_yticks(bounds)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.ticklabel_format(scilimits=(-3, 3), useMathText=True)

    plt.tight_layout()

    if keep_plots:
        i = 0
        while os.path.isfile(filename + str(i) + '.png'):
            i += 1
        filename = filename + str(i)

    plt.savefig(filename + '.pdf')
    plt.savefig(filename + '.png')
    plt.close("all")


def save_plot_2d(x, y, z, xlabel, ylabel, color, filename='plot', zmin=None, zmax=None, levels=20, mark_maximum=False,
                 keep_plots=False, support_points=None):
    import matplotlib.pyplot as plt
    import matplotlib

    if zmin is None:
        zmin = np.amin(z)
    if zmax is None:
        zmax = np.amax(z)
    bounds = nice_interval(start=zmin, stop=zmax, numsteps=levels)

    font = {'family': 'sans-serif', 'weight': '200', 'size': 14}
    matplotlib.rc('font', **font)

    fig, ax = plt.subplots()
    cs = ax.contourf(x, y, z, cmap=color, vmin=bounds[0], vmax=bounds[-1], levels=bounds, extend='both')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.ticklabel_format(scilimits=(-3, 3), useMathText=True)
    fig.colorbar(cs)

    if support_points is not None:
        plt.scatter(support_points[:, 1], support_points[:, 0], s=10, c='k')

    if mark_maximum:
        index = np.unravel_index(z.argmax(), z.shape)
        plt.text(x[index[1]], y[index[0]], 'x', horizontalalignment='center', verticalalignment='center')

    plt.tight_layout()

    if keep_plots:
        i = 0
        while os.path.isfile(filename + str(i) + '.png'):
            i += 1
        filename = filename + str(i)

    plt.savefig(filename + '.pdf')
    plt.savefig(filename + '.png')
    plt.close("all")


# calculates entropy while varying a set of parameters in parlist and keeping others fixed as specified in simpar.dat
# requires a compiled and ready to go fit whose fit parameters are modified and fixed
# avoid_symmetric prevents calculating symmetry-related results by enforcing the indices of varied parameters
# are an ordered list
class Entropy:

    def __init__(self, fitsource, spath, mcmcpath, runfile, mcmcburn=16000, mcmcsteps=5000, deldir=True,
                 convergence=2.0, miniter=1, mode='water', background_rule=None, bFetchMode=False, bClusterMode=False,
                 calc_symmetric=True, upper_info_plotlevel=None, plotlimits_filename='', slurmscript='',
                 configuration=None, optimizer='grid', keep_plots=False, show_support_points=False, qmin=None,
                 qmax=None, qrangefromfile=False, t_total=None, jupyter_clear_output=False, gpcam_iterations=50,
                 gpcam_init_dataset_size=20, gpcam_step=None, acq_func="variance", fitter='MCMC', remove_fit_dir=True,
                 lm_iterations=3):

        self.fitsource = fitsource
        self.spath = spath
        self.mcmcpath = mcmcpath
        self.runfile = runfile
        self.mcmcburn = mcmcburn
        self.mcmcsteps = mcmcsteps
        self.deldir = deldir
        self.convergence = convergence
        self.miniter = miniter
        self.mode = mode
        self.background_rule = background_rule
        self.bFetchMode = bFetchMode
        self.bClusterMode = bClusterMode
        self.calc_symmetric = calc_symmetric
        self.upper_info_plotlevel = upper_info_plotlevel
        self.plotlimits_filename = plotlimits_filename
        self.slurmscript = slurmscript
        self.configuration = configuration
        self.optimizer = optimizer
        self.gpiteration = 0
        self.keep_plots = keep_plots
        self.show_support_points = show_support_points
        self.qmin = qmin
        self.qmax = qmax
        self.qrangefromfile = qrangefromfile
        self.t_total = t_total
        self.jupyter_clear_output = jupyter_clear_output
        self.gpcam_iterations = gpcam_iterations
        self.gpcam_init_dataset_size = gpcam_init_dataset_size
        self.gpcam_step = gpcam_step
        self.acq_func = acq_func
        self.fitter = fitter
        self.remove_fit_dir = remove_fit_dir
        self.lm_iterations = lm_iterations

        self.my_ae = None

        self.molstat = molstat.CMolStat(fitsource=fitsource, spath=spath, mcmcpath=mcmcpath, runfile=runfile)

        # Parse parameters from entropypar.dat
        header_names = ['type', 'dataset', 'configuration', 'par', 'value', 'l_fit', 'u_fit', 'l_sim', 'u_sim',
                        'step_sim']
        self.allpar = pandas.read_csv('entropypar.dat', sep='\s+', header=None, names=header_names,
                                      skip_blank_lines=True, comment='#')

        # define unique names, since instrument parameters might have the same name for different datasets and
        # configurations
        self.allpar['unique_name'] = ''
        for i in range(len(self.allpar['par'])):
            if 'n' in self.allpar['type'].iloc[i]:
                datastring = str(self.allpar['dataset'].iloc[i])
                if datastring == '*':
                    datastring = 'x'
                configstring = str(self.allpar['configuration'].iloc[i])
                if configstring == '*':
                    configstring = 'x'
                unique_name = self.allpar['par'].iloc[i] + '_' + datastring + '_' + configstring
            else:
                unique_name = self.allpar['par'].iloc[i]
            self.allpar.iloc[i, self.allpar.columns.get_loc("unique_name")] = unique_name

        # identify dependent (a), independent (b), and non-parameters in simpar.dat for the calculation of p(a|b,y)
        # later on. It is assumed that all parameters in setup.cc are also specified in simpar.dat in exactly the same
        # order. This might have to be looked at in the future.
        # keys: i: independent (nuisance parameter), d: dependent (parameter of interest), n or otherwise none
        # want to calculate H(d|i,y)

        self.dependent_parameters = []
        self.independent_parameters = []
        self.parlist = []
        self.joblist = []
        self.jobmax = 10

        # fetch mode and cluster mode are exclusive
        if self.bFetchMode and self.bClusterMode:
            self.bClusterMode = False

        i = 0
        for row in self.allpar.itertuples():
            if row.type == 'i' or row.type == 'fi':
                self.independent_parameters.append(row.par)
                self.parlist.append(row.par)
                i += 1
            elif row.type == 'd' or row.type == 'fd':
                self.dependent_parameters.append(row.par)
                self.parlist.append(row.par)
                i += 1

        # only those parameters that will be varied
        self.steppar = self.allpar.dropna(axis=0)

        # create data frame for simpar.dat needed by the data simulation routines
        # non-parameters such as qrange and prefactor will be included in simpar, but eventually ignored
        # when simulating the scattering, as they will find no counterpart in the model
        self.simpar = self.allpar.loc[:, ['par', 'value', 'dataset', 'configuration', 'unique_name']]

        self.steplist = []
        self.axes = []
        for row in self.steppar.itertuples():
            steps = int((row.u_sim - row.l_sim) / row.step_sim) + 1
            self.steplist.append(steps)
            axis = []
            for i in range(steps):
                axis.append(row.l_sim + i * row.step_sim)
            self.axes.append(axis)

        self.priorentropy, self.priorentropy_marginal = self.calc_prior()

        if optimizer == 'grid':
            self.results_mvn = np.full(self.steplist, self.priorentropy)
            self.results_gmm = np.full(self.steplist, self.priorentropy)
            self.results_mvn_marginal = np.full(self.steplist, self.priorentropy_marginal)
            self.results_gmm_marginal = np.full(self.steplist, self.priorentropy_marginal)
            self.n_mvn = np.zeros(self.results_mvn.shape)
            self.n_gmm = np.zeros(self.results_gmm.shape)
            self.n_mvn_marginal = np.zeros(self.results_mvn_marginal.shape)
            self.n_gmm_marginal = np.zeros(self.results_gmm_marginal.shape)
            self.sqstd_mvn = np.zeros(self.results_mvn.shape)
            self.sqstd_gmm = np.zeros(self.results_gmm.shape)
            self.sqstd_mvn_marginal = np.zeros(self.results_mvn_marginal.shape)
            self.sqstd_gmm_marginal = np.zeros(self.results_gmm_marginal.shape)
            self.prediction_gpcam = np.zeros(self.results_gmm_marginal.shape)
            self.par_median = np.zeros((len(self.parlist),) + self.results_mvn.shape)
            self.par_std = np.zeros((len(self.parlist),) + self.results_mvn.shape)

            if path.isfile(path.join(spath, 'results', 'MVN_entropy.npy')):
                self.load_results(spath)

        elif optimizer == 'gpcam' or optimizer == 'gpCAM':
            if path.isfile(path.join(spath, 'results', 'gpCAMstream.pkl')):
                with open(path.join(spath, 'results', 'gpCAMstream.pkl'), 'rb') as file:
                    self.gpCAMstream = pickle.load(file)
                    self.gpiteration = len(self.gpCAMstream['position'])
            else:
                self.gpCAMstream = {'position': [], 'value': [], 'variance': []}

    def calc_entropy(self, molstat=None, cov=False):

        if cov:
            # use covariance matrix indicates that LM fit is present with active molstat, problem and state
            parnames_dict_keys = molstat.Interactor.problem.labels()
        else:
            # MCMC fit
            if molstat is None:
                molstat = self.molstat
            points, parnames_dict_keys, logp = molstat.Interactor.fnLoadMCMCResults()

        independent_pars = []
        dependent_pars = []
        parnames = [key for key in parnames_dict_keys]
        for index in range(len(parnames)):
            if parnames[index] in self.independent_parameters:
                independent_pars.append(index)
            else:
                dependent_pars.append(index)

        if cov:
            mvnentropy = MVN.MVNEntropy(cov=molstat.Interactor.problem.cov())
            gmm_entropy = mvn_entropy = mvnentropy.entropy()
            gmm_entropy_marginal = mvn_entropy_marginal = mvnentropy.marginal_entropy(independent_pars=independent_pars)

            points_median = molstat.Interactor.problem.getp()
            points_std = molstat.Interactor.problem.stderr()[0]

        else:
            N_entropy = 10000  # was 10000
            N_norm = 10000  # was 2500

            # Do statistics over points
            points_median = np.median(points, axis=0)
            points_std = np.std(points, axis=0)

            # Use a random subset to estimate density
            if N_norm >= len(logp):
                norm_points = points
            else:
                idx = permutation(len(points))[:N_norm]
                norm_points = points[idx]

            mvnentropy = MVN.MVNEntropy(x=norm_points)
            mvn_entropy = mvnentropy.entropy()
            mvn_entropy_marginal = mvnentropy.marginal_entropy(independent_pars=independent_pars)

            gmmentropy = GMM.GMMEntropy(norm_points)
            gmm_entropy = gmmentropy.entropy(N_entropy)
            gmmentropymarginal = GMM.GMMEntropy(np.delete(norm_points, independent_pars, 1))
            gmm_entropy_marginal = gmmentropymarginal.entropy(N_entropy)

        '''
        # This is the original implementation of the KDN entropy from Kramer et al.
        # It is not very stable and replaced by the GMM entropy instead.

        # Use a different subset to estimate the scale factor between density
        # and logp.
        if N_entropy >= len(logp):
            entropy_points, eval_logp = points, logp
        else:
            idx = permutation(len(points))[:N_entropy]
            entropy_points, eval_logp = points[idx], logp[idx]

        # Calculate Kramer Normalized Entropy
        gmmrho = gmmentropy.score_samples(entropy_points)
        frac = exp(eval_logp) / exp(gmmrho)
        n_est, n_err = mean(frac), std(frac)
        s_est = log(n_est) - mean(eval_logp)
        # s_err = n_err/n_est
        # print(n_est, n_err, s_est/LN2, s_err/LN2)
        # print(np.median(frac), log(np.median(frac))/LN2, log(n_est)/LN2)
        kdn_entropy = s_est / LN2

        dependent_points = entropy_points[:, dependent_pars]
        kdn_entropy_marginal = (-1) * np.mean(gmmentropymarginal.score_samples(dependent_points)) / LN2
        '''

        # return MVN entropy, GMM entropy, conditional MVN entropy, conditional GMM entropy
        return mvn_entropy, gmm_entropy, mvn_entropy_marginal, gmm_entropy_marginal, points_median, points_std, parnames

    def calc_entropy_for_iteration(self, molstat_iter, itindex=None, cov=False):
        # calculate entropy, dependent parameters == parameters of interest
        # independent parameters == nuisance parameters
        mvn, gmm, mvn_marginal, gmm_marginal, points_median, points_std, parnames = \
            self.calc_entropy(molstat_iter, cov=cov)

        if mvn_marginal is None or gmm_marginal is None:
            bValidResult = False
        else:
            bValidResult = (self.priorentropy_marginal - gmm_marginal > (-0.5) * len(self.dependent_parameters)) and \
                           (self.priorentropy - gmm > (-0.5) * len(self.parlist))

        # no special treatment for first entry necessary, algorithm catches this
        if itindex is not None:
            if bValidResult:
                self.gridsearch_writeout_result(itindex, mvn, gmm, mvn_marginal, gmm_marginal, points_median,
                                                points_std, parnames)
            # save results for every iteration
            self.save_results(self.spath)

        return gmm_marginal

    # calculates prior entropy
    def calc_prior(self):
        priorentropy = 0
        priorentropy_marginal = 0
        for row in self.allpar.itertuples():  # cycle through all parameters
            if row.type == 'd' or row.type == 'fd' or row.type == 'i' or row.type == 'fi':
                priorentropy += np.log(row.u_fit - row.l_fit) / np.log(2)
                # calculate prior entropy for parameters to be marginalized (dependent parameters)
                if row.type == 'd' or row.type == 'fd':
                    priorentropy_marginal += np.log(row.u_fit - row.l_fit) / np.log(2)
        return priorentropy, priorentropy_marginal

    def gpcam_instrument(self, data, Test=False):
        print("This is the current length of the data received by gpCAM: ", len(data))
        print("Suggested by gpCAM: ", data)
        for entry in data:
            if Test:
                # value = np.sin(np.linalg.norm(entry["position"]))
                # value = np.array(entry['position']).sum() / 1000
                value = (entry['position'][0] - entry['position'][1]) ** 2
                time0 = entry['position'][2]
                time1 = entry['position'][3]
                time2 = entry['position'][4]
                tf = 14400 / (time0 + time1 + time2)
                value += np.log10(time0 * tf) * 0.5
                value += np.log10(time1 * tf) * 1.5
                value += np.log10(time2 * tf) * 1
                entry['value'] = value
                print('Value: ', entry['value'])
                variance = None  # 0.01 * np.abs(entry['value'])
                entry['variance'] = variance
            else:
                if self.fitter == 'LM':
                    me = []
                    for _ in range(self.lm_iterations):
                        marginal_entropy = self.work_on_iteration(position=entry['position'],
                                                                  gpiteration=self.gpiteration)
                        me.append(self.priorentropy_marginal - marginal_entropy)
                    value = numpy.mean(me)
                    variance = numpy.var(me)
                else:
                    marginal_entropy = self.work_on_iteration(position=entry['position'],
                                                              gpiteration=self.gpiteration)
                    value = self.priorentropy_marginal - marginal_entropy
                    variance = None
                entry["value"] = value
                entry['variance'] = variance
                # entry["cost"]  = [np.array([0,0]),entry["position"],np.sum(entry["position"])]

            self.gpCAMstream['position'].append(entry['position'])
            self.gpCAMstream['value'].append(value)
            self.gpCAMstream['variance'].append(variance)
            self.save_results_gpcam(self.spath)
            self.gpiteration += 1
        return data

    def gpcam_prediction(self, my_ae):
        # create a flattened array of all positions to be evaluated, maximize the use of numpy
        prediction_positions = np.array(self.axes[0])
        for i in range(1, len(self.axes)):
            a = np.array([prediction_positions] * len(self.axes[i]))
            # transpose the first two axes of a only
            newshape = np.linspace(0, len(a.shape) - 1, len(a.shape), dtype=int)
            newshape[0] = 1
            newshape[1] = 0
            a = np.transpose(a, newshape)
            b = np.array([self.axes[i]] * prediction_positions.shape[0])
            prediction_positions = np.dstack((a, b))
            # now flatten the first two dimensions
            newshape = list(prediction_positions.shape[1:])
            newshape[0] = newshape[0] * prediction_positions.shape[0]
            newshape = tuple(newshape)
            prediction_positions = np.reshape(prediction_positions, newshape)

        res = my_ae.gp_optimizer.posterior_mean(prediction_positions)
        f = res["f(x)"]
        self.prediction_gpcam = f.reshape(self.steplist)

        path1 = path.join(self.spath, 'plots')
        if not path.isdir(path1):
            mkdir(path1)
        # self.plot_arr(self.prediction_gpcam, filename=path.join(path1, 'prediction_gpcam'), mark_maximum=True)

        if self.show_support_points:
            support_points = np.array(self.gpCAMstream['position'])
        else:
            support_points = None

        self.plot_arr(self.prediction_gpcam,
                      filename=path.join(path1, 'prediction_gpcam'), mark_maximum=True,
                      support_points=support_points)

    def gridsearch_iterate_over_all_indices(self, refinement=False):
        bWorkedOnIndex = False
        # the complicated iteration syntax is due the unknown dimensionality of the results space / arrays
        it = np.nditer(self.results_gmm, flags=['multi_index'])
        while not it.finished:
            itindex = it.multi_index
            # parameter grids can be symmetric, and only one of the symmetry-related indices
            # will be calculated unless calc_symmetric is True
            if all(itindex[i] <= itindex[i + 1] for i in range(len(itindex) - 1)) or self.calc_symmetric:
                # run MCMC if it is first time or the value in results is inf
                invalid_result = np.isinf(self.results_gmm[itindex]) or fabs(self.results_gmm[itindex]) > 10000
                insufficient_iterations = self.n_mvn[itindex] < self.miniter

                outlier = False
                if refinement:
                    # during refinement, check whether the GMM value follows that of the MVN with respect to its
                    # nearest neighbors, implemented convolution because of ill-defined origin of scipy convolute
                    conv_MVN = convolute(self.results_mvn)
                    conv_GMM = convolute(self.results_gmm)
                    dMVN = conv_MVN[itindex] - self.results_mvn[itindex]
                    dGMM = conv_GMM[itindex] - self.results_gmm[itindex]
                    if fabs(dMVN - dGMM) > self.convergence:
                        outlier = True

                # Do we need to work on this particular index?
                if outlier or invalid_result or insufficient_iterations or self.bFetchMode:
                    bWorkedOnIndex = True
                    _ = self.work_on_iteration(it=it)

            it.iternext()
        return bWorkedOnIndex

    def gridsearch_writeout_result(self, _itindex, avg_mvn, avg_gmm, avg_mvn_marginal, avg_gmm_marginal,
                                   points_median, points_std, parnames):
        # writes out entropy and parameter results into numpy arrays
        if not self.calc_symmetric:
            # since symmetry-related points in the optimization were not calculated twice, the current
            # result is copied to all permutations without repetition of the parameter index
            # this is convenient for all interchangeable parameters, such as multiple solvent contrasts
            indexlist = []
            permutated = itertools.permutations(_itindex)
            for element in permutated:
                if element not in indexlist:
                    indexlist.append(element)
        else:
            # otherwise copy out to single parameter
            indexlist = [_itindex]

        for index in indexlist:
            self.n_gmm[index] += 1.0
            self.n_mvn[index] += 1.0
            self.n_gmm_marginal[index] += 1.0
            self.n_mvn_marginal[index] += 1.0
            n = self.n_gmm[index]

            old_mvn = self.results_mvn[index]
            old_gmm = self.results_gmm[index]
            old_mvn_marginal = self.results_mvn_marginal[index]
            old_gmm_marginal = self.results_gmm_marginal[index]

            self.results_mvn[index] = running_mean(self.results_mvn[index], n, avg_mvn)
            self.results_gmm[index] = running_mean(self.results_gmm[index], n, avg_gmm)
            self.results_mvn_marginal[index] = running_mean(self.results_mvn_marginal[index], n, avg_mvn_marginal)
            self.results_gmm_marginal[index] = running_mean(self.results_gmm_marginal[index], n, avg_gmm_marginal)

            for i in range(self.par_median.shape[0]):
                # parameter indices from fitting results and entropy module might be different
                j = parnames.index(self.parlist[i])
                self.par_median[(i,) + index] = running_mean(self.par_median[(i,) + index], n, points_median[j])
                # for par std the average is calculated, not a sqstd of par_median
                self.par_std[(i,) + index] = running_mean(self.par_std[(i,) + index], n, points_std[j])

            self.sqstd_mvn[index] = running_sqstd(self.sqstd_mvn[index], n, avg_mvn, old_mvn,
                                                  self.results_mvn[index])
            self.sqstd_gmm[index] = running_sqstd(self.sqstd_gmm[index], n, avg_gmm, old_gmm,
                                                  self.results_gmm[index])
            self.sqstd_mvn_marginal[index] = running_sqstd(self.sqstd_mvn_marginal[index], n, avg_mvn_marginal,
                                                           old_mvn_marginal, self.results_mvn_marginal[index])
            self.sqstd_gmm_marginal[index] = running_sqstd(self.sqstd_gmm_marginal[index], n, avg_gmm_marginal,
                                                           old_gmm_marginal, self.results_gmm_marginal[index])

    def load_results(self, dirname):
        path1 = path.join(dirname, 'results')
        if self.fitter != 'LM':
            self.results_gmm = np.load(path.join(path1, 'GMM_entropy.npy'))
            self.results_gmm_marginal = np.load(path.join(path1, 'GMM_entropy_marginal.npy'))
            self.n_gmm = np.load(path.join(path1, 'GMM_n.npy'))
            self.n_gmm_marginal = np.load(path.join(path1, 'GMM_n_marginal.npy'))
            self.sqstd_gmm = np.load(path.join(path1, 'GMM_sqstd.npy'))
            self.sqstd_gmm_marginal = np.load(path.join(path1, 'GMM_sqstd_marginal.npy'))

        if path.isfile(path.join(path1, 'Prediction_gpcam.npy')):
            self.prediction_gpcam = np.load(path.join(path1, 'Prediction_gpcam.npy'))

        self.results_mvn = np.load(path.join(path1, 'MVN_entropy.npy'))
        self.results_mvn_marginal = np.load(path.join(path1, 'MVN_entropy_marginal.npy'))
        self.n_mvn = np.load(path.join(path1, 'MVN_n.npy'))
        self.n_mvn_marginal = np.load(path.join(path1, 'MVN_n_marginal.npy'))
        self.sqstd_mvn = np.load(path.join(path1, 'MVN_sqstd.npy'))
        self.sqstd_mvn_marginal = np.load(path.join(path1, 'MVN_sqstd_marginal.npy'))
        self.par_median = np.load(path.join(path1, 'par_median.npy'))
        self.par_std = np.load(path.join(path1, 'par_std.npy'))

    def run_fit(self, molstat_instance, iteration, dirname, fulldirname):
        # wait for a job to finish before submitting next cluster job
        if self.bClusterMode:
            self.waitforjob()

        # run MCMC either cluster or local
        if self.bClusterMode:
            # write runscript
            mcmc_iteration = str(iteration)
            mcmc_dir = dirname
            # replaces the placeholders in slurmscript with variables above
            script = self.slurmscript.format(**locals())

            file = open(path.join(fulldirname, 'runscript'), 'w')
            file.writelines(script)
            file.close()

            lCommand = ['sbatch', path.join(fulldirname, 'runscript')]
            Popen(lCommand)
            self.joblist.append(iteration)

        else:
            if self.fitter == 'MCMC':
                molstat_instance.Interactor.fnRunMCMC(burn=self.mcmcburn, steps=self.mcmcsteps, batch=True)
            else:
                # The best-fit is still loaded in problem from during data simulation. Therefore, we do not
                # reload but use this starting point for a LM, hopefully gaining a speed-up.
                molstat_instance.Interactor.fnRunMCMC(fitter='LM', batch=True, reload_problem=False)
        return

    def plot_results(self, mark_maximum=False):

        path1 = path.join(self.spath, 'plots')
        if not path.isdir(path1):
            mkdir(path1)

        if self.fitter != 'LM':
            self.plot_arr(self.results_gmm, arr_variance=np.sqrt(self.sqstd_gmm), vallabel='Entropy [bits]',
                          filename=path.join(path1, 'GMM_entropy'))
            self.plot_arr(self.results_gmm_marginal, arr_variance=np.sqrt(self.sqstd_gmm_marginal),
                          vallabel='Entropy [bits]', filename=path.join(path1, 'GMM_entropy_marginal'))
            self.plot_arr(self.priorentropy - self.results_gmm, arr_variance=np.sqrt(self.sqstd_gmm),
                          vallabel='information gain [bits]', filename=path.join(path1, 'GMM_infocontent'), valmin=0,
                          mark_maximum=mark_maximum)
            self.plot_arr(self.priorentropy - self.results_gmm, arr_variance=np.sqrt(self.sqstd_gmm),
                          vallabel='information gain [bits]', filename=path.join(path1, 'GMM_infocontent'), valmin=0,
                          mark_maximum=mark_maximum)
            self.plot_arr(self.priorentropy_marginal - self.results_gmm_marginal,
                          arr_variance=np.sqrt(self.sqstd_gmm_marginal), vallabel='information gain [bits]',
                          filename=path.join(path1, 'GMM_infocontent_marginal'), valmin=0,
                          valmax=self.upper_info_plotlevel,
                          mark_maximum=mark_maximum)
            self.plot_arr(self.n_gmm, vallabel='computations', filename=path.join(path1, 'GMM_n'), valmin=0)
            self.plot_arr(self.n_gmm_marginal, vallabel='computations', filename=path.join(path1, 'GMM_n_marginal'),
                          valmin=0)

        self.plot_arr(self.results_mvn, arr_variance=np.sqrt(self.sqstd_mvn), vallabel='Entropy [bits]',
                      filename=path.join(path1, 'MVN_entropy'))
        self.plot_arr(self.results_mvn_marginal, arr_variance=np.sqrt(self.sqstd_mvn_marginal),
                      vallabel='Entropy [bits]', filename=path.join(path1, 'MVN_entropy_marginal'))
        self.plot_arr(self.priorentropy - self.results_mvn, arr_variance=np.sqrt(self.sqstd_mvn),
                      vallabel='information gain [bits]', filename=path.join(path1, 'MVN_infocontent'), valmin=0,
                      mark_maximum=mark_maximum)
        self.plot_arr(self.priorentropy_marginal - self.results_mvn_marginal,
                      arr_variance=np.sqrt(self.sqstd_mvn_marginal), vallabel='information gain [bits]',
                      filename=path.join(path1, 'MVN_infocontent_marginal'), valmin=0, valmax=self.upper_info_plotlevel,
                      mark_maximum=mark_maximum)
        self.plot_arr(self.n_mvn,  vallabel='computations', filename=path.join(path1, 'MVN_n'), valmin=0)
        self.plot_arr(self.n_mvn_marginal,  vallabel='computations', filename=path.join(path1, 'MVN_n_marginal'),
                      valmin=0)

        for i, parname in enumerate(self.parlist):
            self.plot_arr(self.par_median[i], arr_variance=self.par_std[i], vallabel=parname,
                          filename=path.join(path1, 'Par_' + parname + '_median'))
            self.plot_arr(self.par_std[i], arr_variance=None, vallabel=parname,
                          filename=path.join(path1, 'Par_' + parname + '_std'))

    def plot_arr(self, arr_value, arr_variance=None, filename='plot', mark_maximum=False, valmin=None, valmax=None,
                 levels=20, niceticks=False, vallabel='z', support_points=None):
        # onecolormaps = [plt.cm.Greys, plt.cm.Purples, plt.cm.Blues, plt.cm.Greens, plt.cm.Oranges, plt.cm.Reds]
        ec = plt.cm.coolwarm

        path1 = path.join(self.spath, 'plots')

        if len(arr_value.shape) == 1:
            ax0 = self.axes[0]
            sp0 = self.steppar['unique_name'].tolist()[0]
            if arr_variance is not None:
                dy = np.sqrt(arr_variance)
            else:
                dy = None
            save_plot_1d(ax0, arr_value, dy=dy, xlabel=sp0, ylabel=vallabel, filename=path.join(path1, filename),
                         ymin=valmin, ymax=valmax, levels=levels, niceticks=niceticks, keep_plots=self.keep_plots)

        elif len(arr_value.shape) == 2:
            # numpy array and plot axes are reversed
            ax1 = self.axes[0]
            ax0 = self.axes[1]
            sp1 = self.steppar['unique_name'].tolist()[0]
            sp0 = self.steppar['unique_name'].tolist()[1]
            save_plot_2d(ax0, ax1, arr_value, xlabel=sp0, ylabel=sp1, color=ec,
                         filename=path.join(path1, filename), zmin=valmin, zmax=valmax, levels=levels,
                         mark_maximum=mark_maximum, keep_plots=self.keep_plots, support_points=support_points)

        elif len(arr_value.shape) == 3 and arr_value.shape[0] < 6:
            ax2 = self.axes[1]
            ax1 = self.axes[2]
            sp2 = self.steppar['unique_name'].tolist()[1]
            sp1 = self.steppar['unique_name'].tolist()[2]
            for slice_n in range(arr_value.shape[0]):
                save_plot_2d(ax1, ax2, arr_value[slice_n], xlabel=sp1, ylabel=sp2, color=ec,
                             filename=path.join(path1, filename+'_'+str(slice_n)), zmin=valmin, zmax=valmax,
                             levels=levels, mark_maximum=mark_maximum, keep_plots=self.keep_plots)

        if len(arr_value.shape) >= 3:
            # plot projections onto two parameters at a time
            for i in range(len(self.steppar)):
                for j in range(i):
                    ax2 = self.axes[i]
                    ax1 = self.axes[j]
                    sp2 = self.steppar['unique_name'].tolist()[i]
                    sp1 = self.steppar['unique_name'].tolist()[j]
                    projection = np.empty((self.steplist[i], self.steplist[j]))
                    for k in range(self.steplist[i]):
                        for ll in range(self.steplist[j]):
                            projection[k, ll] = np.take(np.take(arr_value, indices=k, axis=i), indices=ll, axis=j).max()
                    save_plot_2d(ax1, ax2, projection, xlabel=sp1, ylabel=sp2, color=ec,
                                 filename=path.join(path1, filename+'_'+sp1+'_'+sp2), zmin=valmin, zmax=valmax,
                                 levels=levels, mark_maximum=mark_maximum, keep_plots=self.keep_plots)

    def run_optimization(self):
        if self.optimizer == 'grid':
            self.run_optimization_grid()
        elif self.optimizer == 'gpCAM' or self.optimizer == 'gpcam':
            self.run_optimization_gpcam()

    def run_optimization_gpcam(self):
        # Using the gpCAM global optimizer, follows the example from the gpCAM website

        # initialization
        # feel free to try different acquisition functions, e.g. optional_acq_func, "covariance", "shannon_ig"
        # note how costs are defined in for the autonomous experimenter
        parlimits = []
        for row in self.steppar.iterrows():
            parlimits.append([row[1].l_sim, row[1].u_sim])
        parlimits = np.array(parlimits)
        numpars = len(parlimits)

        x = self.gpCAMstream['position']
        y = self.gpCAMstream['value']
        v = self.gpCAMstream['variance']

        if len(x) >= 1:
            # use any previously computed results
            x = np.array(x)
            y = np.array(y)
            v = np.array(v)
            self.gpiteration = len(x)
            bFirstEval = False
        else:
            x = None
            y = None
            v = None
            self.gpiteration = 0
            bFirstEval = True

        hyperpars = np.ones([numpars + 1])
        # the zeroth hyper bound is associated with a signal variance
        # the others with the length scales of the parameter inputs
        hyper_bounds = np.array([[0.001, 100]] * (numpars + 1))
        for i in range(len(parlimits)):
            delta = parlimits[i][1] - parlimits[i][0]
            hyper_bounds[i + 1] = [delta * 1e-3, delta * 1e1]

        self.my_ae = AutonomousExperimenterGP(parlimits, hyperpars, hyper_bounds,
                                              init_dataset_size=self.gpcam_init_dataset_size,
                                              instrument_func=self.gpcam_instrument,
                                              acq_func=self.acq_func,  # optional_acq_func,
                                              # cost_func = optional_cost_function,
                                              # cost_update_func = optional_cost_update_function,
                                              x=x, y=y, v=v,
                                              # cost_func_params={"offset": 5.0, "slope": 10.0},
                                              kernel_func=None, use_inv=True,
                                              communicate_full_dataset=False, ram_economy=True)

        # save and evaluate initial data set if it has been freshly calculate
        if bFirstEval:
            self.save_results_gpcam(self.spath)
            self.gpcam_prediction(self.my_ae)

        while len(self.my_ae.x) < self.gpcam_iterations:
            print("length of the dataset: ", len(self.my_ae.x))
            self.my_ae.train(method="global", max_iter=10000)  # or not, or both, choose "global","local" and "hgdl"
            # update hyperparameters in case they are optimized asynchronously
            self.my_ae.train(method="local")  # or not, or both, choose between "global","local" and "hgdl"
            # training and client can be killed if desired and in case they are optimized asynchronously
            # self.my_ae.kill_training()
            if self.gpcam_step is not None:
                target_iterations = len(self.my_ae.x) + self.gpcam_step
                retrain_async_at = []
            else:
                target_iterations = self.gpcam_iterations
                retrain_async_at = np.logspace(start=np.log10(len(self.my_ae.x)),
                                               stop=np.log10(self.gpcam_iterations / 2), num=3, dtype=int)
            # run the autonomous loop
            self.my_ae.go(N=target_iterations,
                          retrain_async_at=retrain_async_at,
                          retrain_globally_at=[],
                          retrain_locally_at=[],
                          acq_func_opt_setting=lambda number: "global" if number % 2 == 0 else "local",
                          training_opt_max_iter=20,
                          training_opt_pop_size=10,
                          training_opt_tol=1e-6,
                          acq_func_opt_max_iter=20,
                          acq_func_opt_pop_size=20,
                          acq_func_opt_tol=1e-6,
                          number_of_suggested_measurements=1,
                          acq_func_opt_tol_adjust=0.1)

            # training and client can be killed if desired and in case they are optimized asynchronously
            if self.gpcam_step is None:
                self.my_ae.kill_training()
            self.save_results_gpcam(self.spath)
            self.gpcam_prediction(self.my_ae)

    def run_optimization_grid(self):
        # Grid search
        # every index has at least one result before re-analyzing any data point (refinement)
        bRefinement = False
        while True:
            bWorkedOnAnyIndex = self.gridsearch_iterate_over_all_indices(bRefinement)
            if not bWorkedOnAnyIndex:
                if not bRefinement:
                    # all indices have the minimum number of iterations -> start refinement
                    bRefinement = True
                else:
                    # done with refinement
                    break

            if self.bClusterMode or self.bFetchMode:
                # never repeat iterations on cluster or when just calculating entropies
                break

        # wait for all jobs to finish
        if self.bClusterMode:
            while self.joblist:
                self.waitforjob(bFinish=True)

    def save_results_gpcam(self, dirname):
        path1 = path.join(dirname, 'results')
        if not path.isdir(path1):
            mkdir(path1)
        with open(path.join(self.spath, 'results', 'gpCAMstream.pkl'), 'wb') as file:
            pickle.dump(self.gpCAMstream, file)

    def save_results(self, dirname):
        path1 = path.join(dirname, 'results')
        if not path.isdir(path1):
            mkdir(path1)

        if self.fitter != 'LM':
            np.save(path.join(path1, 'GMM_entropy'), self.results_gmm, allow_pickle=False)
            np.save(path.join(path1, 'GMM_entropy_marginal'), self.results_gmm_marginal, allow_pickle=False)
            np.save(path.join(path1, 'GMM_infocontent'), self.priorentropy - self.results_gmm, allow_pickle=False)
            np.save(path.join(path1, 'GMM_infocontent_marginal'),
                    self.priorentropy_marginal - self.results_gmm_marginal, allow_pickle=False)
            np.save(path.join(path1, 'GMM_sqstd'), self.sqstd_gmm, allow_pickle=False)
            np.save(path.join(path1, 'GMM_sqstd_marginal'), self.sqstd_gmm_marginal, allow_pickle=False)
            np.save(path.join(path1, 'GMM_n'), self.n_gmm, allow_pickle=False)
            np.save(path.join(path1, 'GMM_n_marginal'), self.n_gmm_marginal, allow_pickle=False)

        np.save(path.join(path1, 'MVN_entropy'), self.results_mvn, allow_pickle=False)
        np.save(path.join(path1, 'MVN_entropy_marginal'), self.results_mvn_marginal, allow_pickle=False)
        np.save(path.join(path1, 'MVN_infocontent'), self.priorentropy - self.results_mvn, allow_pickle=False)
        np.save(path.join(path1, 'MVN_infocontent_marginal'), self.priorentropy_marginal - self.results_mvn_marginal,
                allow_pickle=False)
        np.save(path.join(path1, 'MVN_sqstd'), self.sqstd_mvn, allow_pickle=False)
        np.save(path.join(path1, 'MVN_sqstd_marginal'), self.sqstd_mvn_marginal, allow_pickle=False)
        np.save(path.join(path1, 'MVN_n'), self.n_mvn, allow_pickle=False)
        np.save(path.join(path1, 'MVN_n_marginal'), self.n_mvn_marginal, allow_pickle=False)
        np.save(path.join(path1, 'par_median'), self.par_median, allow_pickle=False)
        np.save(path.join(path1, 'par_std'), self.par_std, allow_pickle=False)

        # save to txt when not more than two-dimensional array
        if len(self.steplist) <= 2:
            if self.fitter != 'LM':
                np.savetxt(path.join(path1, 'GMM_entropy.txt'), self.results_gmm - 0)
                np.savetxt(path.join(path1, 'GMM_entropy_marginal.txt'), self.results_gmm_marginal - 0)
                np.savetxt(path.join(path1, 'GMM_infocontent.txt'), self.priorentropy - self.results_gmm)
                np.savetxt(path.join(path1, 'GMM_infocontent_marginal.txt'), self.priorentropy_marginal -
                           self.results_gmm_marginal)
                np.savetxt(path.join(path1, 'GMM_sqstd.txt'), self.sqstd_gmm - 0)
                np.savetxt(path.join(path1, 'GMM_sqstd_marginal.txt'), self.sqstd_gmm_marginal - 0)
                np.savetxt(path.join(path1, 'GMM_n.txt'), self.n_gmm - 0)
                np.savetxt(path.join(path1, 'GMM_n_marginal.txt'), self.n_gmm_marginal - 0)

            np.savetxt(path.join(path1, 'MVN_entropy.txt'), self.results_mvn - 0)
            np.savetxt(path.join(path1, 'MVN_entropy_marginal.txt'), self.results_mvn_marginal - 0)
            np.savetxt(path.join(path1, 'MVN_infocontent.txt'), self.priorentropy - self.results_mvn)
            np.savetxt(path.join(path1, 'MVN_infocontent_marginal.txt'), self.priorentropy_marginal -
                       self.results_mvn_marginal)
            np.savetxt(path.join(path1, 'MVN_sqstd.txt'), self.sqstd_mvn - 0)
            np.savetxt(path.join(path1, 'MVN_sqstd_marginal.txt'), self.sqstd_mvn_marginal - 0)
            np.savetxt(path.join(path1, 'MVN_n.txt'), self.n_mvn - 0)
            np.savetxt(path.join(path1, 'MVN_n_marginal.txt'), self.n_mvn_marginal - 0)
            i = 0
            for parname in self.parlist:
                np.savetxt(path.join(path1, 'Par_' + parname + '_median.txt'), self.par_median[i])
                np.savetxt(path.join(path1, 'Par_' + parname + '_std.txt'), self.par_std[i])
                i += 1

        # save three-dimensional array in slices of the first parameter
        if len(self.steplist) == 3 and self.results_gmm.shape[0] < 6:
            for sl in range(self.results_gmm.shape[0]):
                if self.fitter != 'LM':
                    np.savetxt(path.join(path1, 'GMM_entropy_' + str(sl) + '.txt'), self.results_gmm[sl])
                    np.savetxt(path.join(path1, 'GMM_entropy_marginal_' + str(sl) + '.txt'),
                               self.results_gmm_marginal[sl])
                    np.savetxt(path.join(path1, 'GMM_infocontent_' + str(sl) + '.txt'), self.priorentropy -
                               self.results_gmm[sl])
                    np.savetxt(path.join(path1, 'GMM_infocontent_marginal_' + str(sl) + '.txt'),
                               self.priorentropy_marginal - self.results_gmm_marginal[sl])
                    np.savetxt(path.join(path1, 'GMM_sqstd_' + str(sl) + '.txt'), self.sqstd_gmm[sl])
                    np.savetxt(path.join(path1, 'GMM_sqstd_marginal_' + str(sl) + '.txt'),
                               self.sqstd_gmm_marginal[sl])
                    np.savetxt(path.join(path1, 'GMM_n_' + str(sl) + '.txt'), self.n_gmm[sl])
                    np.savetxt(path.join(path1, 'GMM_n_marginal_' + str(sl) + '.txt'), self.n_gmm_marginal[sl])

                np.savetxt(path.join(path1, 'MVN_entropy_' + str(sl) + '.txt'), self.results_mvn[sl])
                np.savetxt(path.join(path1, 'MVN_entropy_marginal_' + str(sl) + '.txt'),
                           self.results_mvn_marginal[slice])
                np.savetxt(path.join(path1, 'MVN_infocontent_' + str(sl) + '.txt'), self.priorentropy -
                           self.results_mvn[sl])
                np.savetxt(path.join(path1, 'MVN_infocontent_marginal_' + str(sl) + '.txt'),
                           self.priorentropy_marginal - self.results_mvn_marginal[sl])
                np.savetxt(path.join(path1, 'MVN_sqstd_' + str(sl) + '.txt'), self.sqstd_mvn[sl])
                np.savetxt(path.join(path1, 'MVN_sqstd_marginal_' + str(sl) + '.txt'),
                           self.sqstd_mvn_marginal[sl])
                np.savetxt(path.join(path1, 'MVN_n_' + str(sl) + '.txt'), self.n_mvn[sl])
                np.savetxt(path.join(path1, 'MVN_n_marginal_' + str(sl) + '.txt'), self.n_mvn_marginal[sl])

    def set_sim_pars_for_iteration(self, it=None, position=None):
        def _str2int(st):
            if st == '*':
                i = 0
            else:
                i = int(st)
            return i

        def _fill_config(configurations, parname, parvalue, dataset, configuration):
            if dataset == '*':
                for ds in range(len(configurations)):
                    if configuration == '*':
                        for cf in range(len(configurations[ds])):
                            configurations[ds][cf][parname] = parvalue
                    else:
                        cf = _str2int(configuration)
                        configurations[ds][cf][parname] = parvalue
            else:
                ds = _str2int(dataset)
                if configuration == '*':
                    for cf in range(len(configurations[ds])):
                        configurations[ds][cf][parname] = parvalue
                else:
                    cf = _str2int(configuration)
                    configurations[ds][cf][parname] = parvalue
            return configurations

        def _set_background(configurations, dset, config, value):
            # calculate background value
            cb = 0
            if self.mode == 'SANS_linear':
                cb = self.background_rule['y_intercept'] + self.background_rule['slope'] * value
            if self.mode == 'water':
                # hardcoded water cb = 0.07 for D2O, cb = 1.00 for H2O
                cb = 0.9245 - 0.1348 * value

            # search if there is a designated background parameter that takes the cb instead of a configuration
            bFoundBackground = False
            for row in self.allpar.itertuples():
                if row.dataset == 'b' + str(dset) or row.dataset == 'b' or row.dataset == 'b*':
                    self.simpar.loc[self.simpar['par'] == row.par, 'value'] = cb
                    bFoundBackground = True
                    if cb < row.l_fit or cb > row.u_fit:
                        print('Parameter with background rule outside of fit boundaries!')
                        print('Parameter: ', row.par)
                        print('Value: ', cb)
                        print('Lower boundary: ', row.l_fit)
                        print('Upper boundary: ', row.u_fit)
                        sys.exit(1)
                    return configurations
            if bFoundBackground:
                return configurations

            # change buffer crosssection in configurations
            configurations = _fill_config(configurations, 'differential_cross_section_buffer', cb, dset, config)

            return configurations

        # Configurations are imported externally, if not empty configuration initialization here
        # In any case, missing parameters are set to the default in the API simulation routines
        configurations = self.configuration
        if configurations is None:
            configurations = [[{}]]

        # cycle through all parameters
        isim = 0
        for row in self.allpar.itertuples():
            simvalue = None
            # is it a parameter to iterate over?
            if row.unique_name in self.steppar['unique_name'].tolist():
                lsim = self.steppar.loc[self.steppar['unique_name'] == row.unique_name, 'l_sim'].iloc[0]
                stepsim = self.steppar.loc[self.steppar['unique_name'] == row.unique_name, 'step_sim'].iloc[0]
                value = self.steppar.loc[self.steppar['unique_name'] == row.unique_name, 'value'].iloc[0]
                lfit = self.steppar.loc[self.steppar['unique_name'] == row.unique_name, 'l_fit'].iloc[0]
                ufit = self.steppar.loc[self.steppar['unique_name'] == row.unique_name, 'u_fit'].iloc[0]

                if it is not None:
                    # grid mode, calculate value from evaluation grid and index
                    simvalue = lsim + stepsim * it.multi_index[isim]
                    # print(self.steppar['unique_name'], '  ', simvalue, '\n')
                else:
                    # gpcam mode, use position suggested by gpcam
                    simvalue = position[isim]

                if row.type == 'd' or row.type == 'fd' or row.type == 'i' or row.type == 'fi':
                    if row.type == 'fd' or row.type == 'fi':
                        # fixed fit boundaries, not floating, for such things as volfracs between 0 and 1
                        lowersim = lfit
                        uppersim = ufit
                    else:
                        lowersim = simvalue - (value - lfit)
                        uppersim = simvalue + (ufit - value)

                    if 'b' not in row.dataset:
                        # only change simpar when it will not be filled by a background rule
                        self.simpar.loc[self.simpar['unique_name'] == row.unique_name, 'value'] = simvalue
                        self.molstat.Interactor.fnReplaceParameterLimitsInSetup(row.par, lowersim, uppersim)
                    else:
                        print('Parameter with background rule cannot be varied!')
                        print('Parameter: ', row.unique_name)
                        sys.exit(1)
                else:
                    # must be instrument parameter
                    configurations = _fill_config(configurations, row.par, simvalue, row.dataset, row.configuration)
                    for indx in self.simpar.index:
                        if self.simpar['unique_name'][indx] == row.unique_name and \
                                self.simpar['dataset'][indx] == row.dataset and \
                                self.simpar['configuration'][indx] == row.configuration:
                            self.simpar.iloc[indx, self.simpar.columns.get_loc('value')] = simvalue
                            break

                isim += 1
            elif row.type == 'd' or row.type == 'fd' or row.type == 'i' or row.type == 'fi':
                # We also fill parameters with background rule here. They are supposed to have fit boundaries
                # that encompass all possible values, and the particular row.par does not matter.
                self.molstat.Interactor.fnReplaceParameterLimitsInSetup(row.par, row.l_fit, row.u_fit)
            else:
                # must be instrument parameter
                configurations = _fill_config(configurations, row.par, row.value, row.dataset, row.configuration)

            if (row.type != 'n') and (row.dataset != '_') and ('b' not in row.dataset):
                # this is a parameter that will determine one or more isotropic backgrounds
                if simvalue is None:
                    configurations = _set_background(configurations, row.dataset, row.configuration, row.value)
                else:
                    configurations = _set_background(configurations, row.dataset, row.configuration, simvalue)

        simparsave = self.simpar.loc[:, ['par', 'value']]
        simparsave.to_csv(path.join(self.spath, 'simpar.dat'), sep=' ', header=None, index=False)
        return configurations

    def waitforjob(self, bFinish=False):
        # finish flag means that parent is waiting for all jobs to finish and not because of a too long
        # job queue
        if len(self.joblist) >= self.jobmax or bFinish:
            repeat = True
            while repeat:
                for job in self.joblist:
                    if path.isfile(path.join(self.spath, 'iteration_' + str(job), 'save', self.runfile + '-chain.mc')):
                        # wait 2 minutes to allow all output files to be written
                        sleep(180)
                        # zip up finished job
                        # shutil.make_archive('iteration_' + str(job), 'zip', 'iteration_' + str(job))
                        # call(['rm', '-r', 'iteration_' + str(job)])
                        self.joblist.remove(job)
                        repeat = False
                        break
                sleep(60)
        return

    def work_on_iteration(self, it=None, position=None, gpiteration=None):

        if it is not None:
            # grid mode, calculate which iteration we are on from itindex
            # Recover itindex from 'it'. This is not an argument to the function anymore, as work_on_index might
            # be used independently of iterate_over_all_indices
            itindex = it.multi_index
            if self.calc_symmetric:
                iteration = it.iterindex
            else:
                # TODO: This is potentially expensive. Find better method for symmetry-conscious calculation
                it2 = np.nditer(self.results_gmm, flags=['multi_index'])
                iteration = 0
                while not it2.finished:
                    itindex2 = it2.multi_index
                    if all(itindex2[i] <= itindex2[i + 1] for i in range(len(itindex2) - 1)):
                        # iterations are only increased if this index is not dropped because of symmetry
                        iteration += 1
                    it2.iternext()
        else:
            # gpcam mode
            iteration = gpiteration
            itindex = None

        dirname = 'iteration_' + str(iteration)
        fulldirname = path.join(self.spath, dirname)
        path1 = path.join(fulldirname, 'save')
        chainname = path.join(path1, self.runfile+'-chain.mc')

        # most relevant result for a particular index to return for general use of this function
        avg_gmm_marginal = 0

        # fetch mode and cluster mode are exclusive
        if not self.bFetchMode:
            # LM sometimes produces a singular matrix, which we try to avoid
            fit_counter = 0
            fit_success = False
            while not fit_success:
                # run a new fit, preparations are done in the root directory and the new fit is copied into the
                # iteration directory, preparations in the iterations directory are not possible, because it would
                # be lacking a result directory, which is needed for restoring a state/parameters
                self.molstat.Interactor.fnBackup(target=path.join(self.spath, 'simbackup'))
                configurations = self.set_sim_pars_for_iteration(it, position)
                self.molstat.fnSimulateData(mode=self.mode, liConfigurations=configurations, qmin=self.qmin,
                                            qmax=self.qmax, qrangefromfile=self.qrangefromfile, t_total=self.t_total)
                self.molstat.Interactor.fnBackup(origin=self.spath, target=fulldirname)
                # previous save needs to be removed as output serves as flag for HPC job termination
                if path.isdir(path1):
                    shutil.rmtree(path1)
                self.molstat.Interactor.fnRemoveBackup(target=path.join(self.spath, 'simbackup'))

                # changing the working directory became necessary at some point for loading the correct data
                os.chdir(fulldirname)
                molstat_iter = molstat.CMolStat(fitsource=self.fitsource, spath=fulldirname, mcmcpath='save',
                                                runfile=self.runfile, load_state=False)

                if self.fitter == 'LM':
                    # copy best-fit parameters from data simulation instance
                    molstat_iter.Interactor.problem.setp(self.molstat.Interactor.problem.getp())
                    self.run_fit(molstat_iter, iteration, dirname, fulldirname)
                    # use covariance matrix for entropy calculation in case of LM
                    avg_gmm_marginal = self.calc_entropy_for_iteration(molstat_iter, itindex=itindex, cov=True)
                    fit_counter += 1
                    if avg_gmm_marginal is not None:
                        fit_success = True
                    else:
                        if fit_counter > 5:
                            print("Singular matrix encountered for five times in LM.")
                            print("Assume information gain of zero.")
                            avg_gmm_marginal = 0
                            fit_success = True
                else:
                    self.run_fit(molstat_iter, iteration, dirname, fulldirname)
                    fit_counter += 1
                    fit_success = True

                os.chdir(self.spath)


        # Do not run entropy calculation when on cluster, or no valid result, or entropy from covariance via LM.
        if self.fitter != 'LM':
            bPriorResultExists = path.isfile(chainname) or path.isfile(chainname + '.gz')
            if not self.bClusterMode and bPriorResultExists:
                molstat_iter = molstat.CMolStat(fitsource=self.fitsource, spath=fulldirname, mcmcpath='save',
                                                runfile=self.runfile)
                avg_gmm_marginal = self.calc_entropy_for_iteration(molstat_iter, itindex=itindex)

        # delete big files except in Cluster mode. They are needed there for future fetching
        if self.remove_fit_dir and not self.bClusterMode:
            shutil.rmtree(fulldirname)
        elif self.deldir and not self.bClusterMode:
            rm_file(path.join(path1, self.runfile+'-point.mc'))
            rm_file(path.join(path1, self.runfile+'-chain.mc'))
            rm_file(path.join(path1, self.runfile+'-stats.mc'))
            rm_file(path.join(path1, self.runfile+'-point.mc.gz'))
            rm_file(path.join(path1, self.runfile+'-chain.mc.gz'))
            rm_file(path.join(path1, self.runfile+'-stats.mc.gz'))

        if self.jupyter_clear_output:
            clear_output(wait=True)

        return avg_gmm_marginal

