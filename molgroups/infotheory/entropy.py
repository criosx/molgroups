from __future__ import division

__all__ = ["entropy"]

import numpy as np
from molgroups.support import molstat
import pandas
import itertools
from numpy import log, sqrt, log2, pi, e
from numpy.random import permutation
from sklearn.mixture import BayesianGaussianMixture as GMM
import os
from os import path, mkdir
from math import fabs, pow, floor, ceil
from subprocess import Popen
from time import sleep
from IPython.display import clear_output
import shutil

LN2 = log(2)


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
    Paul Kienzle's method to obtain nicely spaced intevals for plot axes
    """

    if step is None:
        step = (stop-start)/numsteps

    sign = 1.0
    if step < 0:
        sign = -1.0
        step = step * (-1)

    if step == 0:
        return [start, stop+1]

    exponent = floor(log(step)/log(10))
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
                 niceticks=False):
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
    plt.savefig(filename + '.pdf')
    plt.close()


def save_plot_2d(x, y, z, xlabel, ylabel, color, filename='plot', zmin=None, zmax=None, levels=20, mark_maximum=False):
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
    cs = ax.contourf(x, y, z, cmap=color, vmin=bounds[0], vmax=bounds[-1], levels=bounds)     # extend='both'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.ticklabel_format(scilimits=(-3, 3), useMathText=True)
    fig.colorbar(cs)

    if mark_maximum:
        index = np.unravel_index(z.argmax(), z.shape)
        plt.text(x[index[1]], y[index[0]], 'x', horizontalalignment='center', verticalalignment='center')

    plt.savefig(path.join('plots', filename) + '.pdf')
    plt.close()


def cov_entropy(C):
    """
    Entropy estimate from covariance matrix C
    """
    return 0.5 * (len(C) * log2(2*pi*e) + log2(abs(np.linalg.det(C))))


class MVNEntropy(object):
    def __init__(self, x, alpha=0.05, max_points=1000):
        # max points was 1000
        # compute Mardia test coefficient
        n, p = x.shape   # num points, num dimensions
        self.mu = np.mean(x, axis=0)
        self.C = np.cov(x.T, bias=True) if p > 1 else np.array([[np.var(x.T, ddof=1)]])

    def entropy(self):
        return cov_entropy(self.C)

    def marginal_entropy(self, independent_pars):
        MC = self.C.copy()
        MC = np.delete(MC, independent_pars, 0)
        MC = np.delete(MC, independent_pars, 1)
        return cov_entropy(MC)


class GMMEntropy(object):
    def __init__(self, x):
        self.n_components = int(5 * sqrt(x.shape[1]))
        self.gmmpredictor = GMM(n_components=self.n_components, max_iter=10000)
        self.gmmpredictor.fit(x)

    def entropy(self, n_est):
        sample, _ = self.gmmpredictor.sample(n_est)
        ll = self.gmmpredictor.score_samples(sample)
        h = -np.mean(ll) / LN2
        return h

    def score_samples(self, x):
        return self.gmmpredictor.score_samples(x)


# calculates entropy while varying a set of parameters in parlist and
# keeping others fixed as specified in simpar.dat
# requires a compiled and ready to go fit whose fit parameters are modified and fixed
# avoid_symmetric prevents calculating symmetry-related results by enforcing the indices of varied parameters
# are an ordered list
class Entropy:

    def __init__(self, fitsource, spath, mcmcpath, runfile, mcmcburn=16000, mcmcsteps=5000, deldir=True,
                 convergence=2.0, miniter=1, mode='water', background_rule=None, bFetchMode=False, bClusterMode=False,
                 calc_symmetric=True, upper_info_plotlevel=None, plotlimits_filename='', slurmscript='',
                 configuration=None):

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

        self.molstat = molstat.CMolStat(fitsource=fitsource, spath=spath, mcmcpath=mcmcpath, runfile=runfile)

        # Parse parameters from entropypar.dat
        header_names = ['type', 'dataset', 'configuration', 'par', 'value', 'l_fit', 'u_fit', 'l_sim', 'u_sim',
                        'step_sim']
        self.allpar = pandas.read_csv('entropypar.dat', sep='\s+', header=None, names=header_names,
                                      skip_blank_lines=True, comment='#')

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
        self.simpar = self.allpar.loc[:, ['par', 'value', 'dataset', 'configuration']]

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

    def calc_entropy(self, molstat=None):
        if molstat is None:
            molstat = self.molstat

        N_entropy = 10000  # was 10000
        N_norm = 10000  # was 2500

        # read MCMC result and save sErr.dat
        points, parnames_dict_keys, logp = molstat.Interactor.fnLoadMCMCResults()
        independent_pars = []
        dependent_pars = []
        parnames = [key for key in parnames_dict_keys]
        for index in range(len(parnames)):
            if parnames[index] in self.independent_parameters:
                independent_pars.append(index)
            else:
                dependent_pars.append(index)

        # Do statistics over points
        points_median = np.median(points, axis=0)
        points_std = np.std(points, axis=0)

        # Use a random subset to estimate density
        if N_norm >= len(logp):
            norm_points = points
        else:
            idx = permutation(len(points))[:N_norm]
            norm_points = points[idx]

        mvnentropy = MVNEntropy(norm_points)
        mvn_entropy = mvnentropy.entropy()
        mvn_entropy_marginal = mvnentropy.marginal_entropy(independent_pars=independent_pars)

        gmmentropy = GMMEntropy(norm_points)
        gmm_entropy = gmmentropy.entropy(N_entropy)
        gmmentropymarginal = GMMEntropy(np.delete(norm_points, independent_pars, 1))
        gmm_entropy_marginal = gmmentropymarginal.entropy(N_entropy)

        # Use a different subset to estimate the scale factor between density
        # and logp.
        if N_entropy >= len(logp):
            entropy_points, eval_logp = points, logp
        else:
            idx = permutation(len(points))[:N_entropy]
            entropy_points, eval_logp = points[idx], logp[idx]

        # This is the original implementation of the KDN entropy from Kramer et al.
        # It is not very stable and replaced by the GMM entropy instead.
        '''
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

    # calculates prior entropy
    def calc_prior(self):
        priorentropy = 0
        priorentropy_marginal = 0
        for row in self.allpar.itertuples():  # cycle through all parameters
            if row.type == 'd' or row.type == 'fd' or row.type == 'i' or row.type == 'fi':
                priorentropy += log(row.u_fit - row.l_fit) / log(2)
                # calculate prior entropy for parameters to be marginalized (dependent parameters)
                if row.type == 'd' or row.type == 'fd':
                    priorentropy_marginal += log(row.u_fit - row.l_fit) / log(2)

        return priorentropy, priorentropy_marginal

    def load_results(self, dirname):
        path1 = path.join(dirname, 'results')
        self.results_mvn = np.load(path.join(path1, 'MVN_entropy.npy'))
        self.results_gmm = np.load(path.join(path1, 'GMM_entropy.npy'))
        self.results_mvn_marginal = np.load(path.join(path1, 'MVN_entropy_marginal.npy'))
        self.results_gmm_marginal = np.load(path.join(path1, 'GMM_entropy_marginal.npy'))
        if path.isfile(path.join(path1, 'Prediction_gpcam.npy')):
            self.prediction_gpcam = np.load(path.join(path1, 'Prediction_gpcam.npy'))
        self.n_mvn = np.load(path.join(path1, 'MVN_n.npy'))
        self.n_gmm = np.load(path.join(path1, 'GMM_n.npy'))
        self.n_mvn_marginal = np.load(path.join(path1, 'MVN_n_marginal.npy'))
        self.n_gmm_marginal = np.load(path.join(path1, 'GMM_n_marginal.npy'))
        self.sqstd_mvn = np.load(path.join(path1, 'MVN_sqstd.npy'))
        self.sqstd_gmm = np.load(path.join(path1, 'GMM_sqstd.npy'))
        self.sqstd_mvn_marginal = np.load(path.join(path1, 'MVN_sqstd_marginal.npy'))
        self.sqstd_gmm_marginal = np.load(path.join(path1, 'GMM_sqstd_marginal.npy'))
        self.par_median = np.load(path.join(path1, 'par_median.npy'))
        self.par_std = np.load(path.join(path1, 'par_std.npy'))

    def runmcmc(self, molstat_instance, iteration, dirname, fulldirname):
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
            molstat_instance.Interactor.fnRunMCMC(burn=self.mcmcburn, steps=self.mcmcsteps, batch=True)
        return

    def plot_results(self, mark_maximum=False):
        import matplotlib.pyplot as plt

        if not path.isdir('plots'):
            mkdir('plots')

        onecolormaps = [plt.cm.Greys, plt.cm.Purples, plt.cm.Blues, plt.cm.Greens, plt.cm.Oranges, plt.cm.Reds]
        ec = plt.cm.coolwarm

        path1 = path.join(self.spath, 'plots')

        if self.plotlimits_filename != '':
            plotlimits = True
            plotlim = pandas.read_csv(self.plotlimits_filename, sep="\s+", header=None)
            plotlim.columns = ['name', 'value']

            if not plotlim.loc[plotlim['name'] == 'info'].empty:
                self.upper_info_plotlevel = float(plotlim.loc[plotlim['name'] == 'info']['value'])
        else:
            plotlimits = False
            plotlim = None

        if len(self.steplist) == 1:
            ax0 = self.axes[0]
            sp0 = self.steppar['par'].tolist()[0]
            save_plot_1d(ax0, self.results_mvn, np.sqrt(self.sqstd_mvn), sp0, 'Entropy [bits]',
                         filename=path.join(path1, 'MVN_entropy'))
            save_plot_1d(ax0, self.results_gmm, np.sqrt(self.sqstd_gmm), sp0, 'Entropy [bits]',
                         filename=path.join(path1, 'GMM_entropy'))
            save_plot_1d(ax0, self.results_mvn_marginal, np.sqrt(self.sqstd_mvn_marginal), sp0, 'Entropy [bits]',
                         filename=path.join(path1, 'MVN_entropy_marginal'))
            save_plot_1d(ax0, self.results_gmm_marginal, np.sqrt(self.sqstd_gmm_marginal), sp0, 'Entropy [bits]',
                         filename=path.join(path1, 'GMM_entropy_marginal'))
            save_plot_1d(ax0, self.priorentropy - self.results_mvn, np.sqrt(self.sqstd_mvn), sp0,
                         'information gain [bits]', filename=path.join(path1, 'MVN_infocontent'), ymin=0)
            save_plot_1d(ax0, self.priorentropy - self.results_gmm, np.sqrt(self.sqstd_gmm), sp0,
                         'information gain [bits]', filename=path.join(path1, 'GMM_infocontent'), ymin=0)
            save_plot_1d(ax0, self.priorentropy_marginal - self.results_mvn_marginal, np.sqrt(self.sqstd_mvn_marginal),
                         sp0, 'information gain [bits]', filename=path.join(path1, 'MVN_infocontent_marginal'), ymin=0,
                         ymax=self.upper_info_plotlevel)
            save_plot_1d(ax0, self.priorentropy_marginal - self.results_gmm_marginal, np.sqrt(self.sqstd_gmm_marginal),
                         sp0, 'information gain [bits]', filename=path.join(path1, 'GMM_infocontent_marginal'), ymin=0,
                         ymax=self.upper_info_plotlevel)
            save_plot_1d(ax0, self.n_mvn, None, sp0, 'computations', filename=path.join(path1, 'MVN_n'), ymin=0)
            save_plot_1d(ax0, self.n_gmm, None, sp0, 'computations', filename=path.join(path1, 'GMM_n'), ymin=0)
            save_plot_1d(ax0, self.n_mvn_marginal, None, sp0, 'computations',
                         filename=path.join(path1, 'MVN_n_marginal'), ymin=0)
            save_plot_1d(ax0, self.n_gmm_marginal, None, sp0, 'computations',
                         filename=path.join(path1, 'GMM_n_marginal'), ymin=0)

            i = 0
            for parname in self.parlist:
                save_plot_1d(ax0, self.par_median[i], self.par_std[i], sp0, parname,
                             filename=path.join(path1, 'Par_' + parname + '_median'))
                i += 1

        elif len(self.steplist) == 2:
            # numpy array and plot axes are reversed
            ax1 = self.axes[0]
            ax0 = self.axes[1]
            sp1 = self.steppar['par'].tolist()[0]
            sp0 = self.steppar['par'].tolist()[1]

            save_plot_2d(ax0, ax1, self.results_mvn, sp0, sp1, ec, filename=path.join(path1, 'MVN_entropy'))
            save_plot_2d(ax0, ax1, self.results_gmm, sp0, sp1, ec, filename=path.join(path1, 'GMM_entropy'))
            save_plot_2d(ax0, ax1, self.results_mvn_marginal, sp0, sp1, ec,
                         filename=path.join(path1, 'MVN_entropy_marginal'))
            save_plot_2d(ax0, ax1, self.results_gmm_marginal, sp0, sp1, ec,
                         filename=path.join(path1, 'GMM_entropy_marginal'))
            save_plot_2d(ax0, ax1, self.priorentropy - self.results_mvn, sp0, sp1, ec,
                         filename=path.join(path1, 'MVN_infocontent'), zmin=0, mark_maximum=mark_maximum)
            save_plot_2d(ax0, ax1, self.priorentropy - self.results_gmm, sp0, sp1, ec,
                         filename=path.join(path1, 'GMM_infocontent'), zmin=0, mark_maximum=mark_maximum)
            save_plot_2d(ax0, ax1, self.priorentropy_marginal - self.results_mvn_marginal, sp0, sp1, ec,
                         filename=path.join(path1, 'MVN_infocontent_marginal'), zmin=0, zmax=self.upper_info_plotlevel
                         , mark_maximum=mark_maximum)
            save_plot_2d(ax0, ax1, self.priorentropy_marginal - self.results_gmm_marginal, sp0, sp1, ec,
                         filename=path.join(path1, 'GMM_infocontent_marginal'), zmin=0, zmax=self.upper_info_plotlevel
                         , mark_maximum=mark_maximum)
            save_plot_2d(ax0, ax1, self.prediction_gpcam, sp0, sp1, ec, filename=path.join(path1, 'Prediction_gpcam'),
                         zmin=0, zmax=self.upper_info_plotlevel, mark_maximum=mark_maximum)
            save_plot_2d(ax0, ax1, self.sqstd_mvn, sp0, sp1, ec, filename=path.join(path1, 'MVN_sqstd'), zmin=0)
            save_plot_2d(ax0, ax1, self.sqstd_gmm, sp0, sp1, ec, filename=path.join(path1, 'GMM_sqstd'), zmin=0)
            save_plot_2d(ax0, ax1, self.sqstd_mvn_marginal, sp0, sp1, ec,
                         filename=path.join(path1, 'MVN_sqstd_marginal'), zmin=0)
            save_plot_2d(ax0, ax1, self.sqstd_gmm_marginal, sp0, sp1, ec,
                         filename=path.join(path1, 'GMM_sqstd_marginal'), zmin=0)
            save_plot_2d(ax0, ax1, self.n_mvn, sp0, sp1, ec, filename=path.join(path1, 'MVN_n'), zmin=0)
            save_plot_2d(ax0, ax1, self.n_gmm, sp0, sp1, ec, filename=path.join(path1, 'GMM_n'), zmin=0)
            save_plot_2d(ax0, ax1, self.n_mvn_marginal, sp0, sp1, ec, filename=path.join(path1, 'MVN_n_marginal'),
                         zmin=0)
            save_plot_2d(ax0, ax1, self.n_gmm_marginal, sp0, sp1, ec, filename=path.join(path1, 'GMM_n_marginal'),
                         zmin=0)

            i = 0
            j = 0
            for parname in self.parlist:
                if plotlimits and not plotlim.loc[plotlim['name'] == parname].empty:
                    zmax = float(plotlim.loc[plotlim['name'] == parname]['value'])
                else:
                    zmax = None
                save_plot_2d(ax0, ax1, self.par_median[i], sp0, sp1, onecolormaps[j],
                             filename=path.join(path1, 'Par_' + parname + '_median'))
                save_plot_2d(ax0, ax1, self.par_std[i], sp0, sp1, onecolormaps[j],
                             filename=path.join(path1, 'Par_' + parname + '_std'), zmin=0, zmax=zmax)
                i += 1
                j += 1
                if j == len(onecolormaps):
                    j = 0

        elif len(self.steplist) == 3 and self.results_gmm.shape[0] < 6:
            ax2 = self.axes[1]
            ax1 = self.axes[2]
            sp2 = self.steppar['par'].tolist()[1]
            sp1 = self.steppar['par'].tolist()[2]
            for slice in range(self.results_gmm.shape[0]):
                save_plot_2d(ax1, ax2, self.results_mvn[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_entropy_' + str(slice)))
                save_plot_2d(ax1, ax2, self.results_gmm[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'GMM_entropy_' + str(slice)))
                save_plot_2d(ax1, ax2, self.results_gmm_marginal[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'GMM_entropy_marginal_' + str(slice)))
                save_plot_2d(ax1, ax2, self.results_mvn_marginal[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_entropy_marginal_' + str(slice)))
                save_plot_2d(ax1, ax2, self.priorentropy - self.results_mvn[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_infocontent_' + str(slice)), zmin=0,
                             mark_maximum=mark_maximum)
                save_plot_2d(ax1, ax2, self.priorentropy - self.results_gmm[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'GMM_infocontent_' + str(slice)), zmin=0,
                             mark_maximum=mark_maximum)
                save_plot_2d(ax1, ax2, self.priorentropy_marginal - self.results_mvn_marginal[slice], sp1, sp2,
                             ec, filename=path.join(path1, 'MVN_infocontent_marginal_' + str(slice)), zmin=0,
                             mark_maximum=mark_maximum)
                save_plot_2d(ax1, ax2, self.priorentropy_marginal - self.results_gmm_marginal[slice], sp1, sp2,
                             ec, filename=path.join(path1, 'GMM_infocontent_marginal_' + str(slice)), zmin=0,
                             mark_maximum=mark_maximum)
                save_plot_2d(ax1, ax2, self.sqstd_mvn[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_sqstd_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.sqstd_gmm[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'GMM_sqstd_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.sqstd_mvn_marginal[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_sqstd_marginal_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.sqstd_gmm_marginal[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'GMM_sqstd_marginal_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.n_mvn[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_n_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.n_gmm[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'GMM_n_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.n_mvn_marginal[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_n_marginal_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.n_gmm_marginal[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'GMM_n_marginal_' + str(slice)), zmin=0)

        if len(self.steplist) >= 3:
            # plot projections onto two parameters at a time
            for i in range(len(self.steppar)):
                for j in range(i):
                    ax2 = self.axes[i]
                    ax1 = self.axes[j]
                    sp2 = self.steppar['par'].tolist()[i]
                    sp1 = self.steppar['par'].tolist()[j]
                    projection_gmm_marginal = np.empty((self.steplist[i], self.steplist[j]))
                    projection_gpcam_prediction = np.empty((self.steplist[i], self.steplist[j]))
                    for k in range(self.steplist[i]):
                        for l in range(self.steplist[j]):
                            projection_gmm_marginal[k, l] = np.take(np.take(self.priorentropy_marginal -
                                                                            self.results_gmm_marginal, indices=k,
                                                                            axis=i), indices=l, axis=j).max()
                            projection_gpcam_prediction[k, l] = np.take(np.take(self.prediction_gpcam, indices=k,
                                                                            axis=i), indices=l, axis=j).max()
                    save_plot_2d(ax1, ax2, projection_gmm_marginal, sp1, sp2, ec,
                                 filename=path.join(path1, 'GMM_n_marginal_projection' + sp1 + '_' + sp2), zmin=0,
                                 mark_maximum=mark_maximum)
                    save_plot_2d(ax1, ax2, projection_gpcam_prediction, sp1, sp2, ec,
                                 filename=path.join(path1, 'gpCAM_prediction_projection' + sp1 + '_' + sp2), zmin=0,
                                 mark_maximum=mark_maximum)

    def save_results(self, dirname):
        path1 = path.join(dirname, 'results')
        if not path.isdir(path1):
            mkdir(path1)
        np.save(path.join(path1, 'GMM_entropy'), self.results_gmm, allow_pickle=False)
        np.save(path.join(path1, 'MVN_entropy'), self.results_mvn, allow_pickle=False)
        np.save(path.join(path1, 'GMM_entropy_marginal'), self.results_gmm_marginal, allow_pickle=False)
        np.save(path.join(path1, 'MVN_entropy_marginal'), self.results_mvn_marginal, allow_pickle=False)
        np.save(path.join(path1, 'GMM_infocontent'), self.priorentropy - self.results_gmm, allow_pickle=False)
        np.save(path.join(path1, 'MVN_infocontent'), self.priorentropy - self.results_mvn, allow_pickle=False)
        np.save(path.join(path1, 'GMM_infocontent_marginal'), self.priorentropy_marginal - self.results_gmm_marginal,
                allow_pickle=False)
        np.save(path.join(path1, 'MVN_infocontent_marginal'), self.priorentropy_marginal - self.results_mvn_marginal,
                allow_pickle=False)
        np.save(path.join(path1, 'Prediction_gpcam'), self.prediction_gpcam, allow_pickle=False)
        np.save(path.join(path1, 'GMM_sqstd'), self.sqstd_gmm, allow_pickle=False)
        np.save(path.join(path1, 'MVN_sqstd'), self.sqstd_mvn, allow_pickle=False)
        np.save(path.join(path1, 'GMM_sqstd_marginal'), self.sqstd_gmm_marginal, allow_pickle=False)
        np.save(path.join(path1, 'MVN_sqstd_marginal'), self.sqstd_mvn_marginal, allow_pickle=False)
        np.save(path.join(path1, 'GMM_n'), self.n_gmm, allow_pickle=False)
        np.save(path.join(path1, 'MVN_n'), self.n_mvn, allow_pickle=False)
        np.save(path.join(path1, 'GMM_n_marginal'), self.n_gmm_marginal, allow_pickle=False)
        np.save(path.join(path1, 'MVN_n_marginal'), self.n_mvn_marginal, allow_pickle=False)
        np.save(path.join(path1, 'par_median'), self.par_median, allow_pickle=False)
        np.save(path.join(path1, 'par_std'), self.par_std, allow_pickle=False)

        # save to txt when not more than two-dimensional array
        if len(self.steplist) <= 2:
            np.savetxt(path.join(path1, 'MVN_entropy.txt'), self.results_mvn)
            np.savetxt(path.join(path1, 'GMM_entropy.txt'), self.results_gmm)
            np.savetxt(path.join(path1, 'MVN_entropy_marginal.txt'), self.results_mvn_marginal)
            np.savetxt(path.join(path1, 'GMM_entropy_marginal.txt'), self.results_gmm_marginal)
            np.savetxt(path.join(path1, 'MVN_infocontent.txt'), self.priorentropy - self.results_mvn)
            np.savetxt(path.join(path1, 'GMM_infocontent.txt'), self.priorentropy - self.results_gmm)
            np.savetxt(path.join(path1, 'MVN_infocontent_marginal.txt'), self.priorentropy_marginal -
                       self.results_mvn_marginal)
            np.savetxt(path.join(path1, 'GMM_infocontent_marginal.txt'), self.priorentropy_marginal -
                       self.results_gmm_marginal)
            np.savetxt(path.join(path1, 'Prediction_gpcam.txt'), self.prediction_gpcam)
            np.savetxt(path.join(path1, 'MVN_sqstd.txt'), self.sqstd_mvn)
            np.savetxt(path.join(path1, 'GMM_sqstd.txt'), self.sqstd_gmm)
            np.savetxt(path.join(path1, 'MVN_sqstd_marginal.txt'), self.sqstd_mvn_marginal)
            np.savetxt(path.join(path1, 'GMM_sqstd_marginal.txt'), self.sqstd_gmm_marginal)
            np.savetxt(path.join(path1, 'MVN_n.txt'), self.n_mvn)
            np.savetxt(path.join(path1, 'GMM_n.txt'), self.n_gmm)
            np.savetxt(path.join(path1, 'MVN_n_marginal.txt'), self.n_mvn_marginal)
            np.savetxt(path.join(path1, 'GMM_n_marginal.txt'), self.n_gmm_marginal)
            i = 0
            for parname in self.parlist:
                np.savetxt(path.join(path1, 'Par_' + parname + '_median.txt'), self.par_median[i])
                np.savetxt(path.join(path1, 'Par_' + parname + '_std.txt'), self.par_std[i])
                i += 1

        # save three-dimensional array in slices of the first parameter
        if len(self.steplist) == 3 and self.results_gmm.shape[0] < 6:
            for slice in range(self.results_gmm.shape[0]):
                np.savetxt(path.join(path1, 'MVN_entropy_' + str(slice) + '.txt'), self.results_mvn[slice])
                np.savetxt(path.join(path1, 'GMM_entropy_' + str(slice) + '.txt'), self.results_gmm[slice])
                np.savetxt(path.join(path1, 'MVN_entropy_marginal_' + str(slice) + '.txt'),
                           self.results_mvn_marginal[slice])
                np.savetxt(path.join(path1, 'GMM_entropy_marginal_' + str(slice) + '.txt'),
                           self.results_gmm_marginal[slice])
                np.savetxt(path.join(path1, 'MVN_infocontent_' + str(slice) + '.txt'), self.priorentropy -
                           self.results_mvn[slice])
                np.savetxt(path.join(path1, 'GMM_infocontent_' + str(slice) + '.txt'), self.priorentropy -
                           self.results_gmm[slice])
                np.savetxt(path.join(path1, 'MVN_infocontent_marginal_' + str(slice) + '.txt'),
                           self.priorentropy_marginal - self.results_mvn_marginal[slice])
                np.savetxt(path.join(path1, 'GMM_infocontent_marginal_' + str(slice) + '.txt'),
                           self.priorentropy_marginal - self.results_gmm_marginal[slice])
                np.savetxt(path.join(path1, 'MVN_sqstd_' + str(slice) + '.txt'), self.sqstd_mvn[slice])
                np.savetxt(path.join(path1, 'GMM_sqstd_' + str(slice) + '.txt'), self.sqstd_gmm[slice])
                np.savetxt(path.join(path1, 'MVN_sqstd_marginal_' + str(slice) + '.txt'),
                           self.sqstd_mvn_marginal[slice])
                np.savetxt(path.join(path1, 'GMM_sqstd_marginal_' + str(slice) + '.txt'),
                           self.sqstd_gmm_marginal[slice])
                np.savetxt(path.join(path1, 'MVN_n_' + str(slice) + '.txt'), self.n_mvn[slice])
                np.savetxt(path.join(path1, 'GMM_n_' + str(slice) + '.txt'), self.n_gmm[slice])
                np.savetxt(path.join(path1, 'MVN_n_marginal_' + str(slice) + '.txt'), self.n_mvn_marginal[slice])
                np.savetxt(path.join(path1, 'GMM_n_marginal_' + str(slice) + '.txt'), self.n_gmm_marginal[slice])

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

    def run_optimization(self, qmin=None, qmax=None, qrangefromfile=False, t_total=None, optimizer='grid',
                         jupyter_clear_output=False, gpcam_iterations=50, retrain_async_at=(15, 30)):

        def writeout_result(_itindex, avg_mvn, avg_gmm, avg_mvn_marginal, avg_gmm_marginal, points_median, points_std,
                            parnames):
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

        def calc_entropy_for_index(molstat_iter, itindex):
            # Calculate Entropy n times and average
            mvn_entropy = []
            gmm_entropy = []
            mvn_entropy_marginal = []
            gmm_entropy_marginal = []
            points_median = []
            points_std = []
            parnames = []

            # GMM and MVN are much more stable than KVN. Averaging over several calculations appears
            # not to be necessary
            for j in range(1):  # was 10
                # calculate entropy, dependent parameters == parameters of interest
                # independent parameters == nuisance parameters
                a, b, c, d, points_median, points_std, parnames = self.calc_entropy(molstat_iter)
                mvn_entropy.append(a)
                gmm_entropy.append(b)
                mvn_entropy_marginal.append(c)
                gmm_entropy_marginal.append(d)

            # remove outliers and average calculated entropies, don't average over parameter stats
            avg_mvn, std_mvn = average(mvn_entropy)
            avg_gmm, std_gmm = average(gmm_entropy)
            avg_mvn_marginal, std_mvn_marginal = average(mvn_entropy_marginal)
            avg_gmm_marginal, std_gmm_marginal = average(gmm_entropy_marginal)

            bValidResult = (std_gmm < self.convergence) and \
                           (self.priorentropy_marginal - avg_gmm_marginal > (-0.5) * len(
                               self.dependent_parameters)) and \
                           (self.priorentropy - avg_gmm > (-0.5) * len(self.parlist))

            # no special treatment for first entry necessary, algorithm catches this
            if bValidResult:
                writeout_result(itindex, avg_mvn, avg_gmm, avg_mvn_marginal, avg_gmm_marginal, points_median,
                                points_std, parnames)

            # save results for every iteration and delete large files
            self.save_results(self.spath)

            return avg_gmm_marginal

        def set_sim_pars_for_index(it):
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

                # search if there is a designated background parameter that takes the cb instead of a configuration
                bFoundBackground = False
                for row in self.allpar.itertuples():
                    if row.dataset == 'b'+str(dset) or row.dataset == 'b' or row.dataset == 'b*':
                        self.simpar.loc[self.simpar['par'] == row.par, 'value'] = cb
                        bFoundBackground = True

                if bFoundBackground:
                    return configurations

                # change buffer crosssection in configurations
                configurations = _fill_config(configurations, 'differential_cross_section_buffer', cb, dset, config)

                return configurations

            # Configurations are imported externally, if not empty configuration initialization here
            # In any case, missing parameters are set to the default in the API simulation routines
            # TODO: Properly implement 'mode' functionality

            configurations = self.configuration
            if configurations is None:
                configurations = [[{}]]

            # cycle through all parameters
            isim = 0
            for row in self.allpar.itertuples():
                # is it a parameter to iterate over?
                if row.par in self.steppar['par'].tolist():

                    lsim = self.steppar.loc[self.steppar['par'] == row.par, 'l_sim'].iloc[0]
                    stepsim = self.steppar.loc[self.steppar['par'] == row.par, 'step_sim'].iloc[0]
                    value = self.steppar.loc[self.steppar['par'] == row.par, 'value'].iloc[0]
                    lfit = self.steppar.loc[self.steppar['par'] == row.par, 'l_fit'].iloc[0]
                    ufit = self.steppar.loc[self.steppar['par'] == row.par, 'u_fit'].iloc[0]
                    # TODO: Here we might implicitly assume an order. Check and resolve.
                    simvalue = lsim + stepsim * it.multi_index[isim]

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
                            self.simpar.loc[self.simpar['par'] == row.par, 'value'] = simvalue
                            self.molstat.Interactor.fnReplaceParameterLimitsInSetup(row.par, lowersim, uppersim)

                    else:
                        # must be instrument parameter
                        configurations = _fill_config(configurations, row.par, simvalue, row.dataset, row.configuration)
                        for indx in self.simpar.index:
                            if self.simpar['par'][indx] == row.par and self.simpar['dataset'][indx] == row.dataset and \
                                    self.simpar['configuration'][indx] == row.configuration:
                                self.simpar.iloc[indx, self.simpar.columns.get_loc('value')] = simvalue
                                break
                    isim += 1

                elif row.type == 'd' or row.type == 'fd' or row.type == 'i' or row.type == 'fi':
                    if 'b' not in row.dataset:
                        # only change boundaries when it will not be filled by a background rule
                        self.molstat.Interactor.fnReplaceParameterLimitsInSetup(row.par, row.l_fit, row.u_fit)

                else:
                    # must be instrument parameter
                    configurations = _fill_config(configurations, row.par, row.value, row.dataset, row.configuration)

                if (row.dataset != 'n') and (row.dataset != '_') and ('b' not in row.dataset):
                    # this is a parameter that will determine one or more isotropic backgrounds
                    configurations = _set_background(configurations, row.dataset, row.configuration, row.value)

            simparsave = self.simpar.loc[:, ['par', 'value']]
            simparsave.to_csv(path.join(self.spath, 'simpar.dat'), sep=' ', header=None, index=False)
            return configurations

        def work_on_index(it):
            # Recover itindex from 'it'. This is not an argument to the function anymore, as work_on_index might
            # be used independently of iterate_over_all_indices
            itindex = it.multi_index

            # calculate which iteration we are on from itindex
            if self.calc_symmetric:
                iteration = it.iterindex
            else:
                # TODO: This is potentially expensive. Find better method for symmetry-concious calculation
                it2 = np.nditer(self.results_gmm, flags=['multi_index'])
                iteration = 0
                while not it2.finished:
                    itindex2 = it2.multi_index
                    if all(itindex2[i] <= itindex2[i + 1] for i in range(len(itindex2) - 1)):
                        # iterations are only increased if this index is not dropped because of symmetry
                        iteration += 1
                    it2.iternext()

            dirname = 'iteration_' + str(iteration)
            fulldirname = path.join(self.spath, dirname)
            path1 = path.join(fulldirname, 'save')
            chainname = path.join(path1, self.runfile+'-chain.mc')

            # most relevant result for a particular index to return for general use of this function
            avg_gmm_marginal = 0

            # fetch mode and cluster mode are exclusive
            if not self.bFetchMode:
                # run a new fit, preparations are done in the root directory and the new fit is copied into the
                # iteration directory, preparations in the iterations directory are not possible, because it would
                # be lacking a result directory, which is needed for restoring a state/parameters
                molstat_iter = molstat.CMolStat(fitsource=self.fitsource, spath=fulldirname, mcmcpath='save',
                                                runfile=self.runfile, load_state=False)
                self.molstat.Interactor.fnBackup(target=path.join(self.spath, 'simbackup'))
                configurations = set_sim_pars_for_index(it)
                self.molstat.fnSimulateData(mode=self.mode, liConfigurations=configurations, qmin=qmin, qmax=qmax,
                                            qrangefromfile=qrangefromfile, t_total=t_total)
                self.molstat.Interactor.fnBackup(origin=self.spath, target=fulldirname)
                # previous save needs to be removed as output serves as flag for HPC job termination
                if path.isdir(path1):
                    shutil.rmtree(path1)
                self.molstat.Interactor.fnRemoveBackup(target=path.join(self.spath, 'simbackup'))
                self.runmcmc(molstat_iter, iteration, dirname, fulldirname)

            # Do not run entropy calculation when on cluster.
            bPriorResultExists = path.isfile(chainname) or path.isfile(chainname + '.gz')
            if not self.bClusterMode and bPriorResultExists:
                molstat_iter = molstat.CMolStat(fitsource=self.fitsource, spath=fulldirname, mcmcpath='save',
                                                runfile=self.runfile)
                avg_gmm_marginal = calc_entropy_for_index(molstat_iter, itindex)

            # delete big files except in Cluster mode. They are needed there for future fetching
            if self.deldir and not self.bClusterMode:
                rm_file(path.join(path1, self.runfile+'-point.mc'))
                rm_file(path.join(path1, self.runfile+'-chain.mc'))
                rm_file(path.join(path1, self.runfile+'-stats.mc'))
                rm_file(path.join(path1, self.runfile+'-point.mc.gz'))
                rm_file(path.join(path1, self.runfile+'-chain.mc.gz'))
                rm_file(path.join(path1, self.runfile+'-stats.mc.gz'))

            if jupyter_clear_output:
                clear_output(wait=True)

            return avg_gmm_marginal

        def iterate_over_all_indices(refinement=False):
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
                        _ = work_on_index(it)

                it.iternext()
            return bWorkedOnIndex
        if optimizer == 'grid':
            # Grid search
            # every index has at least one result before re-analyzing any data point (refinement)
            bRefinement = False
            while True:
                bWorkedOnAnyIndex = iterate_over_all_indices(bRefinement)
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
        elif optimizer == 'gpCAM' or optimizer == 'gpcam':
            # Using the gpCAM global optimizer
            # follows the example from the gpCAM website
            import time
            from gpcam.autonomous_experimenter import AutonomousExperimenterGP

            def instrument(data):
                print("This is the current length of the data received by gpCAM: ", len(data))
                for entry in data:
                    # convert position into grid indices
                    position_index = []
                    for ax in range(len(self.axes)):
                        list1 = [np.array(self.axes[ax]) - entry['position'][ax]]
                        arr1 = np.array(list1)
                        arr1 = np.abs(arr1)
                        position_index.append(np.argmin(arr1))
                    # convert grid indices into iterator
                    # TODO: There must be a more effient way to set an iterator by an index
                    it = np.nditer(self.results_gmm, flags=['multi_index'])
                    while True:
                        if np.array_equal(it.multi_index, position_index):
                            break
                        it.iternext()
                    # We constrain ourselves to the initial grid, but we do not let gpCAM know (might change in future)
                    # entry["position"] = position_index
                    _ = work_on_index(it)
                    entry["value"] = self.priorentropy_marginal - self.results_gmm_marginal[it.multi_index]
                    var = self.sqstd_gmm_marginal[it.multi_index]
                    if var > 0:
                        entry['variance'] = var
                    else:
                        entry['variance'] = 1
                    # entry["cost"]  = [np.array([0,0]),entry["position"],np.sum(entry["position"])]
                    self.plot_results(mark_maximum=True)
                return data

            # initialization
            # feel free to try different acquisition functions, e.g. optional_acq_func, "covariance", "shannon_ig"
            # note how costs are defined in for the autonomous experimenter

            parlimits = []
            for row in self.steppar.iterrows():
                parlimits.append([row[1].l_sim, row[1].u_sim])
            parlimits = np.array(parlimits)
            numpars = len(parlimits)

            it = np.nditer(self.results_gmm_marginal, flags=['multi_index'])
            x = []
            y = []
            v = []
            while not it.finished:
                if (self.priorentropy_marginal - self.results_gmm_marginal[it.multi_index]) != 0.:
                    position = []
                    for i in range(len(self.axes)):
                        position.append(self.axes[i][it.multi_index[i]])
                    x.append(position)
                    y.append(self.priorentropy_marginal - self.results_gmm_marginal[it.multi_index])
                    if self.sqstd_gmm_marginal[it.multi_index] > 0:
                        v.append(self.sqstd_gmm_marginal[it.multi_index])
                    else:
                        v.append(1.0)
                it.iternext()

            if len(x) > 10:
                x = np.array(x)
                y = np.array(y)
                v = np.array(v)
            else:
                x = None
                y = None
                v = None

            hyperpars = np.ones([numpars+1])
            hyper_bounds = np.array([[0.001, 100]] * (numpars+1))

            my_ae = AutonomousExperimenterGP(parlimits, hyperpars, hyper_bounds,
                                             init_dataset_size=20, instrument_func=instrument,
                                             acq_func="variance",  # optional_acq_func,
                                             # cost_func = optional_cost_function,
                                             # cost_update_func = optional_cost_update_function,
                                             x=x, y=y, v=v,
                                             # cost_func_params={"offset": 5.0, "slope": 10.0},
                                             kernel_func=None, use_inv=True,
                                             communicate_full_dataset=False, ram_economy=True)
            # , prior_mean_func = optional_mean_func)

            print("length of the dataset: ", len(my_ae.x))

            # my_ae.train_async()                 #train asynchronously
            my_ae.train(method="global")           # or not, or both, choose between "global","local" and "hgdl"

            # update hyperparameters in case they are optimized asynchronously
            # my_ae.update_hps()

            # training and client can be killed if desired and in case they are optimized asynchronously
            # my_ae.kill_training()

            # here we see how python's help function is used to get info about a function
            # help(my_ae.go)

            # run the autonomous loop
            my_ae.go(N=gpcam_iterations,
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

            # save prediction
            dimension = 1
            for factor in self.steplist:
                dimension *= factor
            prediction_positions = np.zeros((dimension, len(self.steplist)))
            it = np.nditer(self.results_gmm, flags=['multi_index'])
            while not it.finished:
                for i in range(len(self.axes)):
                    prediction_positions[it.iterindex][i] = self.axes[i][it.multi_index[i]]
                it.iternext()

            res = my_ae.gp_optimizer.posterior_mean(prediction_positions)
            f = res["f(x)"]
            self.prediction_gpcam = f.reshape(self.steplist)
            self.plot_results(mark_maximum=True)
            self.save_results(self.spath)
