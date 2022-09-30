from __future__ import division

__all__ = ["entropy"]

import numpy as np
from molgroups.support import molstat
import pandas
import itertools
from numpy import mean, std, exp, log, sqrt, log2, pi, e, ndarray
from numpy.random import permutation
from sklearn.mixture import BayesianGaussianMixture as GMM
import os
from os import path, mkdir
from math import fabs, pow, floor, ceil
from subprocess import Popen
from time import sleep
import shutil

LN2 = log(2)


# static methods
def average(alist):
    notfinished = True
    while notfinished:
        s = np.std(alist)
        result = np.mean(alist)

        maxs = 0
        for element in alist:
            s2 = fabs(result - element)
            if s2 > maxs:
                maxs = s2
                val = element
        if maxs > 3 * s:
            alist.remove(val)
        else:
            notfinished = False

    return result, s


# multidimensional convolution of nearest horizontal, vertical, and diagonal neighbors but not the center
def convolute(a):
    conv_arr = np.zeros(a.shape)
    offsetshape = tuple(3 for i in range(len(a.shape)))
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

    if step is None:
        step = (stop-start)/numsteps

    sign = 1.0
    if step < 0:
        sign = -1.0
        step = step *(-1)

    if step == 0:
        return [start, stop+1]

    exponent = floor(log(step)/log(10))
    mantisse = step / pow(10,exponent)

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


def save_plot_2d(x, y, z, xlabel, ylabel, color, filename='plot', zmin=None, zmax=None, levels=20):
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
        mu = np.mean(x, axis=0)
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
        self.gmmpredictor = GMM(n_components=self.n_components, max_iter=1000)
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

    def __init__(self, fitsource, spath, mcmcpath, runfile, mcmcburn=16000, mcmcsteps=5000, deldir=True, convergence=2.0,
                 miniter=1, mode='water', bFetchMode=False, bClusterMode=False, calc_symmetric=True,
                 upper_info_plotlevel=None, plotlimits_filename='', slurmscript=''):

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
        self.bFetchMode = bFetchMode
        self.bClusterMode = bClusterMode
        self.calc_symmetric = calc_symmetric
        self.upper_info_plotlevel = upper_info_plotlevel
        self.plotlimits_filename = plotlimits_filename
        self.slurmscript = slurmscript

        self.molstat = molstat.CMolStat(fitsource=fitsource, spath=spath, mcmcpath=mcmcpath, runfile=runfile)

        # all parameters from entropypar.dat
        header_names = ['type', 'par', 'value', 'l_fit', 'u_fit', 'l_sim', 'u_sim', 'step_sim']
        self.allpar = pandas.read_csv('entropypar.dat', sep='\s+', header=None, names=header_names,
                                      skip_blank_lines=True, comment='#')

        # identify dependent (a), independent (b), and non-parameters in simpar.dat for the calculation of p(a|b,y)
        # later on
        # it is assumed that all parameters in setup.cc are also specified in simpar.dat in exactly the same order
        # this might have to be looked at in future
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
                self.independent_parameters.append(i)
                self.parlist.append(row.par)
                i += 1
            elif row.type == 'd' or row.type == 'fd':
                self.dependent_parameters.append(i)
                self.parlist.append(row.par)
                i += 1

        # only those parameters that will be varied
        self.steppar = self.allpar.dropna(axis=0)

        # create data frame for simpar.dat needed by the reflectivity simulation routine
        # non-parameters such as qrange and prefactor will be included in simpar, but eventually ignored
        # when simulating the reflectivity, as they will find no counterpart in par.dat
        self.simpar = self.allpar.loc[:, ['par', 'value']]

        self.steplist = []
        self.axes = []
        for row in self.steppar.itertuples():
            steps = int((row.u_sim - row.l_sim) / row.step_sim) + 1
            self.steplist.append(steps)
            axis = []
            for i in range(steps):
                axis.append(row.l_sim + i * row.step_sim)
            self.axes.append(axis)

        if path.isfile(path.join(spath, 'results', 'MVN_entropy.npy')):
            self.load_results(spath)
        else:
            self.results_mvn = np.zeros(self.steplist)
            self.results_kdn = np.zeros(self.steplist)
            self.results_mvn_marginal = np.zeros(self.steplist)
            self.results_kdn_marginal = np.zeros(self.steplist)
            self.n_mvn = np.zeros(self.results_mvn.shape)
            self.n_kdn = np.zeros(self.results_kdn.shape)
            self.n_mvn_marginal = np.zeros(self.results_mvn_marginal.shape)
            self.n_kdn_marginal = np.zeros(self.results_kdn_marginal.shape)
            self.sqstd_mvn = np.zeros(self.results_mvn.shape)
            self.sqstd_kdn = np.zeros(self.results_kdn.shape)
            self.sqstd_mvn_marginal = np.zeros(self.results_mvn_marginal.shape)
            self.sqstd_kdn_marginal = np.zeros(self.results_kdn_marginal.shape)
            self.par_median = np.zeros((len(self.parlist),) + self.results_mvn.shape)
            self.par_std = np.zeros((len(self.parlist),) + self.results_mvn.shape)

        self.priorentropy, self.priorentropy_marginal = self.calc_prior()

    def calc_entropy(self, molstat=None):
        if molstat is None:
            molstat = self.molstat

        N_entropy = 10000  # was 10000
        N_norm = 10000  # was 2500

        # read MCMC result and save sErr.dat
        points, parnames, logp = molstat.Interactor.fnLoadMCMCResults()

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
        mvn_entropy_marginal = mvnentropy.marginal_entropy(independent_pars=self.independent_parameters)

        gmmentropy = GMMEntropy(norm_points)
        gmm_entropy = gmmentropy.entropy(N_entropy)
        gmmentropymarginal = GMMEntropy(np.delete(norm_points, self.independent_parameters, 1))
        gmm_entropy_marginal = gmmentropymarginal.entropy(N_entropy)

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

        dependent_points = entropy_points[:, self.dependent_parameters]
        kdn_entropy_marginal = (-1) * np.mean(gmmentropymarginal.score_samples(dependent_points)) / LN2

        # return MVN entropy, KDN entropy, conditional MVN entropy, conditional KDN entropy
        return gmm_entropy, kdn_entropy, gmm_entropy_marginal, kdn_entropy_marginal, points_median, points_std

    # calculates prior entropy
    def calc_prior(self):
        priorentropy = 0
        priorentropy_marginal = 0
        for row in self.allpar.itertuples():  # cycle through all parameters
            if row.par != 'qrange' and row.par != 'prefactor':
                if 'rho' in row.par:
                    priorentropy += log((row.u_fit - row.l_fit) * 1e6) / log(2)
                else:
                    priorentropy += log(row.u_fit - row.l_fit) / log(2)
                # calculate prior entropy for parameters to be marginalized (dependent parameters)
                if row.type == 'd':
                    if 'rho' in row.par:
                        priorentropy_marginal += log((row.u_fit - row.l_fit) * 1e6) / log(2)
                    else:
                        priorentropy_marginal += log(row.u_fit - row.l_fit) / log(2)

        return priorentropy, priorentropy_marginal

    def load_results(self, dirname):
        path1 = path.join(dirname, 'results')
        self.results_mvn = np.load(path.join(path1, 'MVN_entropy.npy'))
        self.results_kdn = np.load(path.join(path1, 'KDN_entropy.npy'))
        self.results_mvn_marginal = np.load(path.join(path1, 'MVN_entropy_marginal.npy'))
        self.results_kdn_marginal = np.load(path.join(path1, 'KDN_entropy_marginal.npy'))
        self.n_mvn = np.load(path.join(path1, 'MVN_n.npy'))
        self.n_kdn = np.load(path.join(path1, 'KDN_n.npy'))
        self.n_mvn_marginal = np.load(path.join(path1, 'MVN_n_marginal.npy'))
        self.n_kdn_marginal = np.load(path.join(path1, 'KDN_n_marginal.npy'))
        self.sqstd_mvn = np.load(path.join(path1, 'MVN_sqstd.npy'))
        self.sqstd_kdn = np.load(path.join(path1, 'KDN_sqstd.npy'))
        self.sqstd_mvn_marginal = np.load(path.join(path1, 'MVN_sqstd_marginal.npy'))
        self.sqstd_kdn_marginal = np.load(path.join(path1, 'KDN_sqstd_marginal.npy'))
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

    def plot_results(self):
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

        if len(self.steplist) == 1:
            ax0 = self.axes[0]
            sp0 = self.steppar['par'].tolist()[0]
            save_plot_1d(ax0, self.results_mvn, np.sqrt(self.sqstd_mvn), sp0, 'Entropy [bits]',
                         filename=path.join(path1, 'MVN_entropy'))
            save_plot_1d(ax0, self.results_kdn, np.sqrt(self.sqstd_kdn), sp0, 'Entropy [bits]',
                         filename=path.join(path1, 'KDN_entropy'))
            save_plot_1d(ax0, self.results_mvn_marginal, np.sqrt(self.sqstd_mvn_marginal), sp0, 'Entropy [bits]',
                         filename=path.join(path1, 'MVN_entropy_marginal'))
            save_plot_1d(ax0, self.results_kdn_marginal, np.sqrt(self.sqstd_kdn_marginal), sp0, 'Entropy [bits]',
                         filename=path.join(path1, 'KDN_entropy_marginal'))
            save_plot_1d(ax0, self.priorentropy - self.results_mvn, np.sqrt(self.sqstd_mvn), sp0,
                         'information gain [bits]', filename=path.join(path1, 'MVN_infocontent'), ymin=0)
            save_plot_1d(ax0, self.priorentropy - self.results_kdn, np.sqrt(self.sqstd_kdn), sp0,
                         'information gain [bits]', filename=path.join(path1, 'KDN_infocontent'), ymin=0)
            save_plot_1d(ax0, self.priorentropy_marginal - self.results_mvn_marginal, np.sqrt(self.sqstd_mvn_marginal),
                         sp0, 'information gain [bits]', filename=path.join(path1, 'MVN_infocontent_marginal'), ymin=0,
                         ymax=self.upper_info_plotlevel)
            save_plot_1d(ax0, self.priorentropy_marginal - self.results_kdn_marginal, np.sqrt(self.sqstd_kdn_marginal),
                         sp0, 'information gain [bits]', filename=path.join(path1, 'KDN_infocontent_marginal'), ymin=0,
                         ymax=self.upper_info_plotlevel)
            save_plot_1d(ax0, self.n_mvn, None, sp0, 'computations', filename=path.join(path1, 'MVN_n'), ymin=0)
            save_plot_1d(ax0, self.n_kdn, None, sp0, 'computations', filename=path.join(path1, 'KDN_n'), ymin=0)
            save_plot_1d(ax0, self.n_mvn_marginal, None, sp0, 'computations',
                         filename=path.join(path1, 'MVN_n_marginal'), ymin=0)
            save_plot_1d(ax0, self.n_kdn_marginal, None, sp0, 'computations',
                         filename=path.join(path1, 'KDN_n_marginal'), ymin=0)

            i = 0
            for parname in self.parlist:
                save_plot_1d(ax0, self.par_median[i], self.par_std[i], sp0, parname,
                             filename=path.join(path1, 'Par_' +  parname + '_median'))
                i += 1

        elif len(self.steplist) == 2:
            # numpy array and plot axes are reversed
            ax1 = self.axes[0]
            ax0 = self.axes[1]
            sp1 = self.steppar['par'].tolist()[0]
            sp0 = self.steppar['par'].tolist()[1]

            save_plot_2d(ax0, ax1, self.results_mvn, sp0, sp1, ec, filename=path.join(path1, 'MVN_entropy'))
            save_plot_2d(ax0, ax1, self.results_kdn, sp0, sp1, ec, filename=path.join(path1, 'KDN_entropy'))
            save_plot_2d(ax0, ax1, self.results_mvn_marginal, sp0, sp1, ec,
                         filename=path.join(path1, 'MVN_entropy_marginal'))
            save_plot_2d(ax0, ax1, self.results_kdn_marginal, sp0, sp1, ec,
                         filename=path.join(path1, 'KDN_entropy_marginal'))
            save_plot_2d(ax0, ax1, self.priorentropy - self.results_mvn, sp0, sp1, ec,
                         filename=path.join(path1, 'MVN_infocontent'), zmin=0)
            save_plot_2d(ax0, ax1, self.priorentropy - self.results_kdn, sp0, sp1, ec,
                         filename=path.join(path1, 'KDN_infocontent'), zmin=0)
            save_plot_2d(ax0, ax1, self.priorentropy_marginal - self.results_mvn_marginal, sp0, sp1, ec,
                         filename=path.join(path1, 'MVN_infocontent_marginal'), zmin=0, zmax=self.upper_info_plotlevel)
            save_plot_2d(ax0, ax1, self.priorentropy_marginal - self.results_kdn_marginal, sp0, sp1, ec,
                         filename=path.join(path1, 'KDN_infocontent_marginal'), zmin=0, zmax=self.upper_info_plotlevel)
            save_plot_2d(ax0, ax1, self.sqstd_mvn, sp0, sp1, ec, filename=path.join(path1, 'MVN_sqstd'), zmin=0)
            save_plot_2d(ax0, ax1, self.sqstd_kdn, sp0, sp1, ec, filename=path.join(path1, 'KDN_sqstd'), zmin=0)
            save_plot_2d(ax0, ax1, self.sqstd_mvn_marginal, sp0, sp1, ec,
                         filename=path.join(path1, 'MVN_sqstd_marginal'), zmin=0)
            save_plot_2d(ax0, ax1, self.sqstd_kdn_marginal, sp0, sp1, ec,
                         filename=path.join(path1, 'KDN_sqstd_marginal'), zmin=0)
            save_plot_2d(ax0, ax1, self.n_mvn, sp0, sp1, ec, filename=path.join(path1, 'MVN_n'), zmin=0)
            save_plot_2d(ax0, ax1, self.n_kdn, sp0, sp1, ec, filename=path.join(path1, 'KDN_n'), zmin=0)
            save_plot_2d(ax0, ax1, self.n_mvn_marginal, sp0, sp1, ec, filename=path.join(path1, 'MVN_n_marginal'),
                         zmin=0)
            save_plot_2d(ax0, ax1, self.n_kdn_marginal, sp0, sp1, ec, filename=path.join(path1, 'KDN_n_marginal'),
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

        elif len(self.steplist) == 3:
            ax2 = self.axes[1]
            ax1 = self.axes[2]
            sp2 = self.steppar['par'].tolist()[1]
            sp1 = self.steppar['par'].tolist()[2]
            for slice in range(self.results_kdn.shape[0]):
                save_plot_2d(ax1, ax2, self.results_mvn[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_entropy_' + str(slice)))
                save_plot_2d(ax1, ax2, self.results_kdn[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'KDN_entropy_' + str(slice)))
                save_plot_2d(ax1, ax2, self.results_mvn_marginal[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_entropy_marginal_' + str(slice)))
                save_plot_2d(ax1, ax2, self.results_mvn_marginal[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_entropy_marginal_' + str(slice)))
                save_plot_2d(ax1, ax2, self.priorentropy - self.results_mvn[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_infocontent_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.priorentropy - self.results_kdn[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'KDN_infocontent_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.priorentropy_marginal - self.results_mvn_marginal[slice], sp1, sp2,
                             ec, filename=path.join(path1, 'MVN_infocontent_marginal_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.priorentropy_marginal - self.results_kdn_marginal[slice], sp1, sp2,
                             ec, filename=path.join(path1, 'KDN_infocontent_marginal_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.sqstd_mvn[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_sqstd_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.sqstd_kdn[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'KDN_sqstd_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.sqstd_mvn_marginal[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_sqstd_marginal_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.sqstd_kdn_marginal[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'KDN_sqstd_marginal_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.n_mvn[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_n_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.n_kdn[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'KDN_n_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.n_mvn_marginal[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'MVN_n_marginal_' + str(slice)), zmin=0)
                save_plot_2d(ax1, ax2, self.n_kdn_marginal[slice], sp1, sp2, ec,
                             filename=path.join(path1, 'KDN_n_marginal_' + str(slice)), zmin=0)

    def save_results(self, dirname):
        path1 = path.join(dirname, 'results')
        if not path.isdir(path1):
            mkdir(path1)
        np.save(path.join(path1, 'KDN_entropy'), self.results_kdn, allow_pickle=False)
        np.save(path.join(path1, 'MVN_entropy'), self.results_mvn, allow_pickle=False)
        np.save(path.join(path1, 'KDN_entropy_marginal'), self.results_kdn_marginal, allow_pickle=False)
        np.save(path.join(path1, 'MVN_entropy_marginal'), self.results_mvn_marginal, allow_pickle=False)
        np.save(path.join(path1, 'KDN_infocontent'), self.priorentropy - self.results_kdn, allow_pickle=False)
        np.save(path.join(path1, 'MVN_infocontent'), self.priorentropy - self.results_mvn, allow_pickle=False)
        np.save(path.join(path1, 'KDN_infocontent_marginal'), self.priorentropy_marginal - self.results_kdn_marginal,
                allow_pickle=False)
        np.save(path.join(path1, 'MVN_infocontent_marginal'), self.priorentropy_marginal - self.results_mvn_marginal,
                allow_pickle=False)
        np.save(path.join(path1, 'KDN_sqstd'), self.sqstd_kdn, allow_pickle=False)
        np.save(path.join(path1, 'MVN_sqstd'), self.sqstd_mvn, allow_pickle=False)
        np.save(path.join(path1, 'KDN_sqstd_marginal'), self.sqstd_kdn_marginal, allow_pickle=False)
        np.save(path.join(path1, 'MVN_sqstd_marginal'), self.sqstd_mvn_marginal, allow_pickle=False)
        np.save(path.join(path1, 'KDN_n'), self.n_kdn, allow_pickle=False)
        np.save(path.join(path1, 'MVN_n'), self.n_mvn, allow_pickle=False)
        np.save(path.join(path1, 'KDN_n_marginal'), self.n_kdn_marginal, allow_pickle=False)
        np.save(path.join(path1, 'MVN_n_marginal'), self.n_mvn_marginal, allow_pickle=False)
        np.save(path.join(path1, 'par_median'), self.par_median, allow_pickle=False)
        np.save(path.join(path1, 'par_std'), self.par_std, allow_pickle=False)

        # save to txt when not more than two-dimensional array
        if len(self.steplist) <= 2:
            np.savetxt(path.join(path1, 'MVN_entropy.txt'), self.results_mvn)
            np.savetxt(path.join(path1, 'KDN_entropy.txt'), self.results_kdn)
            np.savetxt(path.join(path1, 'MVN_entropy_marginal.txt'), self.results_mvn_marginal)
            np.savetxt(path.join(path1, 'KDN_entropy_marginal.txt'), self.results_kdn_marginal)
            np.savetxt(path.join(path1, 'MVN_infocontent.txt'), self.priorentropy - self.results_mvn)
            np.savetxt(path.join(path1, 'KDN_infocontent.txt'), self.priorentropy - self.results_kdn)
            np.savetxt(path.join(path1, 'MVN_infocontent_marginal.txt'), self.priorentropy_marginal -
                       self.results_mvn_marginal)
            np.savetxt(path.join(path1, 'KDN_infocontent_marginal.txt'), self.priorentropy_marginal -
                       self.results_kdn_marginal)
            np.savetxt(path.join(path1, 'MVN_sqstd.txt'), self.sqstd_mvn)
            np.savetxt(path.join(path1, 'KDN_sqstd.txt'), self.sqstd_kdn)
            np.savetxt(path.join(path1, 'MVN_sqstd_marginal.txt'), self.sqstd_mvn_marginal)
            np.savetxt(path.join(path1, 'KDN_sqstd_marginal.txt'), self.sqstd_kdn_marginal)
            np.savetxt(path.join(path1, 'MVN_n.txt'), self.n_mvn)
            np.savetxt(path.join(path1, 'KDN_n.txt'), self.n_kdn)
            np.savetxt(path.join(path1, 'MVN_n_marginal.txt'), self.n_mvn_marginal)
            np.savetxt(path.join(path1, 'KDN_n_marginal.txt'), self.n_kdn_marginal)
            i = 0
            for parname in self.parlist:
                np.savetxt(path.join(path1, 'Par_' + parname + '_median.txt'), self.par_median[i])
                np.savetxt(path.join(path1, 'Par_' + parname + '_std.txt'), self.par_std[i])
                i += 1

        # save three-dimensional array in slices of the first parameter
        if len(self.steplist) == 3:
            for slice in range(self.results_kdn.shape[0]):
                np.savetxt(path.join(path1, 'MVN_entropy_' + str(slice) + '.txt'), self.results_mvn[slice])
                np.savetxt(path.join(path1, 'KDN_entropy_' + str(slice) + '.txt'), self.results_kdn[slice])
                np.savetxt(path.join(path1, 'MVN_entropy_marginal_' + str(slice) + '.txt'),
                           self.results_mvn_marginal[slice])
                np.savetxt(path.join(path1, 'KDN_entropy_marginal_' + str(slice) + '.txt'),
                           self.results_kdn_marginal[slice])
                np.savetxt(path.join(path1, 'MVN_infocontent_' + str(slice) + '.txt'), self.priorentropy -
                           self.results_mvn[slice])
                np.savetxt(path.join(path1, 'KDN_infocontent_' + str(slice) + '.txt'), self.priorentropy -
                           self.results_kdn[slice])
                np.savetxt(path.join(path1, 'MVN_infocontent_marginal_' + str(slice) + '.txt'),
                           self.priorentropy_marginal - self.results_mvn_marginal[slice])
                np.savetxt(path.join(path1, 'KDN_infocontent_marginal_' + str(slice) + '.txt'),
                           self.priorentropy_marginal - self.results_kdn_marginal[slice])
                np.savetxt(path.join(path1, 'MVN_sqstd_' + str(slice) + '.txt'), self.sqstd_mvn[slice])
                np.savetxt(path.join(path1, 'KDN_sqstd_' + str(slice) + '.txt'), self.sqstd_kdn[slice])
                np.savetxt(path.join(path1, 'MVN_sqstd_marginal_' + str(slice) + '.txt'),
                           self.sqstd_mvn_marginal[slice])
                np.savetxt(path.join(path1, 'KDN_sqstd_marginal_' + str(slice) + '.txt'),
                           self.sqstd_kdn_marginal[slice])
                np.savetxt(path.join(path1, 'MVN_n_' + str(slice) + '.txt'), self.n_mvn[slice])
                np.savetxt(path.join(path1, 'KDN_n_' + str(slice) + '.txt'), self.n_kdn[slice])
                np.savetxt(path.join(path1, 'MVN_n_marginal_' + str(slice) + '.txt'), self.n_mvn_marginal[slice])
                np.savetxt(path.join(path1, 'KDN_n_marginal_' + str(slice) + '.txt'), self.n_kdn_marginal[slice])

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

    def run_optimization(self):

        def writeout_result(_itindex, avg_mvn, avg_kdn, avg_mvn_marginal, avg_kdn_marginal, points_median, points_std):
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
                self.n_kdn[index] += 1.0
                self.n_mvn[index] += 1.0
                self.n_kdn_marginal[index] += 1.0
                self.n_mvn_marginal[index] += 1.0
                n = self.n_kdn[index]

                old_mvn = self.results_mvn[index]
                old_kdn = self.results_kdn[index]
                old_mvn_marginal = self.results_mvn_marginal[index]
                old_kdn_marginal = self.results_kdn_marginal[index]

                self.results_mvn[index] = running_mean(self.results_mvn[index], n, avg_mvn)
                self.results_kdn[index] = running_mean(self.results_kdn[index], n, avg_kdn)
                self.results_mvn_marginal[index] = running_mean(self.results_mvn_marginal[index], n, avg_mvn_marginal)
                self.results_kdn_marginal[index] = running_mean(self.results_kdn_marginal[index], n, avg_kdn_marginal)
                for i in range(self.par_median.shape[0]):
                    self.par_median[(i,) + index] = running_mean(self.par_median[(i,) + index], n, points_median[i])
                    # for par std the average is calculated, not a sqstd of par_median
                    self.par_std[(i,) + index] = running_mean(self.par_std[(i,) + index], n, points_std[i])

                self.sqstd_mvn[index] = running_sqstd(self.sqstd_mvn[index], n, avg_mvn, old_mvn,
                                                      self.results_mvn[index])
                self.sqstd_kdn[index] = running_sqstd(self.sqstd_kdn[index], n, avg_kdn, old_kdn,
                                                      self.results_kdn[index])
                self.sqstd_mvn_marginal[index] = running_sqstd(self.sqstd_mvn_marginal[index], n, avg_mvn_marginal,
                                                               old_mvn_marginal, self.results_mvn_marginal[index])
                self.sqstd_kdn_marginal[index] = running_sqstd(self.sqstd_kdn_marginal[index], n, avg_kdn_marginal,
                                                               old_kdn_marginal, self.results_kdn_marginal[index])

        def calc_entropy_for_index(molstat_iter, itindex):
            # Calculate Entropy n times and average
            mvn_entropy = []
            kdn_entropy = []
            mvn_entropy_marginal = []
            kdn_entropy_marginal = []
            points_median = []
            points_std = []

            for j in range(1):  # was 10
                # calculate entropy, dependent parameters == parameters of interest
                # independent parameters == nuisance parameters
                a, b, c, d, e, f = self.calc_entropy(molstat_iter)
                mvn_entropy.append(a)
                kdn_entropy.append(b)
                mvn_entropy_marginal.append(c)
                kdn_entropy_marginal.append(d)
                points_median = e
                points_std = f

            # remove outliers and average calculated entropies, don't average over parameter stats
            avg_mvn, std_mvn = average(mvn_entropy)
            avg_kdn, std_kdn = average(kdn_entropy)
            avg_mvn_marginal, std_mvn_marginal = average(mvn_entropy_marginal)
            avg_kdn_marginal, std_kdn_marginal = average(kdn_entropy_marginal)

            bValidResult = (std_kdn < self.convergence) and \
                           (self.priorentropy_marginal - avg_kdn_marginal > (-0.5) * len(
                               self.dependent_parameters)) and \
                           (self.priorentropy - avg_kdn > (-0.5) * len(self.parlist))

            # no special treatment for first entry necessary, algorithm catches this
            if bValidResult:
                writeout_result(itindex, avg_mvn, avg_kdn, avg_mvn_marginal, avg_kdn_marginal, points_median,
                                points_std)

            # save results for every iteration and delete large files
            self.save_results(self.spath)

        def set_sim_pars_for_index(it):
            qrange = 0
            pre = 0
            # cycle through all parameters
            isim = 0
            for row in self.allpar.itertuples():
                if row.par in self.steppar['par'].tolist():
                    lsim = self.steppar.loc[self.steppar['par'] == row.par, 'l_sim'].iloc[0]
                    stepsim = self.steppar.loc[self.steppar['par'] == row.par, 'step_sim'].iloc[0]
                    value = self.steppar.loc[self.steppar['par'] == row.par, 'value'].iloc[0]
                    lfit = self.steppar.loc[self.steppar['par'] == row.par, 'l_fit'].iloc[0]
                    ufit = self.steppar.loc[self.steppar['par'] == row.par, 'u_fit'].iloc[0]

                    simvalue = lsim + stepsim * it.multi_index[isim]
                    if row.type == 'fd' or row.type == 'fi':
                        # fixed fit boundaries, not floating, for such things as volfracs between 0 and 1
                        lowersim = lfit
                        uppersim = ufit
                    else:
                        lowersim = simvalue - (value - lfit)
                        uppersim = simvalue + (ufit - value)

                    # catch non-fit parameters in entropy.dat for q-range and counting time prefactor
                    if row.par == 'qrange':
                        qrange = simvalue
                    elif row.par == 'prefactor':
                        pre = simvalue
                    else:
                        self.simpar.loc[self.simpar['par'] == row.par, 'value'] = simvalue
                        self.molstat.Interactor.fnReplaceParameterLimitsInSetup(row.par, lowersim, uppersim)
                    isim += 1
                else:
                    if row.par == 'qrange':
                        qrange = row.value
                    elif row.par == 'prefactor':
                        pre = row.value
                    else:
                        self.molstat.Interactor.fnReplaceParameterLimitsInSetup(row.par, row.l_fit, row.u_fit)

            self.simpar.to_csv(path.join(self.spath, 'simpar.dat'), sep=' ', header=None, index=False)
            return pre, qrange

        def work_on_index(iteration, it, itindex):
            dirname = 'iteration_' + str(iteration)
            fulldirname = path.join(self.spath, dirname)
            path1 = path.join(fulldirname, 'save')
            chainname = path.join(path1, self.runfile+'-chain.mc')

            # fetch mode and cluster mode are exclusive
            if not self.bFetchMode:
                # run a new fit, preparations are done in the root directory and the new fit is copied into the
                # iteration directory, preparations in the iterations directory are not possible, because it would
                # be lacking a results directory, which is needed for restoring a state/parameters
                molstat_iter = molstat.CMolStat(fitsource=self.fitsource, spath=fulldirname, mcmcpath='save',
                                                runfile=self.runfile, load_state=False)
                self.molstat.Interactor.fnBackup(target=path.join(self.spath, 'simbackup'))
                pre, qrange = set_sim_pars_for_index(it)
                self.molstat.fnSimulateData(mode=self.mode, pre=pre, qrange=qrange)
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
                calc_entropy_for_index(molstat_iter, itindex)

            # delete big files except in Cluster mode. They are needed there for future fetching
            if self.deldir and not self.bClusterMode:
                rm_file(path.join(path1, self.runfile+'-point.mc'))
                rm_file(path.join(path1, self.runfile+'-chain.mc'))
                rm_file(path.join(path1, self.runfile+'-stats.mc'))
                rm_file(path.join(path1, self.runfile+'-point.mc.gz'))
                rm_file(path.join(path1, self.runfile+'-chain.mc.gz'))
                rm_file(path.join(path1, self.runfile+'-stats.mc.gz'))

        def iterate_over_all_indices(refinement=False):
            bWorkedOnIndex = False
            # the complicated iteration syntax is due the unknown dimensionality of the results space / arrays
            it = np.nditer(self.results_kdn, flags=['multi_index'])
            iteration = 0
            while not it.finished:
                itindex = tuple(it.multi_index[i] for i in range(len(self.steplist)))
                # parameter grids can be symmetric, and only one of the symmetry-related indices
                # will be calculated unless calc_symmetric is True
                if all(itindex[i] <= itindex[i + 1] for i in range(len(itindex) - 1)) or self.calc_symmetric:
                    # run MCMC if it is first time or the value in results is inf
                    invalid_result = np.isinf(self.results_kdn[itindex]) or fabs(self.results_kdn[itindex]) > 10000
                    insufficient_iterations = self.n_mvn[itindex] < self.miniter

                    outlier = False
                    if refinement:
                        # during refinement, check whether the KDN value follows that of the MVN with respect to its
                        # nearest neighbors, implemented convolution because of ill-defined origin of scipy convolute
                        conv_MVN = convolute(self.results_mvn)
                        conv_KDN = convolute(self.results_kdn)
                        dMVN = conv_MVN[itindex] - self.results_mvn[itindex]
                        dKDN = conv_KDN[itindex] - self.results_kdn[itindex]
                        if fabs(dMVN - dKDN) > self.convergence:
                            outlier = True

                    # Do we need to work on this particular index?
                    if outlier or invalid_result or insufficient_iterations or self.bFetchMode:
                        bWorkedOnIndex = True
                        work_on_index(iteration, it, itindex)

                    # iterations are only increased if this index is not dropped because of symmetry
                    iteration += 1
                it.iternext()
            return bWorkedOnIndex

        # every index has at least one result before re-analyzing any data point (refinement)
        bRefinement = False
        while True:
            if not iterate_over_all_indices(bRefinement):
                if not bRefinement:
                    # all indicees have the minimum number of iterations -> refinement starts
                    bRefinement = True
                else:
                    break
            # never repeat iterations on cluster or when just calculating entropies
            if self.bClusterMode or self.bFetchMode:
                break

        # wait for all jobs to finish
        if self.bClusterMode:
            while self.joblist:
                self.waitforjob(bFinish=True)
