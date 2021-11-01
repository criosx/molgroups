from __future__ import division

__all__ = ["entropy"]

import numpy as np
import rs
import matplotlib.pyplot as plt
import matplotlib
import pandas
import itertools
import zipfile
from numpy import mean, std, exp, log, max, sqrt, log2, pi, e
from numpy.random import permutation
from sklearn.mixture import BayesianGaussianMixture as GMM
from os import path, mkdir
from math import fabs, pow, floor, ceil
from subprocess import call, Popen
from time import sleep
from shutil import rmtree, make_archive
from sys import argv
from decimal import Decimal
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


def running_mean(current_mean, n, new_point):
    # see Tony Finch, Incremental calculation of weighted mean and variance
    return current_mean * (n - 1) / n + new_point * (1 / n)


def running_sqstd(current_sqstd, n, new_point, previous_mean, current_mean):
    # see Tony Finch, Incremental calculation of weighted mean and variance
    return (current_sqstd * (n - 1) + (new_point - previous_mean) * (new_point - current_mean)) / n


def save_plot(x, y, z, xlabel, ylabel, color, filename, zmin=None, zmax=None, levels=20):

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
    plt.savefig('./plots/' + filename + '.pdf')
    plt.close()

def cov_entropy(C):
    """
    Entropy estimate from covariance matrix C
    """
    return 0.5 * (len(C) * log2(2*pi*e) + log2(abs(np.linalg.det(C))))

class MVNEntropy(object):
    def __init__(self, x, alpha=0.05, max_points=1000):  #max points was 1000
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

    def __init__(self, mcmcburn=16000, mcmcsteps=5000, deldir=True, convergence=2.0, miniter=1, mode='water',
                 bFetchMode=False, bClusterMode=False, time=2, calc_symmetric=True, upper_info_plotlevel=None,
                 plotlimits_filename=''):

        self.mcmcburn = mcmcburn
        self.mcmcsteps = mcmcsteps
        self.deldir = deldir
        self.convergence = convergence
        self.miniter = miniter
        self.mode = mode
        self.bFetchMode = bFetchMode
        self.bClusterMode = bClusterMode
        self.time = time
        self.calc_symmetric = calc_symmetric
        self.upper_info_plotlevel = upper_info_plotlevel
        self.plotlimits_filename=plotlimits_filename

        self.ReflPar = rs.CReflectometry()

        # all parameters from entropypar.dat
        self.allpar = pandas.read_csv('entropypar.dat', sep=' ', header=None,
                                 names=['type', 'par', 'value', 'l_fit', 'u_fit', 'l_sim', 'u_sim', 'step_sim'],
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

        if path.isfile('results/MVN_entropy.npy'):
            self.results_mvn = np.load('results/MVN_entropy.npy')
            self.results_kdn = np.load('results/KDN_entropy.npy')
            self.results_mvn_marginal = np.load('results/MVN_entropy_marginal.npy')
            self.results_kdn_marginal = np.load('results/KDN_entropy_marginal.npy')
            self.n_mvn = np.load('results/MVN_n.npy')
            self.n_kdn = np.load('results/KDN_n.npy')
            self.n_mvn_marginal = np.load('results/MVN_n_marginal.npy')
            self.n_kdn_marginal = np.load('results/KDN_n_marginal.npy')
            self.sqstd_mvn = np.load('results/MVN_sqstd.npy')
            self.sqstd_kdn = np.load('results/KDN_sqstd.npy')
            self.sqstd_mvn_marginal = np.load('results/MVN_sqstd_marginal.npy')
            self.sqstd_kdn_marginal = np.load('results/KDN_sqstd_marginal.npy')
            self.par_median = np.load('results/par_median.npy')
            self.par_std = np.load('results/par_std.npy')
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

    def calc_entropy(self, dirname='MCMC'):

        N_entropy = 10000  # was 10000
        N_norm = 10000  # was 2500

        # read MCMC result and save sErr.dat
        fit_interactor = rs.CRefl1DInteractor()
        points, parnames, logp = fit_interactor.fnLoadMCMCResults(dirname=dirname)

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

    def runmcmc(self, iteration):

        # wait for a job to finish before submitting next cluster job
        if self.bClusterMode:
            self.waitforjob()

        # copy garefl/refl1d files into directory
        dirname = 'iteration_' + str(iteration)
        if not path.isdir(dirname):
            call(['mkdir', dirname])
        call(['rm', '-r', dirname + '/save'])
        call('cp *.dat ' + dirname, shell=True)
        call('cp *.cc ' + dirname, shell=True)
        call('cp *.h ' + dirname, shell=True)
        call('cp *.o ' + dirname, shell=True)
        call('cp *.py ' + dirname, shell=True)
        call('cp *.pyc ' + dirname, shell=True)
        call(['cp', 'Makefile', dirname])

        # run MCMC either cluster or local
        self.ReflPar.fnMake(dirname=dirname)
        if self.bClusterMode:
            script = []
            script.append('#!/bin/bash\n')
            script.append('#SBATCH --job-name=entro' + str(iteration) + '\n')
            script.append('#SBATCH -A mc4s9np\n')
            script.append('#SBATCH -p RM\n')
            script.append('#SBATCH -t 0' + str(self.time) + ':00:00\n')
            script.append('#SBATCH -N 4\n')
            script.append('#SBATCH --ntasks-per-node 28\n')
            script.append('\n')
            script.append('set +x\n')
            script.append('cd $SLURM_SUBMIT_DIR\n')
            # script.append('cd '+dirname+'\n')
            script.append('\n')
            script.append('module load python/2.7.11_gcc\n')
            script.append('export PYTHONPATH=/home/hoogerhe/bin/lib/python2.7/site-packages:/home/hoogerhe/src/bumps\n')
            script.append('\n')
            script.append('mpirun -np 112 python /home/hoogerhe/src/refl1d/bin/refl1d_cli.py ' + dirname +
                          '/run.py --fit=dream --mpi --init=lhs --batch --pop=28 --time=' + str(
                float(self.time) - 0.1) + ' --thin=20 --store=' + dirname +
                          '/save --burn=' + str(self.mcmcburn) + ' --steps=' + str(self.mcmcsteps) + '\n')
            # write runscript
            file = open(dirname + '/runscript', 'w')
            file.writelines(script)
            file.close()

            lCommand = ['sbatch', dirname + '/runscript']
            # For testing purposes
            # lCommand = ['refl1d_cli.py', dirname + '/run.py', '--fit=dream', '--parallel', '--init=lhs', '--batch', '--thin=200']
            # lCommand.append('--store='+dirname+'/save')
            # lCommand.append('--burn=' + str(mcmcburn))
            # lCommand.append('--steps=' + str(mcmcsteps))
            Popen(lCommand)
            self.joblist.append(iteration)

        else:
            lCommand = ['refl1d_cli.py', dirname + '/run.py', '--fit=dream', '--parallel=0', '--init=lhs', '--batch']
            lCommand.append('--store=' + dirname + '/save')
            lCommand.append('--burn=' + str(self.mcmcburn))
            lCommand.append('--steps=' + str(self.mcmcsteps))
            call(lCommand)

        return

    def plot_results(self):

        if not path.isdir('plots'):
            mkdir('plots')

        onecolormaps = [plt.cm.Greys, plt.cm.Purples, plt.cm.Blues, plt.cm.Greens, plt.cm.Oranges, plt.cm.Reds]
        ec = plt.cm.coolwarm

        if self.plotlimits_filename != '':
            plotlimits = True
            plotlim = pandas.read_csv(self.plotlimits_filename, sep=" ", header=None)
            plotlim.columns = ['name', 'value']

            if not plotlim.loc[plotlim['name'] == 'info'].empty:
                self.upper_info_plotlevel = float(plotlim.loc[plotlim['name'] == 'info']['value'])

        else:
            plotlimits = False

        if len(self.steplist) <= 2:
            # numpy array and plot axes are reversed
            ax1 = self.axes[0]
            ax0 = self.axes[1]
            sp1 = self.steppar['par'].tolist()[0]
            sp0 = self.steppar['par'].tolist()[1]

            save_plot(ax0, ax1, self.results_mvn, sp0, sp1, ec, 'MVN_entropy')
            save_plot(ax0, ax1, self.results_kdn, sp0, sp1, ec, 'KDN_entropy')
            save_plot(ax0, ax1, self.results_mvn_marginal, sp0, sp1, ec, 'MVN_entropy_marginal')
            save_plot(ax0, ax1, self.results_kdn_marginal, sp0, sp1, ec, 'KDN_entropy_marginal')
            save_plot(ax0, ax1, self.priorentropy - self.results_mvn, sp0, sp1, ec, 'MVN_infocontent', zmin=0)
            save_plot(ax0, ax1, self.priorentropy - self.results_kdn, sp0, sp1, ec, 'KDN_infocontent', zmin=0)
            save_plot(ax0, ax1, self.priorentropy_marginal - self.results_mvn_marginal, sp0, sp1, ec,
                      'MVN_infocontent_marginal', zmin=0, zmax=self.upper_info_plotlevel)
            save_plot(ax0, ax1, self.priorentropy_marginal - self.results_kdn_marginal, sp0, sp1, ec,
                      'KDN_infocontent_marginal', zmin=0, zmax=self.upper_info_plotlevel)
            save_plot(ax0, ax1, self.sqstd_mvn, sp0, sp1, ec, 'MVN_sqstd', zmin=0)
            save_plot(ax0, ax1, self.sqstd_kdn, sp0, sp1, ec, 'KDN_sqstd', zmin=0)
            save_plot(ax0, ax1, self.sqstd_mvn_marginal, sp0, sp1, ec, 'MVN_sqstd_marginal', zmin=0)
            save_plot(ax0, ax1, self.sqstd_kdn_marginal, sp0, sp1, ec, 'KDN_sqstd_marginal', zmin=0)
            save_plot(ax0, ax1, self.n_mvn, sp0, sp1, ec, 'MVN_n', zmin=0)
            save_plot(ax0, ax1, self.n_kdn, sp0, sp1, ec, 'KDN_n', zmin=0)
            save_plot(ax0, ax1, self.n_mvn_marginal, sp0, sp1, ec, 'MVN_n_marginal', zmin=0)
            save_plot(ax0, ax1, self.n_kdn_marginal, sp0, sp1, ec, 'KDN_n_marginal', zmin=0)

            i = 0
            j = 0
            for parname in self.parlist:
                if plotlimits and not plotlim.loc[plotlim['name'] == parname].empty:
                    zmax = float(plotlim.loc[plotlim['name'] == parname]['value'])
                else:
                    zmax = None
                save_plot(ax0, ax1, self.par_median[i], sp0, sp1, onecolormaps[j], 'Par_' + parname + '_median')
                save_plot(ax0, ax1, self.par_std[i], sp0, sp1, onecolormaps[j], 'Par_' + parname + '_std', zmin=0,
                          zmax=zmax)
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
                save_plot(ax1, ax2, self.results_mvn[slice], sp1, sp2, ec, 'MVN_entropy_'+str(slice))
                save_plot(ax1, ax2, self.results_kdn[slice], sp1, sp2, ec, 'KDN_entropy_'+str(slice))
                save_plot(ax1, ax2, self.results_mvn_marginal[slice], sp1, sp2, ec, 'MVN_entropy_marginal_'+str(slice))
                save_plot(ax1, ax2, self.results_mvn_marginal[slice], sp1, sp2, ec, 'MVN_entropy_marginal_'+str(slice))
                save_plot(ax1, ax2, self.priorentropy - self.results_mvn[slice], sp1, sp2, ec,
                          'MVN_infocontent_'+str(slice), zmin=0)
                save_plot(ax1, ax2, self.priorentropy - self.results_kdn[slice], sp1, sp2, ec,
                          'KDN_infocontent_'+str(slice), zmin=0)
                save_plot(ax1, ax2, self.priorentropy_marginal - self.results_mvn_marginal[slice], sp1, sp2,
                          ec, 'MVN_infocontent_marginal_'+str(slice), zmin=0)
                save_plot(ax1, ax2, self.priorentropy_marginal - self.results_kdn_marginal[slice], sp1, sp2,
                          ec, 'KDN_infocontent_marginal_'+str(slice), zmin=0)
                save_plot(ax1, ax2, self.sqstd_mvn[slice], sp1, sp2, ec, 'MVN_sqstd_'+str(slice), zmin=0)
                save_plot(ax1, ax2, self.sqstd_kdn[slice], sp1, sp2, ec, 'KDN_sqstd_'+str(slice), zmin=0)
                save_plot(ax1, ax2, self.sqstd_mvn_marginal[slice], sp1, sp2, ec, 'MVN_sqstd_marginal_'+str(slice),
                          zmin=0)
                save_plot(ax1, ax2, self.sqstd_kdn_marginal[slice], sp1, sp2, ec, 'KDN_sqstd_marginal_'+str(slice),
                          zmin=0)
                save_plot(ax1, ax2, self.n_mvn[slice], sp1, sp2, ec, 'MVN_n_'+str(slice), zmin=0)
                save_plot(ax1, ax2, self.n_kdn[slice], sp1, sp2, ec, 'KDN_n_'+str(slice), zmin=0)
                save_plot(ax1, ax2, self.n_mvn_marginal[slice], sp1, sp2, ec, 'MVN_n_marginal_'+str(slice), zmin=0)
                save_plot(ax1, ax2, self.n_kdn_marginal[slice], sp1, sp2, ec, 'KDN_n_marginal_'+str(slice), zmin=0)


    def save_results(self):

        if not path.isdir('results'):
            mkdir('results')

        np.save('results/KDN_entropy', self.results_kdn, allow_pickle=False)
        np.save('results/MVN_entropy', self.results_mvn, allow_pickle=False)
        np.save('results/KDN_entropy_marginal', self.results_kdn_marginal, allow_pickle=False)
        np.save('results/MVN_entropy_marginal', self.results_mvn_marginal, allow_pickle=False)
        np.save('results/KDN_infocontent', self.priorentropy - self.results_kdn, allow_pickle=False)
        np.save('results/MVN_infocontent', self.priorentropy - self.results_mvn, allow_pickle=False)
        np.save('results/KDN_infocontent_marginal', self.priorentropy_marginal - self.results_kdn_marginal,
                allow_pickle=False)
        np.save('results/MVN_infocontent_marginal', self.priorentropy_marginal - self.results_mvn_marginal,
                allow_pickle=False)
        np.save('results/KDN_sqstd', self.sqstd_kdn, allow_pickle=False)
        np.save('results/MVN_sqstd', self.sqstd_mvn, allow_pickle=False)
        np.save('results/KDN_sqstd_marginal', self.sqstd_kdn_marginal, allow_pickle=False)
        np.save('results/MVN_sqstd_marginal', self.sqstd_mvn_marginal, allow_pickle=False)
        np.save('results/KDN_n', self.n_kdn, allow_pickle=False)
        np.save('results/MVN_n', self.n_mvn, allow_pickle=False)
        np.save('results/KDN_n_marginal', self.n_kdn_marginal, allow_pickle=False)
        np.save('results/MVN_n_marginal', self.n_mvn_marginal, allow_pickle=False)
        np.save('results/par_median', self.par_median, allow_pickle=False)
        np.save('results/par_std', self.par_std, allow_pickle=False)

        # save to txt when not more than two-dimensional array
        if len(self.steplist) <= 2:
            np.savetxt('results/MVN_entropy.txt', self.results_mvn)
            np.savetxt('results/KDN_entropy.txt', self.results_kdn)
            np.savetxt('results/MVN_entropy_marginal.txt', self.results_mvn_marginal)
            np.savetxt('results/KDN_entropy_marginal.txt', self.results_kdn_marginal)
            np.savetxt('results/MVN_infocontent.txt', self.priorentropy - self.results_mvn)
            np.savetxt('results/KDN_infocontent.txt', self.priorentropy - self.results_kdn)
            np.savetxt('results/MVN_infocontent_marginal.txt', self.priorentropy_marginal - self.results_mvn_marginal)
            np.savetxt('results/KDN_infocontent_marginal.txt', self.priorentropy_marginal - self.results_kdn_marginal)
            np.savetxt('results/MVN_sqstd.txt', self.sqstd_mvn)
            np.savetxt('results/KDN_sqstd.txt', self.sqstd_kdn)
            np.savetxt('results/MVN_sqstd_marginal.txt', self.sqstd_mvn_marginal)
            np.savetxt('results/KDN_sqstd_marginal.txt', self.sqstd_kdn_marginal)
            np.savetxt('results/MVN_n.txt', self.n_mvn)
            np.savetxt('results/KDN_n.txt', self.n_kdn)
            np.savetxt('results/MVN_n_marginal.txt', self.n_mvn_marginal)
            np.savetxt('results/KDN_n_marginal.txt', self.n_kdn_marginal)
            i = 0
            for parname in self.parlist:
                np.savetxt('results/Par_' + parname + '_median.txt', self.par_median[i])
                np.savetxt('results/Par_' + parname + '_std.txt', self.par_std[i])
                i += 1

        # save three-dimensional array in slices of the first parameter
        if len(self.steplist) == 3:
            for slice in range(self.results_kdn.shape[0]):
                np.savetxt('results/MVN_entropy_' + str(slice) + '.txt', self.results_mvn[slice])
                np.savetxt('results/KDN_entropy_' + str(slice) + '.txt', self.results_kdn[slice])
                np.savetxt('results/MVN_entropy_marginal_' + str(slice) + '.txt', self.results_mvn_marginal[slice])
                np.savetxt('results/KDN_entropy_marginal_' + str(slice) + '.txt', self.results_kdn_marginal[slice])
                np.savetxt('results/MVN_infocontent_' + str(slice) + '.txt', self.priorentropy - self.results_mvn[slice])
                np.savetxt('results/KDN_infocontent_' + str(slice) + '.txt', self.priorentropy - self.results_kdn[slice])
                np.savetxt('results/MVN_infocontent_marginal_' + str(slice) + '.txt', self.priorentropy_marginal
                              - self.results_mvn_marginal[slice])
                np.savetxt('results/KDN_infocontent_marginal_' + str(slice) + '.txt', self.priorentropy_marginal
                              - self.results_kdn_marginal[slice])
                np.savetxt('results/MVN_sqstd_' + str(slice) + '.txt', self.sqstd_mvn[slice])
                np.savetxt('results/KDN_sqstd_' + str(slice) + '.txt', self.sqstd_kdn[slice])
                np.savetxt('results/MVN_sqstd_marginal_' + str(slice) + '.txt', self.sqstd_mvn_marginal[slice])
                np.savetxt('results/KDN_sqstd_marginal_' + str(slice) + '.txt', self.sqstd_kdn_marginal[slice])
                np.savetxt('results/MVN_n_' + str(slice) + '.txt', self.n_mvn[slice])
                np.savetxt('results/KDN_n_' + str(slice) + '.txt', self.n_kdn[slice])
                np.savetxt('results/MVN_n_marginal_' + str(slice) + '.txt', self.n_mvn_marginal[slice])
                np.savetxt('results/KDN_n_marginal_' + str(slice) + '.txt', self.n_kdn_marginal[slice])

    def waitforjob(self, bFinish=False):
        # finish flag means that parent is waiting for all jobs to finish and not because of a too long
        # job queue
        if len(self.joblist) >= self.jobmax or bFinish:
            repeat = True
            while repeat:
                for job in self.joblist:
                    if path.isfile('iteration_' + str(job) + '/save/run-chain.mc'):
                        # wait 2 minutes to allow all output files to be written
                        sleep(180)
                        # zip up finished job
                        make_archive('iteration_' + str(job), 'zip', 'iteration_' + str(job))
                        call(['rm', '-r', 'iteration_' + str(job)])
                        self.joblist.remove(job)
                        repeat = False
                        break
                sleep(60)
        return

    def writeoutresult(self, _itindex, avg_mvn, avg_kdn, avg_mvn_marginal, avg_kdn_marginal,
                       points_median, points_std):
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

            self.sqstd_mvn[index] = running_sqstd(self.sqstd_mvn[index], n, avg_mvn, old_mvn, self.results_mvn[index])
            self.sqstd_kdn[index] = running_sqstd(self.sqstd_kdn[index], n, avg_kdn, old_kdn, self.results_kdn[index])
            self.sqstd_mvn_marginal[index] = running_sqstd(self.sqstd_mvn_marginal[index], n, avg_mvn_marginal,
                                                           old_mvn_marginal, self.results_mvn_marginal[index])
            self.sqstd_kdn_marginal[index] = running_sqstd(self.sqstd_kdn_marginal[index], n, avg_kdn_marginal,
                                                           old_kdn_marginal, self.results_kdn_marginal[index])

    def run_optimization(self):

        repeats = True
        # indicator whether every systematic variation has at least one result
        # this should be achieved before re-analyzing any data point
        no_zeros = False
        while repeats:

            repeats = False
            it = np.nditer(self.results_kdn, flags=['multi_index'])
            iteration = 0

            while not it.finished:

                itindex = tuple(it.multi_index[i] for i in range(len(self.steplist)))
                donotdropindex = self.calc_symmetric or all(itindex[i] <= itindex[i + 1]
                                                            for i in range(len(itindex) - 1))

                if donotdropindex:

                    # run MCMC if it is first time or the value in results is inf
                    invalid_result = np.isinf(self.results_kdn[itindex]) or fabs(self.results_kdn[itindex]) > 10000
                    insufficient_iterations = self.n_mvn[itindex] < self.miniter
                    outlier = False

                    if no_zeros and (not insufficient_iterations):
                        # if a valid result exists, check whether the KDN value follows that of the MVN
                        # with respect to its nearest neighbors

                        # implemented own convolution because of ill-defined origin of scipy convolute
                        conv_MVN = convolute(self.results_mvn)
                        conv_KDN = convolute(self.results_kdn)

                        dMVN = conv_MVN[itindex] - self.results_mvn[itindex]
                        dKDN = conv_KDN[itindex] - self.results_kdn[itindex]

                        if fabs(dMVN - dKDN) > self.convergence:
                            outlier = True

                    if outlier or invalid_result or insufficient_iterations or (self.n_mvn[itindex] == 0) \
                            or self.bFetchMode:
                        repeats = True
                        # set up fit limits and simulation parameters including parameters that are varied and
                        # fixed during the process
                        isim = 0
                        qrange = 0
                        pre = 0
                        bUnzippedPriorResult = path.isfile('iteration_' + str(iteration) + '/save/run-chain.mc')
                        bPriorResultExists = bUnzippedPriorResult or path.isfile('iteration_' + str(iteration) + '.zip')
                        for row in self.allpar.itertuples():  # cycle through all parameters
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
                                    if not bPriorResultExists or not self.bFetchMode:
                                        self.ReflPar.fnReplaceParameterLimitsInSetup(row.par, lowersim, uppersim)
                                isim += 1
                            else:
                                if row.par == 'qrange':
                                    qrange = row.value
                                elif row.par == 'prefactor':
                                    pre = row.value
                                else:
                                    if not bPriorResultExists or not self.bFetchMode:
                                        self.ReflPar.fnReplaceParameterLimitsInSetup(row.par, row.l_fit, row.u_fit)

                        if ((not bPriorResultExists and not self.bFetchMode) or
                                (self.bClusterMode and bUnzippedPriorResult)):
                            # fetch mode and cluster mode are exclusive
                            # there should be no unzipped results in cluster mode, as cluster mode performs only one
                            # single iteration over the parameter space, because the entropy has to be calculated
                            # offsite. therefore, if an unzipped result is found, it is ignored
                            self.simpar.to_csv('simpar.dat', sep=' ', header=None, index=False)
                            self.ReflPar.fnSimulateReflectivity(mode=self.mode, pre=pre, qrange=qrange)
                            self.runmcmc(iteration)
                            bPriorResultExists = True

                        # Entropy calculation
                        # do not run entropy calculation when on cluster
                        if not self.bClusterMode and bPriorResultExists:

                            # check if result directory is zipped. If yes, unzip.
                            if path.isfile('iteration_' + str(iteration) + '.zip'):
                                if not path.isfile('iteration_' + str(iteration) + '/save/run-chain.mc'):
                                    # if a zip file and unzipped file exists -> prefer the unzipped file
                                    File = zipfile.ZipFile('iteration_' + str(iteration) + '.zip', 'r')
                                    call(['mkdir', 'iteration_' + str(iteration)])
                                    File.extractall('iteration_' + str(iteration))
                                    File.close()
                                call(['rm', 'iteration_' + str(iteration) + '.zip'])

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
                                a, b, c, d, e, f = self.calc_entropy(dirname='iteration_' + str(iteration) + '/save')
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
                                self.writeoutresult(itindex, avg_mvn, avg_kdn, avg_mvn_marginal, avg_kdn_marginal,
                                                    points_median, points_std)

                            # save results for every iteration and delete large files
                            self.save_results()
                            if self.deldir:
                                call(['rm', 'iteration_' + str(iteration) + '/save/run-point.mc'])
                                call(['rm', 'iteration_' + str(iteration) + '/save/run-chain.mc'])
                                call(['rm', 'iteration_' + str(iteration) + '/save/run-stats.mc'])

                    iteration += 1
                    it.iternext()

            # Repeats is False while no_zeros is False means that the algorithm went over all iterations but has not
            # found one with an insufficient number of results or invalid results. This means now, one can check for
            # outliers and improve the statistics.
            if repeats is False and no_zeros is False:
                repeats = True
                no_zeros = True

            # Never repeat iterations in cluster or when just calculating entropies
            if self.bClusterMode or self.bFetchMode:
                repeats = False

        # wait for all jobs to finish
        if self.bClusterMode:
            while self.joblist:
                self.waitforjob(bFinish=True)


if __name__ == "__main__":
    i = 1
    burn = 16000
    steps = 5000
    convergence = 2.0
    miniter = 1
    mode = 'water'
    bClusterMode = False
    bFetchMode = False
    time = 2
    bcalcsymmetric = True
    plotonly = False
    upper_info_plotlevel = None
    filename = ''
    calcsingle = False

    while i < len(argv):
        if argv[i] == '-burn':
            burn = int(argv[i + 1])
            i += 2
        elif argv[i] == '-steps':
            steps = int(argv[i + 1])
            i += 2
        elif argv[i] == '-convergence':
            convergence = int(argv[i + 1])
            i += 2
        elif argv[i] == '-miniter':
            miniter = int(argv[i + 1])
            i += 2
        elif argv[i] == '-mode':
            mode = argv[i + 1]
            i += 2
        elif argv[i] == '-cluster':
            bClusterMode = True
            i += 1
        elif argv[i] == '-fetch':
            bFetchMode = True
            i += 1
        elif argv[i] == '-time':
            time = int(argv[i + 1])
            i += 2
        elif argv[i] == '-symmetry':
            bcalcsymmetric = False
            i += 1
        elif argv[i] == '-plotonly':
            plotonly = True
            i += 1
        elif argv[i] == '-upperinfo':
            upper_info_plotlevel = float(argv[i+1])
            i += 2
        elif argv[i] == '-plotlimits':
            filename = argv[i+1]
            i += 2
        elif argv[i] == '-calcsingle':
            calcsingle = True
            filename = argv[i+1]
            i += 2



    entropy = Entropy(mcmcburn=burn, mcmcsteps=steps, convergence=convergence, miniter=miniter, mode=mode,
                      bClusterMode=bClusterMode, bFetchMode=bFetchMode, time=time, calc_symmetric=bcalcsymmetric,
                      upper_info_plotlevel=upper_info_plotlevel, plotlimits_filename=filename)

    if plotonly:
        entropy.plot_results()
    elif calcsingle:
        a, b, c, d, e, f = entropy.calc_entropy(dirname=filename)
        print('MVN entropy '+str(a)+'\n')
        print('KDN entropy '+str(b)+'\n')
        print('MVN entropy marginal '+str(c)+'\n')
        print('KDN entropy marginal '+str(d)+'\n')

    else:
        entropy.run_optimization()
