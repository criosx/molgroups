import molgroups as mol
from bumps.names import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import scipy.fft as fft
import sys
sys.path.append(".")
# from bumps.fitters import fit
# from bumps.formatnum import format_uncertainty
# from bumps.mapper import MPMapper
# from bumps import varplot
from bumps.stats import var_stats, format_vars, save_vars
# from pylab import figure, savefig, suptitle, rcParams

# Define model function


def modelformfactor(lq, l_lipid, sigma, bulknsld, prefactor, dq, rel_pos, hg_thickness, methyl_sigma):
    # TODO: Think about how to set those variables more conveniently
    maxarea = 100
    stepsize = 0.5
    dimension = 300
    startz = 50

    aArea = np.zeros(dimension).tolist()
    aSL = np.zeros(dimension).tolist()
    
    bilayer.headgroup1.fnSet(hg_thickness, rel_pos)
    bilayer.headgroup2.fnSet(hg_thickness, rel_pos)
    bilayer.methyl_sigma = methyl_sigma
    bilayer.fnSet(sigma, bulknsld, startz, l_lipid, l_lipid, vf_bilayer)
    normArea, aArea, aSL = bilayer.fnWriteProfile(aArea, aSL, dimension, stepsize, maxarea)
    aSLD = getSLD(aArea, aSL, dimension, stepsize, normArea, bulknsld)
    modelform = computeFormFactor(lq, aSLD, dimension, stepsize, bulknsld, prefactor, dq)

    return modelform

def getSLD(aArea, aSL, dimension, stepsize, normArea, bulknsld):
    aSLD = np.zeros(dimension).tolist()
    for i in range(dimension):
        if aArea[i] != 0:
            aSLD[i] = aSL[i] / (aArea[i]*stepsize) * aArea[i] / \
                normArea + bulknsld * (1 - aArea[i]/normArea)
        else:
            aSLD[i] = bulknsld
    return aSLD

def computeFormFactor(lq, aSLD, dimension, stepsize, bulknsld, prefactor, dq):
    center = bilayer.fnGetCenter()//stepsize
    canvas_center = dimension//2
    n = int(canvas_center - center)
    centered_bilayer = np.roll(aSLD, n)
    symmetrized_bilayer = np.add(centered_bilayer, centered_bilayer[::-1])*0.5
    symmetrized_bilayer -= bulknsld
    half_bilayer = symmetrized_bilayer[int(dimension/2):]

    # TODO: Make sure that lq and x are roughly comparable
    dct_dimension = 5000
    F = fft.dct(half_bilayer, n=dct_dimension)
    F = np.abs(F)
    x = np.array([np.pi/(2*dct_dimension*stepsize)*(2*i) +
                 dq for i in range(int(dct_dimension))])

    # interpolate (x, F) onto lq -> (lq, modelform)
    return np.interp(lq, x, F, left=None, right=None, period=None)*prefactor

# Load experimental data
F2 = np.loadtxt("Experimental_form_factors/dopc.dat", skiprows=1)
F2 = np.abs(F2)
q_exp = F2[:,0]
form_exp = F2[:,1]
#constant error bar estimate of .05 Å
dform_exp = [0.05]*len(form_exp)

# Initialize bilayer model
# -----------------------------------------------------------------------------------------
bilayer = mol.BLM_quaternary()
na1, nh1, nm1, va1, vm1, vh1, lh1 = 7.5978E-03, 4.6150E-03, 5.0652E-04, 972.00, 98, 331.00, 9.56
na2, nh2, nm2, va2, vm2, vh2, lh2 = 0, 0, 0, 0, 0, 0, 0
na3, nh3, nm3, va3, vm3, vh3, lh3 = 0, 0, 0, 0, 0, 0, 0
vc, nc = 0, 0
bilayer.fnInit(va1, na1, vm1, nm1, vh1, nh1, lh1, va2, na2, vm2,
               nm2, vh2, nh2, lh2, va3, na3, vm3, nm3, vh3, nh3, lh3, vc, nc)


# Define variables
# ----------------------------------------------------------------------------------------
l_lipid = 11.6
vf_bilayer = 1.0
sigma = 2.0
bulknsld = 9.4114E-06
prefactor = 15000
dq = 0.
rel_pos = .5
methyl_sigma = 2
# ----------------------------------------------------------------------------------------

M1 = Curve(modelformfactor, q_exp, form_exp, dform_exp, l_lipid=l_lipid, sigma=sigma, bulknsld=bulknsld,
           prefactor=prefactor, dq=dq, rel_pos=rel_pos, hg_thickness=lh1, methyl_sigma=methyl_sigma)
M1.l_lipid.range(9, 13)
M1.sigma.range(1.0, 4.0)
M1.bulknsld.range(9e-6, 10e-6)
M1.prefactor.range(5000, 30000)
M1.dq.range(-0.02, 0.02)
M1.hg_thickness.range(8, 14)
M1.rel_pos.range(0, 1)
M1.methyl_sigma.range(0, 4)

model = M1
problem = FitProblem(model)
# # mapper = MPMapper.start_mapper(problem, None, cpus=0)
result = fit(problem, method='dream', samples=10, burn=10, steps=10, thin=1, alpha=0, outliers='none', trim = 'none')
# draw = result.state.draw(portion=1)
# all_vstats = var_stats(draw)
# figure(figsize=varplot.var_plot_size(len(all_vstats)))
# varplot.plot_vars(draw, all_vstats, nbins=nbins)
# print("final chisq", problem.chisq_str())
# for k, v, dv in zip(problem.labels(), result.x, result.dx):
#     print(k, ":", format_uncertainty(v, dv))
