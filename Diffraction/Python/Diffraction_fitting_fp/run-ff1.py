from bumps.names import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import scipy.fft
import sys
sys.path.append(".")
import molgroups as mol

# Define model function
def modelformfactor(lq, l_lipid, vf_bilayer, sigma, bulknsld, prefactor, dq, rel_pos, methyl_sigma):
    #TODO: Think about how to set those variables more conveniently
    maxarea = 100
    stepsize = 0.5
    dimension = 300
    startz = 50

    aArea = np.zeros(dimension).tolist()
    anSL  = np.zeros(dimension).tolist()
    anSLD = np.zeros(dimension).tolist()

    bilayer.headgroup1.fnSet(lh1, rel_pos)
    # bilayer.headgroup2.fnSet(lh1, rel_pos) 
    bilayer.fnSet(sigma, bulknsld, startz, l_lipid, l_lipid, vf_bilayer)
    dMaxArea, aArea, anSL = bilayer.fnWriteProfile(aArea, anSL, dimension, stepsize, maxarea)

    #TODO: speedup
    for i in range(len(aArea)):
        if aArea[i] != 0:
            anSLD[i] = anSL[i] / (aArea[i]*stepsize) * aArea[i]/dMaxArea + bulknsld * (1 - aArea[i]/dMaxArea)
        else:
            anSLD[i] = bulknsld

    center = bilayer.fnGetCenter()
    center = center//stepsize
    canvas_center = dimension//2
    n = int(canvas_center - center)
    centered_bilayer = np.roll(anSLD, n)
    symmetrized_bilayer = np.add(centered_bilayer,centered_bilayer[::-1])*0.5
    symmetrized_bilayer -= bulknsld
    half_bilayer = symmetrized_bilayer[int(dimension/2):]

    #TODO: Make sure that lq and x are roughly comparable
    dct_dimension = 5000
    F = scipy.fft.dct(half_bilayer, n=dct_dimension)
    F = np.abs(F)
    x = np.array([np.pi/(2*dct_dimension*stepsize)*(2*i)+dq  for i in range(int(dct_dimension))])

    #interpolate (x, F) onto lq -> (lq, modelform)
    modelform =  np.interp(lq, x, F, left=None, right=None, period=None)*prefactor

    return modelform

# Load experimental data

F2 = np.loadtxt("Experimental_form_factors/dopc.dat", skiprows=1)
F2 = np.abs(F2)
form_exp = F2[:,1]
dform_exp = np.zeros(len(form_exp))
#constant error bar estimate of .05 Å
for i in range(len(form_exp)):
    dform_exp[i] = 0.05
q_exp = F2[:,0]

# Initialize bilayer model
#-----------------------------------------------------------------------------------------
bilayer = mol.BLM_quaternary()
na1, nh1, nm1, va1, vm1, vh1, lh1 = 0.00760, 0.00461, 0.000468, 972.00, 98, 331.00, 9.56
na2, nh2, nm2, va2, vm2, vh2, lh2 = 0, 0, 0, 0, 0, 0, 0
na3, nh3, nm3, va3, vm3, vh3, lh3 = 0, 0, 0, 0, 0, 0, 0
vc, nc = 0, 0
bilayer.fnInit(va1, na1, vm1, nm1, vh1,nh1, lh1, va2, na2, vm2, nm2, vh2, nh2, lh2, va3, na3, vm3, nm3, vh3, nh3,lh3, vc, nc)


# Define variables
#----------------------------------------------------------------------------------------
l_lipid = 11.6
vf_bilayer = 1.0
sigma = 2.0
bulknsld = 9.4114E-06
prefactor = 15000
dq = 0.
rel_pos = .5
methyl_sigma = 2
#----------------------------------------------------------------------------------------

M1 = Curve(modelformfactor, q_exp, form_exp, dform_exp, l_lipid=l_lipid, vf_bilayer=vf_bilayer, sigma=sigma, bulknsld=bulknsld, prefactor=prefactor, dq=dq, rel_pos=rel_pos, methyl_sigma= methyl_sigma)
M1.l_lipid.range(9,13)
# M1.vf_bilayer.range(0.95,1.0)
M1.sigma.range(1.0, 4.0)
M1.bulknsld.range(9e-6,10e-6)
M1.prefactor.range(10000,20000)
M1.dq.range(-0.01, 0.01)
# M1.rel_pos.range(0, 1)

model = M1
problem = FitProblem(model)
