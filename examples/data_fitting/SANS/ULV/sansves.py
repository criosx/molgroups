
from bumps.names import *
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasmodels.data import load_data, plot_data

from molgroups import mol
from molgroups import components as cmp
from molgroups import lipids

import numpy

# constants
dimension = 100
stepsize = 1
z = numpy.linspace(0, 99, 100, endpoint=True)

def bilayer(z, bulknsld, l_lipid1, l_lipid2, sigma=3.0, vf_bilayer=1.0):

    # Scale all SLDs from Refl1D units (1e-6 Ang^-2) to molgroups units (Ang^-2)
    bulknsld *= 1e-6

    blm.fnSet(sigma=sigma, bulknsld=bulknsld, startz=0, l_lipid1=l_lipid1, l_lipid2=l_lipid2, vf_bilayer=vf_bilayer)

    print(z)

    # Calculate scattering properties of volume occupied by bilayer
    normarea, area, nsl = blm.fnWriteProfile(z)

    # Fill in the remaining volume with buffer of appropriate nSLD
    nsld = nsl / (normarea * np.gradient(z)) + (1.0 - area / normarea) * bulknsld

    # export objects for post analysis, needs to be from this function
    #problem.bilayers = [blm]
    #problem.dimension = dimension
    #problem.stepsize = stepsize
    #problem.moldat = blm.fnWritePar2Dict({}, 'bilayer', np.arange(dimension) * stepsize)

    # Return nSLD profile in Refl1D units
    return nsld*1e6

### Define bilayer object
DOPC = cmp.Lipid(name='DOPC', headgroup=lipids.PC, tails=[cmp.oleoyl, cmp.oleoyl], methyls=cmp.methyl)
blm = mol.BLM(lipids=[DOPC], lipid_nf=[1.0])


# IMPORT THE DATA USED
data0 = load_data('sim0.dat')
data1 = load_data('sim1.dat')
#data2 = load_data('sim2.dat')

#setattr(data, 'qmin', 0.01)
#setattr(data, 'qmax', 10.0)

# DEFINE THE MODEL
# copied and initialized the custom model with up to 100 sld and thickness parameters in sasmodels/models
# a second hard-coded limit in sasmodels/modelinfo.py line 594 was manually increased from 20 to 120.
kernel = load_model('my_core_multi_shell')

pars0 = dict(scale=0.002, background=0.15, sld_core=3.0, sld_solvent=2.4, radius=60.0, radius_pd=0.3, n=100)
pars1 = dict(scale=0.002, background=0.15, sld_core=3.0, sld_solvent=3.4, radius=60.0, radius_pd=0.3, n=100)
#pars2 = dict(scale=0.002, background=0.15, sld=3.0, sld_solvent=4.4, radius=60.0, radius_pd=0.3)

model0 = Model(kernel, **pars0)
model1 = Model(kernel, **pars1)
#model2 = Model(kernel, **pars2)

# PARAMETER RANGES (ONLY THOSE PARAMETERS ARE FITTED)

background0 = Parameter(name='background0', value=0.4).range(0.01, 0.8)
background1 = Parameter(name='background1', value=0.4).range(0.01, 0.8)
#background2 = Parameter(name='background2', value=0.4).range(0.01, 0.8)
sld_solvent0 = sld_core0 =Parameter(name='sld_solvent0', value=2.0).range(-0.56, 6.4)
sld_solvent1 = sld_core1 = Parameter(name='sld_solvent1', value=2.0).range(-0.56, 6.4)
#sld_solvent2 = Parameter(name='sld_solvent2', value=2.0).range(-0.56, 6.4)

# bilayer specific parameters
l_lipid = Parameter(name='l_lipid', value=13).range(10., 16.)

model0.scale.range(0.0001, 0.05)
model0.background = background0
model0.sld_solvent = sld_solvent0
model0.sld_core = sld_core0
model0.radius.range(40., 120.)
model0.radius_pd.range(0.05, 0.7)

model1.scale = model0.scale
model1.background = background1
model1.sld_solvent = sld_solvent1
model1.sld_core = sld_core1
model1.radius = model0.radius
model1.radius_pd=model0.radius_pd

# bilayer update
# have to figure out how to do a dynamic update
# for the moment, that's static here
sldarr0 = bilayer(z, model0.sld_solvent, l_lipid, l_lipid)
sldarr1 = bilayer(z, model1.sld_solvent, l_lipid, l_lipid)
for i in range(pars0['n']):
    s0 = getattr(model0, 'sld'+str(i+1))
    s0 = sldarr0[i]
    s1 = getattr(model1, 'sld'+str(i+1))
    s1 = sldarr1[i]
    t0 = getattr(model0, 'thickness'+str(i+1))
    t0 = 1.
    t1 = getattr(model1, 'thickness'+str(i+1))
    t1 = 1.



M0 = Experiment(data=data0, model=model0)
M1 = Experiment(data=data1, model=model1)
#M2 = Experiment(data=data2, model=model2)

problem = MultiFitProblem([M0, M1])
#problem = MultiFitProblem([M0, M1, M2])
