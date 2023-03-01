## === Import section ===

import sys
import numpy as np
from molgroups import mol
from molgroups import components as cmp
from molgroups import lipids
from refl1d.names import load4, SLD, Slab, Experiment, FitProblem, Parameter
from refl1d.flayer import FunctionalProfile
import time

## === Film structure definition section ===

"""
def make_samples(func, substrate, contrasts, **kwargs):

    samples = []

    for contrast in contrasts:
        mollayer = FunctionalProfile(dimension*stepsize, 0, profile=func, bulknsld=contrast.rho, **kwargs)
        layer_contrast = Slab(material=contrast, thickness=0.0000, interface=5.0000)
        samples.append(substrate | mollayer | layer_contrast)

    return samples
"""

### Define bilayer parameters
vf_bilayer = Parameter(name='volume fraction bilayer', value=0.9).range(0.0, 1.0)
l_lipid1 = Parameter(name='inner acyl chain thickness', value=10.0).range(8, 16)
l_lipid2 = Parameter(name='outer acyl chain thickness', value=10.0).range(8, 16)
sigma = Parameter(name='bilayer roughness', value=5).range(2, 9)
global_rough = Parameter(name ='tiox roughness', value=5).range(2, 9)
l_tiox = Parameter(name='total tiox thickness', value=120).range(50, 150)
l_submembrane = Parameter(name='submembrane thickness', value=10).range(0, 50)
dopc_nf = Parameter(name='dopc nf', value = 0.8).range(0., 1.)
do_volume = Parameter(name='dioleoyl volume', value=800).range(500, 1000)

### Define molgroups space.
dimension=300       # Number of steps
stepsize=0.5        # Length of steps

## === Stack ===
##
## First, we create a 'material' for each bulk layer, which has an real and imaginary
## scattering length density, stored in a Refl1d object called 'SLD'
d2o = SLD(name='d2o', rho=6.3000, irho=0.0000)
h2o = SLD(name='h2o', rho=-0.56, irho=0.0000)
tiox = SLD(name='tiox', rho=2.1630, irho=0.0000)
siox = SLD(name='siox', rho=4.1000, irho=0.0000)
silicon = SLD(name='silicon', rho=2.0690, irho=0.0000)

## Then bulk layers are created, each with its own 'material'.  If you want to force
## two layers to always match SLD you can use the same material in multiple layers.
## The roughnesses of each layer are set to zero to begin with:

layer_siox = Slab(material=siox, thickness=7.5804, interface=5.000)
layer_silicon = Slab(material=silicon, thickness=0.0000, interface=0.0000)

## Use the bilayer definition function to generate the bilayer SLD profile, passing in the relevant parameters.
## Note that substrate and bulk SLDs are linked to their respective materials.

## Stack the layers into individual samples, using common layer objects for layers that are unchanged between samples
## As a convention, always build the sample from the substrate up. If the neutron beam is incident from the substrate side,
## set back_reflectivity = True in the probe definition later.

### Define bilayer object
DOPC = cmp.Lipid(name='DOPC', headgroup=lipids.PC, tails=[cmp.oleoyl, cmp.oleoyl], methyls=cmp.methyl)

def make_samples(substrate, contrasts):

    samples = []

    for contrast in contrasts:
        blm = mol.ssBLM(lipids=[DOPC, lipids.DOPS], lipid_nf=[0.5, 0.5], name='bilayer')
        blm.l_siox = 0.0
        blm.rho_siox = 0.0
        
        layer_tiox = Slab(material=tiox, thickness=l_tiox - (blm.substrate.z + 0.5 * blm.substrate.length), interface=0.0)

        link_parameters = [(blm, 'global_rough', global_rough),
                           (blm, 'rho_substrate', tiox.rho * 1e-6),
                           (blm, 'l_submembrane', sigma),
                           (blm, 'bulknsld', contrast.rho*1e-6),
                           (blm, 'l_lipid1', l_lipid1),
                           (blm, 'l_lipid2', l_lipid2),
                           (blm, 'vf_bilayer', vf_bilayer),
                           (blm, 'inner_lipid_nf', [dopc_nf, 1 - dopc_nf]),
                           (blm, 'outer_lipid_nf', [dopc_nf, 1 - dopc_nf]),
                           (blm.methylene1_1, 'vol', do_volume)]

        mollayer = mol.MolLayer(contrast.rho*1e-6, link_parameters, base_groups=[blm], thickness=dimension*stepsize, interface=0, name=f'{contrast.name} mollayer')

        layer_contrast = Slab(material=contrast, thickness=0.0000, interface=5.0000)

        samples.append(substrate | layer_tiox | mollayer | layer_contrast)

    return samples

substrate = layer_silicon | layer_siox
sample, sampleh = make_samples(substrate, [d2o, h2o])

#sampleh = layer_silicon | layer_siox | layer_tiox | mollayerh | layer_h2o

## Set sample parameter ranges and constraints between layer properties

# nSLD parameters
d2o.rho.range(5.3000, 6.5000)
h2o.rho.range(-0.6, 0.6)
tiox.rho.range(1.1630, 3.1630)
siox.rho.range(3.1000, 5.1000)

# layer thickness parameters
#layer_tiox.thickness.range(66.379, 266.38)
layer_siox.thickness.range(5, 40)

# layer roughness parameters
###################################################################
## the 'interface' associated with layer0 is the boundary between #
## layer0 and layer1, and similarly for layer(N) and layer(N+1)   #
###################################################################
layer_siox.interface.range(2.0000, 9.000)

# Si and SiOx roughnesses are the same
layer_silicon.interface = layer_siox.interface

## === Data files ===
probe = load4('../noBLM/ch061.refl', back_reflectivity=True)
probeh = load4('../noBLM/ch060.refl', back_reflectivity=True)

# Set instrumental (probe) parameters
probe.background.range(-1e-7, 1e-5)
probeh.background.range(-1e-7, 1e-5)
probe.intensity.range(0.9, 1.05)
probeh.intensity = probe.intensity
probe.theta_offset.range(-0.015, 0.005)
probeh.theta_offset = probe.theta_offset
probe.sample_broadening.range(-0.005, 0.02)
probeh.sample_broadening = probe.sample_broadening

# Define critical edge oversampling for samples that require it
probe.critical_edge(substrate=silicon, surface=d2o)

## === Problem definition ===
## a model object consists of a sample and a probe.

## step = True corresponds to a calculation of the reflectivity from an actual profile
## with microslabbed interfaces.  When step = False, the Nevot-Croce
## approximation is used to account for roughness.  This approximation speeds up
## the calculation tremendously, and is reasonably accuarate as long as the
## roughness is much less than the layer thickness
step = False

model = Experiment(sample=sample, probe=probe, dz=stepsize, step_interfaces = step)
modelh = Experiment(sample=sampleh, probe=probeh, dz=stepsize, step_interfaces = step)

problem = FitProblem([model, modelh])

## === Export objects for post analysis ===
problem.name = "DOPC bilayer on TiOx substrate"
#problem.bilayers = [blm]
problem.dimension = dimension
problem.stepsize = stepsize
