# === Import section ===
import sys
import numpy as np
from molgroups import mol
from refl1d.names import load4, Parameter, SLD, Slab, Experiment, FitProblem
from refl1d.flayer import FunctionalProfile

# === Film structure definition section ===

# molecular layer profile definition function
# traditionally called bilayer, although there is no bilayer in this example


def bilayer(z, sigma, bulknsld, global_rough, rho_substrate, l_surfasil, vf_surfasil):
    """ Fairly generic bilayer. This assumes a stack of materials already existing because siox.l is set to zero """

    # Set unused parameters

    # Scale all SLDs from Refl1D units (1e-6 Ang^-2) to molgroups units (Ang^-2)
    bulknsld *= 1e-6
    rho_substrate *= 1e-6

    substrate.fnSet(nSL=rho_substrate*substrate.vol, sigma=global_rough, position=0)

    volume = l_surfasil * 100 * vf_surfasil
    surfasil.fnSet(length=l_surfasil, position=20+0.5*l_surfasil, nf=1, volume=volume,
                   nSL=0.24e-6 * volume)
    surfasil.fnSetSigma(sigma1=global_rough, sigma2=sigma)

    # Calculate scattering properties of volume occupied by bilayer
    normarea, area, nsl = substrate.fnWriteProfile(z)
    area, nsl = surfasil.fnOverlayProfile(z, area, nsl, normarea)

    # Fill in the remaining volume with buffer of appropriate nSLD
    nsld = nsl / (normarea * np.gradient(z)) + (1.0 - area / normarea) * bulknsld

    # === Export objects for post analysis ===
    problem.name = "Surfasil on SiOx"
    problem.groups = [substrate, surfasil]
    problem.dimension = dimension
    problem.stepsize = stepsize
    dict1 = substrate.fnWriteGroup2Dict({}, 'substrate', np.arange(dimension) * stepsize)
    dict2 = substrate.fnWriteGroup2Dict({}, 'surfasil', np.arange(dimension) * stepsize)
    problem.moldat = {**dict1, **dict2}


    # Return nSLD profile in Refl1D units
    return nsld*1e6


# Define bilayer parameters
l_surfasil = Parameter(name='surfasil thickness', value=60).range(20, 70)
vf_surfasil = Parameter(name='surfasil volfrac', value=0.9).range(0.8, 1.)
sigma = Parameter(name='surfasil roughness', value=5).range(2, 10)
# global_rough = Parameter(name ='global roughness', value=5).range(2, 3)

# Define bilayer object
substrate = mol.Box2Err(dz=10, dsigma1=0, dsigma2=2, dlength=40, dvolume=4000, dnSL=0, dnumberfraction=1, name='siox')
surfasil = mol.Box2Err(dz=50, dsigma1=2, dsigma2=2, dlength=60, dvolume=100, dnSL=0, dnumberfraction=1, name='surfasil')

# Define molgroups space.
dimension = 300       # Number of steps
stepsize = 0.5        # Length of steps

# === Stack ===
#
# First, we create a 'material' for each bulk layer, which has an real and imaginary
# scattering length density, stored in a Refl1d object called 'SLD'
d2o = SLD(name='d2o', rho=6.3000, irho=0.0000)
h2o = SLD(name='h2o', rho=-0.56, irho=0.0000)
siox = SLD(name='siox', rho=3.4000, irho=0.0000)
silicon = SLD(name='silicon', rho=2.0690, irho=0.0000)

# Then bulk layers are created, each with its own 'material'.  If you want to force
# two layers to always match SLD you can use the same material in multiple layers.
# The roughnesses of each layer are set to zero to begin with:

layer_d2o = Slab(material=d2o, thickness=0.0000, interface=2.0000)
layer_h2o = Slab(material=h2o, thickness=0.0000, interface=2.0000)
layer_siox = Slab(material=siox, thickness=7.5804, interface=2.000)
layer_silicon = Slab(material=silicon, thickness=0.0000, interface=0.0000)

# Use the bilayer definition function to generate the bilayer SLD profile, passing in the relevant parameters.
# Note that substrate and bulk SLDs are linked to their respective materials.
mollayer = FunctionalProfile(dimension*stepsize, 0, profile=bilayer, sigma=sigma, bulknsld=d2o.rho,
                             global_rough=layer_siox.interface, rho_substrate=siox.rho, l_surfasil=l_surfasil,
                             vf_surfasil=vf_surfasil)
mollayerh = FunctionalProfile(dimension*stepsize, 0, profile=bilayer, sigma=sigma, bulknsld=h2o.rho,
                              global_rough=layer_siox.interface, rho_substrate=siox.rho, l_surfasil=l_surfasil,
                             vf_surfasil=vf_surfasil)

# Stack the layers into individual samples, using common layer objects for layers that are unchanged between samples
# As a convention, always build the sample from the substrate up. If the neutron beam is incident from the substrate
# side, set back_reflectivity = True in the probe definition later.

sample = layer_silicon | layer_siox | mollayer | layer_d2o
sampleh = layer_silicon | layer_siox | mollayerh | layer_h2o

# Set sample parameter ranges and constraints between layer properties

# nSLD parameters
d2o.rho.range(5.6, 6.4)
h2o.rho.range(4.0, 5.6)
siox.rho.range(3.1000, 3.7000)

# layer thickness parameters
layer_siox.thickness.range(20, 80)

# layer roughness parameters
###################################################################
# the 'interface' associated with layer0 is the boundary between #
# layer0 and layer1, and similarly for layer(N) and layer(N+1)   #
###################################################################
layer_siox.interface.range(2.0000, 3.000)

# Si and SiOx roughnesses are the same
layer_silicon.interface = layer_siox.interface

# === Data files ===
probe = load4('data0.txt', back_reflectivity=True)
probeh = load4('data1.txt', back_reflectivity=True)
#probe = load4('ch060.refl', back_reflectivity=True)
#probeh = load4('ch061.refl', back_reflectivity=True)

# Set instrumental (probe) parameters
probe.background.range(-1e-7, 1e-5)
probeh.background.range(-1e-7, 1e-5)
probe.intensity.range(0.9, 1.05)
probeh.intensity = probe.intensity
probe.theta_offset.range(-0.015, 0.005)
probeh.theta_offset = probe.theta_offset
# probe.sample_broadening.range(-0.005, 0.02)
# probeh.sample_broadening = probe.sample_broadening

# Define critical edge oversampling for samples that require it
probe.critical_edge(substrate=silicon, surface=d2o)

# === Problem definition ===
# a model object consists of a sample and a probe.

# step = True corresponds to a calculation of the reflectivity from an actual profile
# with microslabbed interfaces.  When step = False, the Nevot-Croce
# approximation is used to account for roughness.  This approximation speeds up
# the calculation tremendously, and is reasonably accuarate as long as the
# roughness is much less than the layer thickness
step = True

model = Experiment(sample=sample, probe=probe, dz=stepsize, step_interfaces = step)
modelh = Experiment(sample=sampleh, probe=probeh, dz=stepsize, step_interfaces = step)

problem = FitProblem([model, modelh])

