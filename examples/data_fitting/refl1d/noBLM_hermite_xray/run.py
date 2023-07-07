# === Import section ===
import numpy as np
from molgroups import mol
from refl1d.names import load4, Parameter, SLD, Slab, Experiment, FitProblem
from refl1d.flayer import FunctionalProfile

CONTROLPOINTS = 8
SPACING = 15.0
dDp = [None] * CONTROLPOINTS
dVf = [None] * CONTROLPOINTS


# === Film structure definition section ===
protein = mol.Hermite()
protein.numberofcontrolpoints = CONTROLPOINTS
protein.xray = True

# molecular layer profile definition function
# traditionally called bilayer, although there is no bilayer in this example
def canvas(z, sigma, rho_bulk, irho_bulk, rho_substrate, irho_substrate, rho_protein, irho_protein, nf_protein, penetration, dDp, dVf):
    # Set unused parameters

    # Scale all SLDs from Refl1D units (1e-6 Ang^-2) to molgroups units (Ang^-2)
    rho_bulk *= 1e-6
    irho_bulk *= 1e-6
    rho_substrate *= 1e-6
    irho_substrate *= 1e-6
    rho_protein *= 1e-6
    irho_protein *= 1e-6

    # protein start position on canvas
    startz = 20.
    substrate.fnSet(nSL=rho_substrate*substrate.vol, sigma=sigma, position=0)
    protein.fnSetRelative(SPACING, startz-penetration, dDp, dVf, rho_protein, nf_protein)

    # Calculate scattering properties of volume occupied by bilayer
    normarea, area, nsl = substrate.fnWriteProfile(z)
    area1 = area
    insl = (normarea * np.gradient(z)) * (area / normarea) * irho_substrate

    area, nsl = protein.fnOverlayProfile(z, area, nsl, normarea)
    insl += (normarea * np.gradient(z)) * ((area-area1) / normarea) * irho_protein

    # Fill in the remaining volume with buffer of appropriate nSLD
    nsld = nsl / (normarea * np.gradient(z)) + (1.0 - area / normarea) * rho_bulk
    irho = insl / (normarea * np.gradient(z)) * (1.0 - area / normarea) * irho_bulk

    # === Export objects for post analysis ===
    dict1 = substrate.fnWriteGroup2Dict({}, 'substrate', np.arange(dimension) * stepsize)
    dict2 = protein.fnWriteGroup2Dict({}, 'protein', np.arange(dimension) * stepsize)
    problem.moldat = {**dict1, **dict2}
    dict3 = substrate.fnWriteResults2Dict({}, 'substrate')
    dict4 = protein.fnWriteResults2Dict({}, 'protein')
    problem.results = {**dict3, **dict4}

    # Return nSLD profile in Refl1D units
    return nsld*1e6 + irho*1e6 * 1j


# Define molecule parameters
rho_protein = Parameter(name='rho_protein', value=10).range(8., 15.)
irho_protein = Parameter(name='irho_protein', value=0).range(0., 2.)
h2o_roughness = Parameter(name='h2o_roughnes', value=2.5).range(2., 4.)
penetration = Parameter(name='penetration', value=0.0).range(-15.0, 10.0)

nf_protein = 1
protexchratio = 0.8
for i in range(len(dDp)):
    dDp[i] = Parameter(name='dDp'+str(i), value=0.0).range(-1 * SPACING / 3., SPACING / 3.)
for i in range(1, len(dVf)-1):
    dVf[i] = Parameter(name='dVf'+str(i), value=0.001).range(-0.001, 0.95)
# first and last controlpoint is at zero volume fraction
dVf[0] = dVf[-1] = 0.0

# Define molecule object
substrate = mol.Box2Err(dz=10, dsigma1=0, dsigma2=2, dlength=40, dvolume=4000, dnSL=0, dnumberfraction=1, name='air')
substrate.xray = True

# Define molgroups space.
dimension = 300       # Number of steps
stepsize = 0.5        # Length of steps

# === Stack ===
#
# First, we create a 'material' for each bulk layer, which has a real and imaginary
# scattering length density, stored in a Refl1d object called 'SLD'
air = SLD(name='air', rho=0.00, irho=0.00)
h2o = SLD(name='h2o', rho=9.5, irho=0.00)

# Then bulk layers are created, each with its own 'material'.  If you want to force
# two layers to always match SLD you can use the same material in multiple layers.
# The roughnesses of each layer are set to zero to begin with:

layer_air = Slab(material=air, thickness=0.0000, interface=2.0000)
layer_h2o = Slab(material=h2o, thickness=0.0000, interface=2.0000)

# Use the canvas definition function to generate the bilayer SLD profile, passing in the relevant parameters.
# Note that substrate and bulk SLDs are linked to their respective materials.
mollayer = FunctionalProfile(dimension*stepsize,
                             0,
                             profile=canvas,
                             sigma=h2o_roughness,
                             rho_bulk=h2o.rho,
                             irho_bulk=h2o.irho,
                             rho_substrate=air.rho,
                             irho_substrate=air.irho,
                             rho_protein=rho_protein,
                             irho_protein=irho_protein,
                             nf_protein=nf_protein,
                             penetration=penetration,
                             dDp=dDp,
                             dVf=dVf
                             )

# Stack the layers into individual samples, using common layer objects for layers that are unchanged between samples
# As a convention, always build the sample from the substrate up. If the neutron beam is incident from the substrate
# side, set back_reflectivity = True in the probe definition later.

sample = layer_air | mollayer | layer_h2o

# Set sample parameter ranges and constraints between layer properties

# nSLD parameters
h2o.rho.range(9.0, 10.0)
h2o.irho.range(0.0, 2.0)

# layer thickness parameters
# layer_siox.thickness.range(20, 80)

# layer roughness parameters
###################################################################
# the 'interface' associated with layer0 is the boundary between #
# layer0 and layer1, and similarly for layer(N) and layer(N+1)   #
###################################################################
# layer_siox.interface.range(2.0000, 3.000)

# Si and SiOx roughnesses are the same
# layer_silicon.interface = layer_siox.interface

# === Data files ===
probe = load4('data0.txt', back_reflectivity=True)

# Set instrumental (probe) parameters
probe.background.range(-1e-10, 1e-5)
probe.intensity.range(0.8, 1.2)
probe.theta_offset.range(-0.015, 0.005)
# probe.sample_broadening.range(-0.005, 0.02)
# probeh.sample_broadening = probe.sample_broadening

# Define critical edge oversampling for samples that require it
probe.critical_edge(substrate=air, surface=h2o)

# === Problem definition ===
# a model object consists of a sample and a probe.

# step = True corresponds to a calculation of the reflectivity from an actual profile
# with microslabbed interfaces.  When step = False, the Nevot-Croce
# approximation is used to account for roughness.  This approximation speeds up
# the calculation tremendously, and is reasonably accuarate as long as the
# roughness is much less than the layer thickness
step = True

model = Experiment(sample=sample, probe=probe, dz=stepsize, step_interfaces = step)

problem = FitProblem([model])
