# === Import section ===
import numpy as np
from molgroups import mol
from refl1d.names import load4, Parameter, SLD, Slab, Experiment, FitProblem
from refl1d.flayer import FunctionalProfile

CONTROLPOINTS = 8
SPACING = 15.0
dDp = [0.0] * CONTROLPOINTS
dVf = [0.0] * CONTROLPOINTS


# === Film structure definition section ===
protein = mol.Hermite()
protein.numberofcontrolpoints = CONTROLPOINTS

# molecular layer profile definition function
# traditionally called bilayer, although there is no bilayer in this example
def canvas(z, sigma, rho_bulk, rho_substrate, rho_protein, rho_protein_exchanged, protexchratio, nf_protein, penetration, dDp, dVf):
    # Set unused parameters

    # Scale all SLDs from Refl1D units (1e-6 Ang^-2) to molgroups units (Ang^-2)
    rho_bulk *= 1e-6
    rho_substrate *= 1e-6
    rho_protein *= 1e-6
    rho_protein_exchanged *= 1e-6

    # protein start position on canvas
    substrate.fnSet(nSL=rho_substrate*substrate.vol, sigma=sigma, position=0)
    startz = substrate.z + 0.5 * substrate.length
    rho_protein_mod = (6.4e-6 - rho_bulk) / 7.0e-6 * rho_protein + (rho_bulk + 0.6e-6) / 7.0e-6 * rho_protein_exchanged
    rho_protein_mod = (rho_protein_mod - rho_protein) * protexchratio + rho_protein
    protein.fnSetRelative(SPACING, startz+penetration, dDp, dVf, rho_protein_mod, nf_protein)

    # Calculate scattering properties of volume occupied by bilayer
    normarea, area, nsl = substrate.fnWriteProfile(z)

    protein.fnSetNormarea(normarea)
    area, nsl = protein.fnOverlayProfile(z, area, nsl, normarea)

    # Fill in the remaining volume with buffer of appropriate nSLD
    nsld = nsl / (normarea * np.gradient(z)) + (1.0 - area / normarea) * rho_bulk

    # === Export objects for post analysis ===
    dict1 = substrate.fnWriteGroup2Dict({}, 'substrate', np.arange(dimension) * stepsize)
    dict2 = protein.fnWriteGroup2Dict({}, 'protein', np.arange(dimension) * stepsize)
    problem.moldat = {**dict1, **dict2}
    dict3 = substrate.fnWriteResults2Dict({}, 'substrate')
    dict4 = protein.fnWriteResults2Dict({}, 'protein')
    problem.results = {**dict3, **dict4}

    # Return nSLD profile in Refl1D units
    return nsld*1e6


# Define molecule parameters
nf_protein = 1
protexchratio = 0.8
rho_protein = 1.8
rho_protein_exchanged = 3.4
h2o_roughness = Parameter(name='h2o_roughnes', value=2.5).range(2., 4.)
penetration = Parameter(name='penetration', value=0.0).range(-15.0, 10.0)
nf_protein_h2o = Parameter(name='nf_protein_h2o', value=1.0).range(0.8, 1.0)


for i in range(len(dDp)):
    dDp[i] = Parameter(name='dDp'+str(i), value=0.0).range(-1 * SPACING / 3., SPACING / 3.)
for i in range(1, len(dVf)-1):
    dVf[i] = Parameter(name='dVf'+str(i), value=0.001).range(-0.001, 0.95)
# first and last controlpoint is at zero volume fraction
dVf[0] = dVf[-1] = 0.0

# Define molecule object
substrate = mol.Box2Err(dz=0, dsigma1=0, dsigma2=2, dlength=40, dvolume=4000, dnSL=0, dnumberfraction=1, name='air')

# Define molgroups space.
dimension = 300       # Number of steps
stepsize = 0.5        # Length of steps

# === Stack ===
#
# First, we create a 'material' for each bulk layer, which has a real and imaginary
# scattering length density, stored in a Refl1d object called 'SLD'
air = SLD(name='air', rho=0.00, irho=0.00)
d2o = SLD(name='d2o', rho=6.4, irho=0.00)
h2o = SLD(name='h2o', rho=-0.56, irho=0.00)

# Then bulk layers are created, each with its own 'material'.  If you want to force
# two layers to always match SLD you can use the same material in multiple layers.
# The roughnesses of each layer are set to zero to begin with:

layer_air = Slab(material=air, thickness=0.0000, interface=2.0000)
layer_d2o = Slab(material=d2o, thickness=0.0000, interface=2.0000)
layer_h2o = Slab(material=h2o, thickness=0.0000, interface=2.0000)

# Use the canvas definition function to generate the bilayer SLD profile, passing in the relevant parameters.
# Note that substrate and bulk SLDs are linked to their respective materials.
mollayer0 = FunctionalProfile(dimension*stepsize,
                              0,
                              profile=canvas,
                              sigma=h2o_roughness,
                              rho_bulk=d2o.rho,
                              rho_substrate=air.rho,
                              rho_protein=rho_protein,
                              rho_protein_exchanged = rho_protein_exchanged,
                              protexchratio = protexchratio,
                              nf_protein=nf_protein,
                              penetration=penetration,
                              dDp=dDp,
                              dVf=dVf
                              )

mollayer1 = FunctionalProfile(dimension*stepsize,
                              0,
                              profile=canvas,
                              sigma=h2o_roughness,
                              rho_bulk=h2o.rho,
                              rho_substrate=air.rho,
                              rho_protein=rho_protein,
                              rho_protein_exchanged = rho_protein_exchanged,
                              protexchratio=protexchratio,
                              nf_protein=nf_protein_h2o,
                              penetration=penetration,
                              dDp=dDp,
                              dVf=dVf
                              )

# Stack the layers into individual samples, using common layer objects for layers that are unchanged between samples
# As a convention, always build the sample from the substrate up. If the neutron beam is incident from the substrate
# side, set back_reflectivity = True in the probe definition later.

sample0 = layer_air | mollayer0 | layer_d2o
sample1 = layer_air | mollayer1 | layer_h2o

# Set sample parameter ranges and constraints between layer properties

# nSLD parameters
d2o.rho.range(5.0, 6.4)
h2o.rho.range(-0.56, 0.5)
# h2o.irho.range(0.0, 2.0)

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
probe0 = load4('data0.txt', back_reflectivity=True)
probe1 = load4('data1.txt', back_reflectivity=True)

# Set instrumental (probe) parameters
probe0.background.range(-1e-9, 1e-6)
probe0.intensity.range(0.8, 1.2)
probe0.theta_offset.range(-0.015, 0.005)
probe1.background.range(-1e-9, 1e-5)
probe1.intensity = probe0.intensity
probe1.theta_offset = probe0.theta_offset
# probe.sample_broadening.range(-0.005, 0.02)
# probeh.sample_broadening = probe.sample_broadening

# Define critical edge oversampling for samples that require it
probe0.critical_edge(substrate=air, surface=d2o)
probe1.critical_edge(substrate=air, surface=h2o)

# === Problem definition ===
# a model object consists of a sample and a probe.

# step = True corresponds to a calculation of the reflectivity from an actual profile
# with microslabbed interfaces.  When step = False, the Nevot-Croce
# approximation is used to account for roughness.  This approximation speeds up
# the calculation tremendously, and is reasonably accuarate as long as the
# roughness is much less than the layer thickness
step = True

model0 = Experiment(sample=sample0, probe=probe0, dz=stepsize, step_interfaces = step)
model1 = Experiment(sample=sample1, probe=probe1, dz=stepsize, step_interfaces = step)

problem = FitProblem([model0, model1])
