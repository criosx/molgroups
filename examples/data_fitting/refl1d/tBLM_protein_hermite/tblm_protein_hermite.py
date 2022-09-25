# === Import section ===
import numpy
from molgroups import mol
from molgroups import components as cmp
from molgroups import lipids
from refl1d.names import load4, Parameter, SLD, Slab, Experiment, FitProblem
from refl1d.flayer import FunctionalProfile

# === Constant definition section ===
# Canvas
DIMENSION = 450
STEPSIZE = 0.5

# Hermite Spline
CONTROLPOINTS = 8
SPACING = 15.0
PENETRATION = 42
dDp = [None] * CONTROLPOINTS
dVf = [None] * CONTROLPOINTS

# SLDS
PROTDEUT = 1.67e-6
PROTNONDEUT = 3.14e-6
NSLDH2O = -0.5666e-6
NSLDD2O = 6.36e-6

# Define bilayer and protein objects
blm = mol.tBLM(tether=lipids.HC18SAc, filler=cmp.bmeSAc, lipids=[lipids.POPC], lipid_nf=[1.0])
protein = mol.Hermite(10)
protein.numberofcontrolpoints = CONTROLPOINTS


# Bilayer profile definition function
def bilayer(z, sigma, bulknsld, global_rough, rho_substrate, nf_tether, mult_tether, l_tether, l_lipid1, l_lipid2,
            vf_bilayer):
    """ Generic tethered bilayer """

    # Scale all SLDs from Refl1D units (1e-6 Ang^-2) to molgroups units (Ang^-2)
    bulknsld = bulknsld * 1e-6
    rho_substrate = rho_substrate * 1e-6

    blm.fnSet(sigma=sigma, bulknsld=bulknsld, global_rough=global_rough, rho_substrate=rho_substrate,
              nf_tether=nf_tether, mult_tether=mult_tether, l_tether=l_tether, l_lipid1=l_lipid1, l_lipid2=l_lipid2,
              vf_bilayer=vf_bilayer)

    # Calculate scattering properties of volume occupied by bilayer
    normarea, area, nsl = blm.fnWriteProfile(z)

    # Fill in the remaining volume with buffer of appropriate nSLD
    nsld = nsl / (normarea * numpy.gradient(z)) + (1.0 - area / normarea) * bulknsld

    # export objects for post analysis, needs to be from this function
    problem.bilayers = [blm]
    problem.dimension = DIMENSION
    problem.stepsize = STEPSIZE
    problem.moldat = blm.fnWritePar2Dict({}, 'bilayer', numpy.arange(DIMENSION) * STEPSIZE)

    # Return nSLD profile in Refl1D units
    return nsld * 1e6


def bilayer_prot(z, sigma, bulknsld, global_rough, rho_substrate, nf_tether, mult_tether, l_tether, l_lipid1,
                    l_lipid2, vf_bilayer, nf_protein, protexchratio, dDp, dVf):
    """ Generic tethered bilayer """

    # Scale all SLDs from Refl1D units (1e-6 Ang^-2) to molgroups units (Ang^-2)
    bulknsld = bulknsld * 1e-6
    rho_substrate = rho_substrate * 1e-6

    blm.fnSet(sigma=sigma, bulknsld=bulknsld, global_rough=global_rough, rho_substrate=rho_substrate,
              nf_tether=nf_tether, mult_tether=mult_tether, l_tether=l_tether, l_lipid1=l_lipid1, l_lipid2=l_lipid2,
              vf_bilayer=vf_bilayer)

    # Calculate scattering properties of volume occupied by bilayer
    normarea, area, nsl = blm.fnWriteProfile(z)

    protein.fnSetNormarea(normarea)
    protSLD = PROTNONDEUT + protexchratio * (bulknsld-NSLDH2O) / (NSLDD2O-NSLDH2O) * (PROTDEUT-PROTNONDEUT)
    protein.fnSetRelative(SPACING, blm.headgroups2[0].fnGetZ() + 0.5 * 9.56 - PENETRATION, dDp, dVf, protSLD,
                          nf_protein)

    z1 = blm.methylenes1[0].z - 0.5 * blm.methylenes1[0].l
    z2 = blm.methyls1[0].z + 0.5 * blm.methyls1[0].l
    lipidvol = 0
    for methylene in blm.methylenes1:
        lipidvol += methylene.vol
    for methyl in blm.methyls1:
        lipidvol += methyl.vol
    lipidvol *= vf_bilayer
    v1 = protein.fnGetVolume(z1, z2) / lipidvol

    lipidvol = 0
    for methylene in blm.methylenes2:
        lipidvol += methylene.vol
    for methyl in blm.methyls2:
        lipidvol += methyl.vol
    lipidvol *= vf_bilayer
    z1 = blm.methyls2[0].z - 0.5 * blm.methyls2[0].l
    z2 = blm.methylenes2[0].z + 0.5 * blm.methylenes2[0].l
    v2 = protein.fnGetVolume(z1, z2) / lipidvol

    blm.fnSet(sigma=sigma, bulknsld=bulknsld, global_rough=global_rough, rho_substrate=rho_substrate,
              nf_tether=nf_tether, mult_tether=mult_tether, l_tether=l_tether, l_lipid1=l_lipid1, l_lipid2=l_lipid2,
              vf_bilayer=vf_bilayer, hc_substitution_1=v1, hc_substitution_2=v2)

    normarea, area, nsl = blm.fnWriteProfile(z)
    area, nsl = protein.fnOverlayProfile(z, area, nsl, normarea)
    # Fill in the remaining volume with buffer of appropriate nSLD
    nsld = nsl / (normarea * numpy.gradient(z)) + (1.0 - area / normarea) * bulknsld

    # export objects for post analysis, needs to be in this function

    # for plotting best-fits
    problem.bilayers = [blm]
    problem.dimension = DIMENSION
    problem.stepsize = STEPSIZE
    # for statistical analysis of molgroups
    moldict1 = blm.fnWritePar2Dict({}, 'bilayer', numpy.arange(DIMENSION) * STEPSIZE)
    moldict2 = protein.fnWritePar2Dict({}, 'protein', numpy.arange(DIMENSION) * STEPSIZE)
    problem.moldat = {**moldict1, **moldict2}

    # Return nSLD profile in Refl1D units
    return nsld * 1e6


# bilayer parameters
vf_bilayer = Parameter(name='vf_bilayer', value=0.9).range(0.0, 1.0)                #volume fraction bilayer
l_lipid1 = Parameter(name='l_lipid1', value=10.0).range(8, 30)                      #inner methylenes
l_lipid2 = Parameter(name='l_lipid2', value=10.0).range(8, 16)                      #outer methylenes
sigma = Parameter(name='sigma_blm', value=5).range(2, 8)                            #bilayer roughness
global_rough = Parameter(name='sigma_substrate', value=5).range(2, 9)               #substrate roughness
d_oxide = Parameter(name='d_oxide', value=10).range(5, 40)                          #silicon oxide thickness
d_Cr = Parameter(name='d_cr', value=40).range(10, 60)                              #chromium thickness
d_gold = Parameter(name='d_gold', value=100).range(130, 150)                        #gold thickness
rough_cr_au = Parameter(name='sigma_cr_au', value=10).range(1, 15.0)                #chromium-gold interface roughness
nf_tether = Parameter(name='nf_tether', value=0.7).range(0.2, 1.0)                  #number fraction tether
mult_tether = Parameter(name='mult_tether', value=2).range(0.1, 4)                  #bme-to-tether ratio
l_tether = Parameter(name='l_tether', value=10).range(3, 18)

# protein parameters
dl_lipid = Parameter(name='dl_lipid', value=0.0).range(-3., 3.0)                    #change in methylene thickness
vf_bilayer_prot = Parameter(name='vf_bilayer_prot', value=0.9).range(0.0, 1.0)      #volume fraction bilayer w/ protein
nf_protein = 1
protexchratio = 0.8
for i in range(len(dDp)):
    dDp[i] = Parameter(name='dDp'+str(i), value=0.0).range(-1 * SPACING / 3., SPACING / 3.)
for i in range(len(dVf)):
    dVf[i] = Parameter(name='dVf'+str(i), value=0.001).range(-0.001, 0.4)

# === Stack ===
# First, we create a 'material' for each bulk layer, which has an real and imaginary
# scattering length density, stored in a Refl1d object called 'SLD'
d2o = SLD(name='d2o', rho=6.3000, irho=0.0000)
h2o = SLD(name='h2o', rho=-0.56, irho=0.0000)
d2o_prot = SLD(name='d2o_prot', rho=6.3000, irho=0.0000)
h2o_prot = SLD(name='h2o_prot', rho=-0.56, irho=0.0000)
siox = SLD(name='siox', rho=4.1000, irho=0.0000)
silicon = SLD(name='silicon', rho=2.0690, irho=0.0000)
cr = SLD(name='chromium', rho=2.7, irho=0.0)
gold = SLD(name='gold', rho=4.4, irho=0.0)

# Then bulk layers are created, each with its own 'material'.  If you want to force
# two layers to always match SLD you can use the same material in multiple layers.
# The roughnesses of each layer are set to zero to begin with:

layer_d2o = Slab(material=d2o, thickness=0.0000, interface=5.0000)
layer_h2o = Slab(material=h2o, thickness=0.0000, interface=5.0000)
layer_d2o_prot = Slab(material=d2o_prot, thickness=0.0000, interface=5.0000)
layer_h2o_prot = Slab(material=h2o_prot, thickness=0.0000, interface=5.0000)
layer_siox = Slab(material=siox, thickness=d_oxide, interface=global_rough)
layer_silicon = Slab(material=silicon, thickness=0.0000, interface=global_rough)
layer_cr = Slab(material=cr, thickness=d_Cr, interface=rough_cr_au)
layer_gold = Slab(material=gold, thickness=d_gold - (blm.substrate.z + 0.5 * blm.substrate.l), interface=0.0000)

# Use the bilayer definition function to generate the bilayer SLD profile, passing in the relevant parameters.
# Note that substrate and bulk SLDs are linked to their respective materials.
mollayer = FunctionalProfile(DIMENSION * STEPSIZE, 0, profile=bilayer, sigma=sigma, bulknsld=d2o.rho,
                             global_rough=global_rough, rho_substrate=gold.rho, nf_tether=nf_tether,
                             mult_tether=mult_tether, l_tether=l_tether, l_lipid1=l_lipid1, l_lipid2=l_lipid2,
                             vf_bilayer=vf_bilayer)

mollayerh = FunctionalProfile(DIMENSION * STEPSIZE, 0, profile=bilayer, sigma=sigma, bulknsld=h2o.rho,
                              global_rough=global_rough, rho_substrate=gold.rho, nf_tether=nf_tether,
                              mult_tether=mult_tether, l_tether=l_tether, l_lipid1=l_lipid1, l_lipid2=l_lipid2,
                              vf_bilayer=vf_bilayer)
mollayer_prot = FunctionalProfile(DIMENSION * STEPSIZE, 0, profile=bilayer_prot, sigma=sigma, bulknsld=d2o_prot.rho,
                                  global_rough=global_rough, rho_substrate=gold.rho, nf_tether=nf_tether,
                                  mult_tether=mult_tether, l_tether=l_tether, l_lipid1=l_lipid1+dl_lipid,
                                  l_lipid2=l_lipid2+dl_lipid, vf_bilayer=vf_bilayer_prot, nf_protein=nf_protein,
                                  protexchratio=protexchratio, dDp=dDp, dVf=dVf)

mollayerh_prot = FunctionalProfile(DIMENSION * STEPSIZE, 0, profile=bilayer_prot, sigma=sigma, bulknsld=h2o_prot.rho,
                                   global_rough=global_rough, rho_substrate=gold.rho, nf_tether=nf_tether,
                                   mult_tether=mult_tether, l_tether=l_tether, l_lipid1=l_lipid1+dl_lipid,
                                   l_lipid2=l_lipid2+dl_lipid, vf_bilayer=vf_bilayer_prot, nf_protein=nf_protein,
                                   protexchratio=protexchratio, dDp=dDp, dVf=dVf)

# Stack the layers into individual samples, using common layer objects for layers that are unchanged between samples
# Always build the sample from the substrate up. If the neutron beam is incident from the substrate side,
# set back_reflectivity = True in the probe definition later.

sample = layer_silicon | layer_siox | layer_cr | layer_gold | mollayer | layer_d2o
sampleh = layer_silicon | layer_siox | layer_cr | layer_gold | mollayerh | layer_h2o
sample_prot = layer_silicon | layer_siox | layer_cr | layer_gold | mollayer_prot | layer_d2o_prot
sampleh_prot = layer_silicon | layer_siox | layer_cr | layer_gold | mollayerh_prot | layer_h2o_prot

# Set sample parameter ranges and constraints between layer properties, if these are not set using parameters previously

# nSLD parameters
d2o.rho.range(5.8000, 6.4000)
h2o.rho.range(-0.6, 0.6)
d2o_prot.rho.range(5.8000, 6.4000)
h2o_prot.rho.range(-0.6, 0.6)
siox.rho.range(2.7000, 3.80)
cr.rho.range(2.7000, 4.15)
gold.rho.range(4.2000, 4.60)

# === Data files ===
probe = load4('kr095_4column.refl', back_reflectivity=True)
probeh = load4('kr096_4column.refl', back_reflectivity=True)
probe_prot = load4('kr097_4column.refl', back_reflectivity=True)
probeh_prot = load4('kr098_4column.refl', back_reflectivity=True)

# Background parameter
# probe.background.value = 0.0000
probe.background.range(-1e-9, 1e-5)
probeh.background.range(-1e-9, 1e-5)
probe_prot.background.range(-1e-9, 1e-5)
probeh_prot.background.range(-1e-9, 1e-5)
probe.intensity.range(0.95, 1.05)
probeh.intensity = probe.intensity
probe_prot.intensity = probe.intensity
probeh_prot.intensity = probe.intensity
probe.theta_offset.range(-0.001, 0.001)
probeh.theta_offset = probe.theta_offset
probe_prot.theta_offset = probe.theta_offset
probeh_prot.theta_offset = probe.theta_offset
probe.sample_broadening.range(-0.005, 0.02)
probeh.sample_broadening = probe.sample_broadening
probe_prot.sample_broadening = probe.sample_broadening
probeh_prot.sample_broadening = probe.sample_broadening

# Define critical edge oversampling for samples that require it (typically D2O only)
# probe.critical_edge(substrate=silicon, surface=d2o)

# === Problem definition ===
# a model object consists of a sample and a probe,

# step = True corresponds to a calculation of the reflectivity from an actual profile
# with microslabbed interfaces.  When step = False, the Nevot-Croce
# approximation is used to account for roughness.  This approximation speeds up
# the calculation tremendously, and is reasonably accuarate as long as the
# roughness is much less than the layer thickness
step = False

model = Experiment(sample=sample, probe=probe, dz=STEPSIZE, step_interfaces=step)
modelh = Experiment(sample=sampleh, probe=probeh, dz=STEPSIZE, step_interfaces=step)
model_prot = Experiment(sample=sample_prot, probe=probe_prot, dz=STEPSIZE, step_interfaces=step)
modelh_prot = Experiment(sample=sampleh_prot, probe=probeh_prot, dz=STEPSIZE, step_interfaces=step)
problem = FitProblem([model, modelh, model_prot, modelh_prot])
