# === Import section ===
import numpy
from molgroups import mol
from molgroups import components as cmp
from molgroups import lipids
from refl1d.names import load4, Parameter, SLD, Slab, Experiment, FitProblem
from refl1d.flayer import FunctionalProfile

# === Constant definition section ===
# Canvas
DIMENSION = 700
STEPSIZE = 0.5

# Hermite Spline
CONTROLPOINTS = 6
SPACING = 15.0
PENETRATION = 12
dDp = [None] * CONTROLPOINTS
dVf = [None] * CONTROLPOINTS

# SLDS
PROTDEUT = 12e-6
PROTNONDEUT = 12e-6
NSLDH2O = 9.41
NSLDD2O = 9.30

XRAY = 1.5406

# Define bilayer and protein objects
blm = mol.ssBLM(lipids=[lipids.POPC, lipids.PAPC, lipids.POPE, lipids.DLPE, lipids.DPSM, lipids.PAPS, lipids.POPIP2, lipids.chol, lipids.DOPE], 
               lipid_nf=[0.113, 0.061, 0.044, 0.131, 0.088, 0.161, 0.022, 0.28, 0.10],
               xray_wavelength=XRAY)

protein = mol.Hermite()
protein.numberofcontrolpoints = CONTROLPOINTS
blm_prot = mol.BLMProteinComplex(blms=[blm], proteins=[protein])

dense_polymer = mol.Box2Err()
sparse_polymer = mol.Box2Err()

# Bilayer profile definition function
def bilayer(z, sigma, bulknsld, global_rough, rho_substrate, rho_siox, l_siox, l_dense_polymer, l_sparse_polymer, vf_dense_polymer, vf_sparse_polymer, rho_polymer, 
            l_lipid1, l_lipid2, vf_bilayer):
    """ Generic tethered bilayer """

    # Scale all SLDs from Refl1D units (1e-6 Ang^-2) to molgroups units (Ang^-2)
    bulknsld *= 1e-6
    rho_substrate *= 1e-6
    rho_siox *= 1e-6
    rho_polymer *= 1e-6

    blm.fnSet(sigma=sigma, bulknsld=bulknsld, global_rough=global_rough, rho_substrate=rho_substrate, 
              rho_siox=rho_siox, l_siox=l_siox, l_submembrane=l_dense_polymer+l_sparse_polymer,
              l_lipid1=l_lipid1, l_lipid2=l_lipid2, vf_bilayer=vf_bilayer)
    # Calculate scattering properties of volume occupied by bilayer
    maxarea, area, nsl = blm.fnWriteProfile(z)
    # normarea equals maxarea because of substrate layers, otherwise normarea = blm.normarea
    normarea = maxarea
    
    # add in polymer cushion, protein penetration into cushion is not explicitely treated assuming that the
    # polymer density is always low
    vol_dense_polymer = normarea*l_dense_polymer
    dense_polymer.fnSet(volume=vol_dense_polymer, length=l_dense_polymer, position=blm.siox.z+0.5*blm.siox.length,
                        nSL=rho_polymer*vol_dense_polymer, sigma=sigma, nf=vf_dense_polymer)
    vol_sparse_polymer = normarea*l_sparse_polymer
    sparse_polymer.fnSet(volume=vol_sparse_polymer, length=l_sparse_polymer, position=dense_polymer.z+0.5*dense_polymer.length,
                         nSL=rho_polymer*vol_sparse_polymer, sigma=sigma, nf=vf_sparse_polymer)
    area, nsl = dense_polymer.fnOverlayProfile(z, area, nsl, maxarea)
    area, nsl = sparse_polymer.fnOverlayProfile(z, area, nsl, maxarea)
    
    # Fill in the remaining volume with buffer of appropriate nSLD
    nsld = nsl / (normarea * numpy.gradient(z)) + (1.0 - area / normarea) * bulknsld

    # export objects for post analysis, needs to be from this function
    # problem.bilayers = [blm]
    # problem.dimension = DIMENSION
    # problem.stepsize = STEPSIZE
    # problem.moldat = blm.fnWriteGroup2Dict({}, 'bilayer', numpy.arange(DIMENSION) * STEPSIZE)

    # Return nSLD profile in Refl1D units
    return nsld * 1e6


def bilayer_prot(z, sigma, bulknsld, global_rough, rho_substrate, rho_siox, l_siox, l_dense_polymer, l_sparse_polymer, 
                 vf_dense_polymer, vf_sparse_polymer, rho_polymer, l_lipid1, l_lipid2, vf_bilayer, nf_protein, 
                 protexchratio, penetration, dDp, dVf):
    """ Generic tethered bilayer """

    # Scale all SLDs from Refl1D units (1e-6 Ang^-2) to molgroups units (Ang^-2)
    bulknsld *= 1e-6
    rho_substrate *= 1e-6
    rho_siox *= 1e-6
    rho_polymer *= 1e-6

    # Calculate scattering properties of volume occupied by bilayer
    protSLD = PROTNONDEUT + protexchratio * (bulknsld-NSLDH2O) / (NSLDD2O-NSLDH2O) * (PROTDEUT-PROTNONDEUT)
    blm_prot.proteins[0].fnSetRelative(SPACING, blm.headgroups2[0].fnGetZ() + 0.5 * 9.56 - penetration, dDp, dVf, protSLD, nf_protein)
    blm_prot.blms[0].fnSet(sigma=sigma, bulknsld=bulknsld, global_rough=global_rough, rho_substrate=rho_substrate, rho_siox=rho_siox, l_siox=l_siox, 
                           l_submembrane=l_dense_polymer+l_sparse_polymer, l_lipid1=l_lipid1, l_lipid2=l_lipid2, vf_bilayer=vf_bilayer)
    blm_prot.fnAdjustBLMs()
    maxarea, area, nsl = blm_prot.fnWriteProfile(z)
    # normarea equals maxarea because of substrate layers, otherwise normarea = blm_prot.blms[0].normarea, for example
    normarea = maxarea
    
    # add in polymer cushion, protein penetration into cushion is not explicitely treated assuming that the
    # polymer density is always low
    vol_dense_polymer = normarea*l_dense_polymer
    dense_polymer.fnSet(volume=vol_dense_polymer, length=l_dense_polymer, position=blm.siox.z+0.5*blm.siox.length,
                        nSL=rho_polymer*vol_dense_polymer, sigma=sigma, nf=vf_dense_polymer)
    vol_sparse_polymer = normarea*l_sparse_polymer
    sparse_polymer.fnSet(volume=vol_sparse_polymer, length=l_sparse_polymer, position=dense_polymer.z+0.5*dense_polymer.length,
                         nSL=rho_polymer*vol_sparse_polymer, sigma=sigma, nf=vf_sparse_polymer)
    area, nsl = dense_polymer.fnOverlayProfile(z, area, nsl, maxarea)
    area, nsl = sparse_polymer.fnOverlayProfile(z, area, nsl, maxarea)    
    
    # Fill in the remaining volume with buffer of appropriate nSLD
    nsld = nsl / (normarea * numpy.gradient(z)) + (1.0 - area / normarea) * bulknsld

    # export objects for post analysis, needs to be in this function
    # for statistical analysis of molgroups
    moldict1 = blm_prot.blms[0].fnWriteGroup2Dict({}, 'bilayer', numpy.arange(DIMENSION) * STEPSIZE)
    moldict2 = blm_prot.proteins[0].fnWriteGroup2Dict({}, 'protein', numpy.arange(DIMENSION) * STEPSIZE)
    moldict3 = dense_polymer.fnWriteGroup2Dict({}, 'dense_polymer', numpy.arange(DIMENSION) * STEPSIZE)
    moldict4 = sparse_polymer.fnWriteGroup2Dict({}, 'sparse_polymer', numpy.arange(DIMENSION) * STEPSIZE)
    problem.moldat = {**moldict1, **moldict2, **moldict3, **moldict4}
    dict1 = blm_prot.blms[0].fnWriteResults2Dict({}, 'bilayer')
    dict2 = blm_prot.proteins[0].fnWriteResults2Dict({}, 'protein')
    dict3 = dense_polymer.fnWriteResults2Dict({}, 'dense_polymer')
    dict4 = sparse_polymer.fnWriteResults2Dict({}, 'sparse_polymer')
    problem.results = {**dict1, **dict2, **dict3, **dict4}
                        
    # Return nSLD profile in Refl1D units
    return nsld * 1e6

# substrate parameters
global_rough = Parameter(name='sigma_substrate', value=5)              #substrate roughness
d_oxide = Parameter(name='d_oxide', value=10)                          #silicon oxide thickness

# polymer parameters
l_dense_polymer = Parameter(name='l_dense_polymer', value=20).range(10., 100.)
l_sparse_polymer = Parameter(name='l_sparse_polymer', value=50).range(10., 100.)
vf_dense_polymer = Parameter(name='vf_dense_polymer', value=0.4).range(0.02, 0.5)
vf_sparse_polymer = Parameter(name='vf_sparse_polymer', value=0.1).range(0.01, 0.5)

# bilayer parameters
vf_bilayer = Parameter(name='vf_bilayer', value=1.00)                   #volume fraction bilayer
l_lipid1 = Parameter(name='l_lipid1', value=14.0)                       #inner methylenes
l_lipid2 = Parameter(name='l_lipid2', value=14.0)                       #outer methylenes
sigma = Parameter(name='sigma_blm', value=5)                            #bilayer roughness

# protein parameters
protein_penetration = Parameter(name='protein_penetration', value=0).range(-50, 10)
nf_protein = Parameter(name='nf_protein', value=1.).range(0., 1.)
protexchratio = 0.8

for i in range(len(dDp)):
    dDp[i] = Parameter(name='dDp'+str(i), value=0.0)
for i in range(1, len(dVf)-1):
    dVf[i] = Parameter(name='dVf'+str(i), value=0.01).range(0.0, 0.4)
# first and last controlpoint is at zero volume fraction
dVf[0] = dVf[-1] = 0.0

# === Stack ===
# First, we create a 'material' for each bulk layer, which has an real and imaginary
# scattering length density, stored in a Refl1d object called 'SLD'

h2o = SLD(name='h2o', rho=NSLDH2O, irho=0.2)
h2o_prot = SLD(name='h2o_prot', rho=NSLDH2O, irho=0.2)

silicon = SLD(name='silicon', rho=19.67, irho=0.5)
siox = SLD(name='siox', rho=17.73, irho=0.294)
polymer = SLD(name='polymer', rho=11, irho=0.5)

# Then bulk layers are created, each with its own 'material'.  If you want to force
# two layers to always match SLD you can use the same material in multiple layers.
# The roughnesses of each layer are set to zero to begin with:

layer_h2o = Slab(material=h2o, thickness=0.0000, interface=5.0000)
layer_h2o_prot = Slab(material=h2o_prot, thickness=0.0000, interface=5.0000)
layer_siox = Slab(material=siox, thickness=d_oxide, interface=global_rough)
layer_silicon = Slab(material=silicon, thickness=0.0000, interface=global_rough)

# Use the bilayer definition function to generate the bilayer SLD profile, passing in the relevant parameters.
# Note that substrate and bulk SLDs are linked to their respective materials.
mollayer = FunctionalProfile(DIMENSION * STEPSIZE, 0, profile=bilayer, sigma=sigma, bulknsld=h2o.rho, global_rough=global_rough, 
                             rho_substrate=silicon.rho, rho_siox=siox.rho, l_siox=d_oxide, l_dense_polymer=l_dense_polymer, 
                             l_sparse_polymer=l_sparse_polymer, vf_dense_polymer=vf_dense_polymer, vf_sparse_polymer=vf_sparse_polymer, 
                             rho_polymer=polymer.rho, l_lipid1=l_lipid1, l_lipid2=l_lipid2, vf_bilayer=vf_bilayer)
mollayer_prot = FunctionalProfile(DIMENSION * STEPSIZE, 0, profile=bilayer_prot, sigma=sigma, bulknsld=h2o_prot.rho, global_rough=global_rough, 
                                  rho_substrate=silicon.rho, rho_siox=siox.rho, l_siox=d_oxide, l_dense_polymer=l_dense_polymer, 
                                  l_sparse_polymer=l_sparse_polymer, vf_dense_polymer=vf_dense_polymer, vf_sparse_polymer=vf_sparse_polymer, 
                                  rho_polymer=polymer.rho, l_lipid1=l_lipid1, l_lipid2=l_lipid2, vf_bilayer=vf_bilayer, nf_protein=nf_protein, 
                                  protexchratio=protexchratio, penetration=protein_penetration, dDp=dDp, dVf=dVf)

# Stack the layers into individual samples, using common layer objects for layers that are unchanged between samples
# Always build the sample from the substrate up. If the neutron beam is incident from the substrate side,
# set back_reflectivity = True in the probe definition later.

sample = layer_silicon | mollayer | layer_h2o
sample_prot = layer_silicon | mollayer_prot | layer_h2o_prot

# Set sample parameter ranges and constraints between layer properties, if these are not set using parameters previously

# nSLD parameters
# h2o.rho.range(-0.6, 0.6)
# h2o_prot.rho.range(-0.6, 0.6)
# siox.rho.range(2.7000, 3.80)

# === Data files ===
probe = load4('IvsQ_76529+76544_76530_76531.dat', back_reflectivity=True)
probe_prot = load4('IvsQ_76548_76549_76550.dat', back_reflectivity=True)

# Background parameter
# probe.background.value = 0.0000
probe.background.range(-1e-9, 1e-5)
probe_prot.background.range(-1e-9, 1e-5)

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
model_prot = Experiment(sample=sample_prot, probe=probe_prot, dz=STEPSIZE, step_interfaces=step)
problem = FitProblem([model, model_prot])
