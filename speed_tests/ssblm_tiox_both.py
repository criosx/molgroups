import sys
# append path to your molgroups, or just link molgroups to your same directory
#sys.path.append('G:\\My Drive\\software\\nr\\molgroups\\Diffraction\\Python\\Diffraction_fitting_fp')
sys.path.append('../Reflectometry')
import molgroups as mol
from refl1d.names import *
from refl1d.flayer import FunctionalProfile
import importlib
importlib.reload(mol)

#print(mol.__file__)

def bilayer(z, sigma, bulknsld, global_rough, rho_substrate, l_submembrane, l_lipid1, l_lipid2, vf_bilayer):
    """ Fairly generic bilayer. This assumes a stack of materials already existing because siox.l is set to zero """
    
    #global blm                     # Alternative to reinstantiating ssBLM_quaternary every time
    #blm = mol.ssBLM_quaternary()    # Default bilayer is DOPC, so we don't need to change any of the molecular volumes
    dimension=len(z)
    stepsize = z[1]-z[0]
    # define canvas
    l_siox = 0.0 # could make a parameter in the future
    rho_siox = 0.0

    bulknsld = bulknsld*1e-6
    rho_siox = rho_siox * 1e-6
    rho_substrate = rho_substrate * 1e-6

    blm.fnSet(sigma, bulknsld, global_rough, rho_substrate, rho_siox, l_siox, l_submembrane, l_lipid1, l_lipid2, vf_bilayer)
    
    normarea, area, nsl = blm.fnWriteProfile(np.zeros_like(z), np.zeros_like(z), dimension, stepsize, 1.0)

    # this replaces fnWriteCanvas2Model
    nsld = nsl / (normarea * stepsize) + (1.0 - area / normarea) * bulknsld

    return nsld * 1e6

## === Data files ===
probe = load4('ch061.refl', back_reflectivity=True)
probeh = load4('ch060.refl', back_reflectivity=True)

# Background parameter
# probe.background.value = 0.0000
probe.background.range(-1e-7, 1e-5)
probeh.background.range(-1e-7, 1e-5)
probe.intensity.range(0.9, 1.05)
probeh.intensity = probe.intensity
probe.theta_offset.range(-0.02, 0.02)
probeh.theta_offset = probe.theta_offset
probe.sample_broadening.range(-0.005, 0.02)
probeh.sample_broadening = probe.sample_broadening

# Define bilayer parameters
vf_bilayer = Parameter(name='volume fraction bilayer', value=0.98).range(0.0, 1.0)
l_lipid1 = Parameter(name='inner acyl chain thickness', value=10.0).range(8, 16)
l_lipid2 = Parameter(name='outer acyl chain thickness', value=10.0).range(8, 16)
sigma = Parameter(name='bilayer roughness', value=5).range(2, 9)
#sigmah = Parameter(name='bilayer roughness h2o', value=5).range(2, 9)
global_rough = Parameter(name ='substrate roughness', value=5).range(2, 9)
l_tiox = Parameter(name='total tiox thickness', value=120).range(50, 150)
l_submembrane = Parameter(name='submembrane thickness', value=10).range(0, 50)
#l_submembraneh = Parameter(name='submembrane thickness h2o', value=10).range(0, 20)

blm = mol.ssBLM_quaternary()        # required to subtract the bilayer length in layer_tiox definition; only really necessary if using "global blm" in bilayer function
blm.volacyllipid=975
dimension=300
stepsize=0.5

## === Stack ===
##
## First, we create a 'material' for each layer, which has an real and imaginary
## scattering length density, stored in a Refl1d object called 'SLD'
d2o = SLD(name='d2o', rho=6.3000, irho=0.0000)
h2o = SLD(name='h2o', rho=-0.56, irho=0.0000)
tiox = SLD(name='tiox', rho=2.1630, irho=0.0000)
siox = SLD(name='siox', rho=4.1000, irho=0.0000)
silicon = SLD(name='silicon', rho=2.0690, irho=0.0000)

## Then layers are created, each with its own 'material'.  If you want to force
## two layers to always match SLD you can use the same material in multiple layers.
## The roughnesses of each layer are set to zero to begin with:
mollayer = FunctionalProfile(dimension*stepsize, 0, profile=bilayer, sigma=sigma,
                                bulknsld=d2o.rho, global_rough=global_rough, rho_substrate=tiox.rho,
                                l_submembrane=l_submembrane, l_lipid1=l_lipid1, l_lipid2=l_lipid2,
                                vf_bilayer=vf_bilayer)
mollayerh = FunctionalProfile(dimension*stepsize, 0, profile=bilayer, sigma=sigma,
                                bulknsld=h2o.rho, global_rough=global_rough, rho_substrate=tiox.rho,
                                l_submembrane=l_submembrane, l_lipid1=l_lipid1, l_lipid2=l_lipid2,
                                vf_bilayer=vf_bilayer)

layer_d2o = Slab(material=d2o, thickness=0.0000, interface=5.0000)
layer_h2o = Slab(material=h2o, thickness=0.0000, interface=5.0000)
layer_tiox = Slab(material=tiox, thickness=l_tiox - blm.substrate.l/2, interface=0.0)
layer_siox = Slab(material=siox, thickness=7.5804, interface=10.000)
layer_silicon = Slab(material=silicon, thickness=0.0000, interface=0.0000)

sample = Stack()
sample.add(layer_silicon)
sample.add(layer_siox)
sample.add(layer_tiox)
sample.add(mollayer)
sample.add(layer_d2o)

sampleh = Stack()
sampleh.add(layer_silicon)
sampleh.add(layer_siox)
sampleh.add(layer_tiox)
sampleh.add(mollayerh)
sampleh.add(layer_h2o)

"""
# speed tests
import time
starttime = time.time()
for _ in range(1):
    bilayer(np.arange(dimension)*stepsize, sigma.value, d2o.rho.value, global_rough.value, tiox.rho.value, l_submembrane.value, l_lipid1.value, l_lipid2.value, vf_bilayer.value)
print(time.time()-starttime)
"""
probe.critical_edge(substrate=silicon, surface=d2o)

## can also be specified as:
# sample = layer0 | layer1 | layer2 | layer3
  
## === Constraints ===
## thickness, interface (roughness) etc. are parameters and
## can be constrained, e.g.
layer_silicon.interface = layer_siox.interface
## (to tie the first layer to have exactly the same thickness as the third layer)
# layer1.interface = layer2.interface
## (to make the roughness between layer1 and layer2 the same as between layer2 and layer3)
# layer0.material = layer4.material
## (make their sld properties match, real and imaginary)
# sld0.rho = sld1.rho
## (to force only the real rho to match for two materials)

## === Fit parameters ===
## "range" specifies a fitting range in terms of min/max value
## "pmp" specifies fitting range in terms of +/-  %
## "pm" specifies fitting range in terms of +/- value

## THETA OFFSET
## this parameter accounts for theta misalignment
## probe.theta_offset.range(-.01,.01)

## INTENSITY
#probe.intensity.range(0.95,1.05)

## LAYER RHOs
d2o.rho.range(5.3000, 6.5000)
h2o.rho.range(-0.6, 0.6)
#ohg.rho.range(-0.6, 6.4)
#tails.rho.range(-0.5, 3.0)
#ihg.rho.range(-0.6, 6.4)
tiox.rho.range(1.1630, 3.1630)
siox.rho.range(3.1000, 5.1000)
#silicon.rho.range(1.0690, 3.0690)

## LAYER ABSORPTIONS (imaginary rho)
#sld0.irho.range(-1.0000, 1.0000)
#sld1.irho.range(-1.0000, 1.0000)
#sld2.irho.range(-1.0000, 1.0000)
#sld3.irho.range(-1.0000, 1.0000)

## LAYER THICKNESSES
#layer_ipa.thickness.range(0.0000, 100.00)
layer_tiox.thickness.range(66.379, 266.38)
layer_siox.thickness.range(5, 40)
#layer3.thickness.range(0.0000, 100.00)

## LAYER ROUGHNESSES
###################################################################
## the 'interface' associated with layer0 is the boundary between #
## layer0 and layer1, and similarly for layer(N) and layer(N+1)   #
###################################################################
#layer_tiox.interface.range(0.0000, 15.000)
layer_siox.interface.range(2.0000, 9.000)

## === Problem definition ===
## a model object consists of a sample and a probe,
## zed is the step size in Angstroms to be used for rendering the profile
## increase zed to speed up the calculation
zed = stepsize

## step = True corresponds to a calculation of the reflectivity from an actual profile
## with microslabbed interfaces.  When step = False, the Nevot-Croce
## approximation is used to account for roughness.  This approximation speeds up
## the calculation tremendously, and is reasonably accuarate as long as the
## roughness is much less than the layer thickness
step = False

model = Experiment(sample=sample, probe=probe, dz=zed, step_interfaces = step)
modelh = Experiment(sample=sampleh, probe=probeh, dz=zed, step_interfaces = step)
## simultaneous fitting: if you define two models
# models = model1, model2
# problem = MultiFitProblem(models=models)

# fitting a single model:
problem = FitProblem([model, modelh])
#problem = FitProblem([model])

problem.name = "tiox_dopc_both"


