
#use %%writefile command in 1st line of cell

#%load tiox_dopc_both.py 



#Added something
from refl1d.names import *
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import os.path

## === Data files ===
probe = load4('ch061.refl', back_reflectivity=False)
probeh = load4('ch060.refl', back_reflectivity=False)
#probe = Probe(T=numpy.linspace(0.18240, 5.7224, 251), L=5.0000)

# Background parameter
# probe.background.value = 0.0000
probe.background.range(-1e-7, 1e-5)
probeh.background.range(-1e-7, 1e-5)

vf_tails = Parameter(name='volume fraction bilayer', value=0.9).range(0.0, 1.0)
vf_ihg = Parameter(name='volume fraction inner headgroups', value=0.4).range(0.0, 1.0)
vf_ohg = Parameter(name='volume fraction outer headgroups', value=0.4).range(0.0, 1.0)
l_hg = Parameter(name='headgroup thickness', value=10.0)

## === Stack ===
##
## First, we create a 'material' for each layer, which has an real and imaginary
## scattering length density, stored in a Refl1d object called 'SLD'
d2o = SLD(name='d2o', rho=6.3000, irho=0.0000)
h2o = SLD(name='h2o', rho=-0.56, irho=0.0000)
hg = SLD(name='headgroups', rho=1.8131, irho=0.000)
tails = SLD(name='lipid tails', rho=-0.2145, irho=0.000)
tiox = SLD(name='tiox', rho=2.1630, irho=0.0000)
siox = SLD(name='siox', rho=4.1000, irho=0.0000)
silicon = SLD(name='silicon', rho=2.0690, irho=0.0000)

buffer = d2o
#buffer = h2o

## Then layers are created, each with its own 'material'.  If you want to force
## two layers to always match SLD you can use the same material in multiple layers.
## The roughnesses of each layer are set to zero to begin with:
layer_water = Slab(material=buffer, thickness=0.0000, interface=5.0000)
layer_ohg = Slab(material=SLD(name='ohg + buffer', rho=vf_ohg*hg.rho + (1 - vf_ohg)*buffer.rho, irho=0.0), thickness=l_hg, interface=5.000)
layer_tails = Slab(material=SLD(name='tails + buffer', rho=vf_tails*tails.rho + (1 - vf_tails)*buffer.rho, irho=0.0), thickness=30.0, interface=5.000)
layer_ihg = Slab(material=SLD(name='ihg + buffer', rho=vf_ihg*hg.rho + (1 - vf_ihg)*buffer.rho, irho=0.0), thickness=l_hg, interface=5.000)
layer_subwater = Slab(material=buffer, thickness=15, interface=10.000)
layer_tiox = Slab(material=tiox, thickness=166.38, interface=10.000)
layer_siox = Slab(material=siox, thickness=7.5804, interface=10.000)
layer_silicon = Slab(material=silicon, thickness=0.0000, interface=0.0000)

sample = Stack()
sample.add(layer_water)
sample.add(layer_ohg)
sample.add(layer_tails)
sample.add(layer_ihg)
sample.add(layer_subwater)
sample.add(layer_tiox)
sample.add(layer_siox)
sample.add(layer_silicon)

bufferh = h2o
#buffer = h2o

## Then layers are created, each with its own 'material'.  If you want to force
## two layers to always match SLD you can use the same material in multiple layers.
## The roughnesses of each layer are set to zero to begin with:
layerh_water = Slab(material=bufferh, thickness=0.0000, interface=5.0000)
layerh_ohg = Slab(material=SLD(name='ohg + bufferh', rho=vf_ohg*hg.rho + (1 - vf_ohg)*bufferh.rho, irho=0.0), thickness=l_hg, interface=5.000)
layerh_tails = Slab(material=SLD(name='tails + bufferh', rho=vf_tails*tails.rho + (1 - vf_tails)*bufferh.rho, irho=0.0), thickness=30.0, interface=5.000)
layerh_ihg = Slab(material=SLD(name='ihg + bufferh', rho=vf_ihg*hg.rho + (1 - vf_ihg)*bufferh.rho, irho=0.0), thickness=l_hg, interface=5.000)
layerh_subwater = Slab(material=bufferh, thickness=15, interface=10.000)

sampleh = Stack()
sampleh.add(layerh_water)
sampleh.add(layerh_ohg)
sampleh.add(layerh_tails)
sampleh.add(layerh_ihg)
sampleh.add(layerh_subwater)
sampleh.add(layer_tiox)
sampleh.add(layer_siox)
sampleh.add(layer_silicon)

## can also be specified as:
# sample = layer0 | layer1 | layer2 | layer3
  
## === Constraints ===
## thickness, interface (roughness) etc. are parameters and
## can be constrained, e.g.
layer_tiox.interface = layer_siox.interface
layer_water.interface = layer_ohg.interface = layer_tails.interface = layer_ihg.interface
layerh_water.interface = layerh_ohg.interface = layerh_tails.interface = layerh_ihg.interface = layer_ihg.interface
layerh_subwater.thickness = layer_subwater.thickness
layerh_tails.thickness = layer_tails.thickness
probeh.intensity = probe.intensity
#ihg.rho = ohg.rho
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
probe.intensity.range(0.95,1.05)

## LAYER RHOs
d2o.rho.range(5.3000, 6.5000)
h2o.rho.range(-0.6, 0.6)
#ohg.rho.range(0, 3)
#tails.rho.range(-0.5, 3.0)
#ihg.rho.range(0, 3)
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
layer_tails.thickness.range(20, 50)
layer_subwater.thickness.range(0, 50)
#layerh_subwater.thickness.range(0, 50)
layer_tiox.thickness.range(66.379, 266.38)
layer_siox.thickness.range(5, 40)
#layer3.thickness.range(0.0000, 100.00)

## LAYER ROUGHNESSES
###################################################################
## the 'interface' associated with layer0 is the boundary between #
## layer0 and layer1, and similarly for layer(N) and layer(N+1)   #
###################################################################
layer_ihg.interface.range(0.0000, 15.000)
layer_subwater.interface.range(0.0000, 15.000)
layer_tiox.interface.range(0.0000, 15.000)
#layer_siox.interface.range(0.0000, 20.000)

## === Problem definition ===
## a model object consists of a sample and a probe,
## zed is the step size in Angstroms to be used for rendering the profile
## increase zed to speed up the calculation
zed = 1    

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

problem.name = "tiox_dopc_both"
