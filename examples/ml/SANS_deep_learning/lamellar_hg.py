from bumps.names import *
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasmodels.data import load_data

# IMPORT THE DATA USED
data = load_data('sim.dat')

# DEFINE THE MODEL
kernel = load_model('lamellar_hg')

pars = dict(scale=1.0, background=0.0005, length_tail=10.0, length_head=15.0, sld=3.4, sld_head=3.4, sld_solvent=6.4)

model = Model(kernel, **pars)

# PARAMETER RANGES (ONLY THOSE PARAMETERS ARE FITTED)
model.scale.range(0.001, 0.1)
# model.background.range(0, 1)
model.length_tail.range(5., 30.)
model.length_head.range(5., 15.)
model.sld.range(-1, 7)
model.sld_head.range(-1, 7)
model.sld_solvent.range(-0.6, 6.4)

M = Experiment(data=data, model=model)
problem = FitProblem(M)
