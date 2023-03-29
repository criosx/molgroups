from bumps.names import *
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasmodels.data import load_data

# IMPORT THE DATA USED
data = load_data('sim.dat')

# DEFINE THE MODEL
kernel = load_model('cylinder')

pars = dict(scale=1.0, background=0.0005, sld=3.4, sld_solvent=6.4, radius=14.0,
            length=12.0)

model = Model(kernel, **pars)

# PARAMETER RANGES (ONLY THOSE PARAMETERS ARE FITTED)
model.scale.range(0.001, 0.1)
# model.background.range(0, 1)
model.sld.range(-1, 7)
model.sld_solvent.range(-0.6, 6.4)
model.radius.range(10, 200)
model.length.range(10, 200)

M = Experiment(data=data, model=model)
problem = FitProblem(M)
