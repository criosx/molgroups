from bumps.names import *
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasmodels.data import load_data

# IMPORT THE DATA USED
data = load_data('sim.dat')

# DEFINE THE MODEL
kernel = load_model('core_shell_cylinder')

pars = dict(scale=1.0, background=0.0005, sld_core=3.4, sld_shell=3.4, sld_solvent=6.4, radius=14.0, thickness=10.0,
            length=12.0)

model = Model(kernel, **pars)

# PARAMETER RANGES (ONLY THOSE PARAMETERS ARE FITTED)
model.scale.range(0.001, 1.)
# model.background.range(0, 1)
model.sld_core.range(-1, 7)
model.sld_shell.range(-1, 7)
model.sld_solvent.range(-0.6, 6.4)
model.radius.range(10, 200)
model.thickness.range(10,200)
model.length.range(10, 200)

M = Experiment(data=data, model=model)
problem = FitProblem(M)
