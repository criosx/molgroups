from bumps.names import *
from sasmodels.core import load_model
from sasmodels.bumps_model import Model, Experiment
from sasmodels.data import load_data, plot_data

# IMPORT THE DATA USED
data = load_data('sim.dat')

#setattr(data, 'qmin', 0.01)
#setattr(data, 'qmax', 10.0)

# DEFINE THE MODEL
kernel = load_model('ellipsoid@hayter_msa')

pars = dict(scale=1.0, background=0.0005, sld=3.4, sld_solvent=6.4, radius_polar=14.0,
            radius_equatorial=12.0, volfraction=0.075, charge=5, temperature=298.0,
            concentration_salt=0.150, dielectconst=71.08)

model = Model(kernel, **pars)

# PARAMETER RANGES (ONLY THOSE PARAMETERS ARE FITTED)
# model.scale.range(0.1, 7.)
# model.background.range(0, 1)
# model.sld.range(-2, 10)
# model.sld_solvent.range(0, 5)
model.radius_polar.range(17, 35)
model.radius_equatorial.range(5, 17)
model.volfraction.range(0.0001,0.02)
# model.charge.range(0, 20)
# model.temperature.range(0, 1000)
# model.concentration_salt.range(0, 1)
# model.dielectconst.range(0,100)

M = Experiment(data=data, model=model)
problem = FitProblem(M)
