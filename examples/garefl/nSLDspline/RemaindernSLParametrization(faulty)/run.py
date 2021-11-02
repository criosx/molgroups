from bumps.fitproblem import FitProblem
from refl1d import garefl
from refl1d.names import *

problem = garefl.load('model.so')

problem.penalty_limit = 50

