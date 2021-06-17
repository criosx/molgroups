import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import scipy.fft 
import molgroups as mol
import plot

maxarea = 100
stepsize = 0.5
dimension = 300

DOPC = mol.BLM_quaternary()
na1, nh1, nm1, va1, vm1, vh1, lh1 = 0.00760, 0.00461, 0.000468, 972.00, 98, 331.00, 9.56 
na2, nh2, nm2, va2, vm2, vh2, lh2 = 0, 0, 0, 0, 0, 0, 0 
na3, nh3, nm3, va3, vm3, vh3, lh3 = 0, 0, 0, 0, 0, 0, 0
vc, nc = 0, 0

DOPC.fnInit(va1, na1, vm1, nm1, vh1,nh1, lh1, va2, na2, vm2, nm2, 
                        vh2, nh2, lh2, va3, na3, vm3, nm3, vh3, nh3,lh3, vc, nc)


sigma, bulknsld, startz, l_lipid1, l_lipid2, vf_bilayer = 2.0, 9.4114E-06, 50, 11.6, 11.6, 1
rel_pos = .2
DOPC.headgroup1.fnSet(lh1, rel_pos)
DOPC.headgroup2.fnSet(lh1, rel_pos)
DOPC.fnSet(sigma, bulknsld, startz, l_lipid1, l_lipid2, vf_bilayer)

bilayer2 = mol.BLM_quaternary()
na1, nh1, nm1, va1, vm1, vh1, lh1 = 7.2038E-03, 0.00461, 4.7211E-04, 925, 98.8, 331.00, 9.56
na2, nh2, nm2, va2, vm2, vh2, lh2 = 6.1908e-3, 7.6286e-3, 4.7211E-04, 1025, 98.8, 500, 12.0
na3, nh3, nm3, va3, vm3, vh3, lh3 = 0, 0, 0, 0, 0, 0, 0
vc, nc = 0, 0
bilayer2.fnSet(sigma, bulknsld, startz, l_lipid1, l_lipid2, vf_bilayer, .1)
bilayer2.fnInit(va1, na1, vm1, nm1, vh1,nh1, lh1, va2, na2, vm2, nm2, vh2, nh2, lh2, va3, na3, vm3, nm3, vh3, nh3,lh3, vc, nc)

plot.graphGroups(bilayer2, dimension, stepsize, maxarea, True)
# plot.graphBilayer(bilayer, dimension, stepsize, maxarea, bulknsld, True)
# plot.graphGroups(bilayer.headgroup1, dimension, stepsize, maxarea, True)

