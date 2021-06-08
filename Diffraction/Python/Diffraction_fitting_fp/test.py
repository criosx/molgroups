import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import scipy.fft 
import molgroups as mol

maxarea = 100
stepsize = 0.5
dimension = 300
aArea = np.zeros(dimension)
anSL  = np.zeros(dimension)
anSLD = np.zeros(dimension)

bilayer = mol.BLM_quaternary()
na1, nh1, nm1, va1, vm1, vh1, lh1 = 0.00760, 0.00461, 0.000468, 972.00, 98, 331.00, 9.56 
na2, nh2, nm2, va2, vm2, vh2, lh2 = 0, 0, 0, 0, 0, 0, 0 
na3, nh3, nm3, va3, vm3, vh3, lh3 = 0, 0, 0, 0, 0, 0, 0
vc, nc = 0, 0

bilayer.fnInit(va1, na1, vm1, nm1, vh1,nh1, lh1, va2, na2, vm2, nm2, 
                        vh2, nh2, lh2, va3, na3, vm3, nm3, vh3, nh3,lh3, vc, nc)

sigma, bulknsld, startz, l_lipid1, l_lipid2, vf_bilayer = 2.0, 9.4114E-06, 50, 11.6, 11.6, 1

bilayer.fnSet(sigma, bulknsld, startz, l_lipid1, l_lipid2, vf_bilayer)

dd, aArea, anSL = bilayer.fnWriteProfile(aArea, anSLD, dimension, stepsize, maxarea)

for i in range(len(aArea)):
    if aArea[i] != 0:
        anSLD[i] = anSL[i] / (aArea[i]*stepsize) * aArea[i]/dd + bulknsld * (1 - aArea[i]/dd)
    else:
        anSLD[i] = bulknsld

x = [stepsize * i for i in range(dimension)]

fig, axes = plt.subplots(3)

axes[0].plot(x, anSL)
axes[0].legend(['SL'])
axes[0].set_ylabel('SL (Å)')
axes[0].set_xlabel('z (Å)')
axes[0].tick_params(axis='both', which='major')

axes[1].plot(x, aArea)
axes[1].legend(['Area'])
axes[1].set_ylabel('Area (Å^2')
axes[1].set_xlabel('z (Å)')
axes[1].tick_params(axis='both', which='major')

axes[2].plot(x, anSLD)
axes[2].set_ylabel('SLD (1/Å^2)')
axes[2].set_xlabel('z (Å)')
axes[2].tick_params(axis='both', which='major')
plt.show()


fig2, ax = plt.subplots(1, 2)
counter = 0 
bm = {"headgroup" : ("blm_headgroup1", "blm_headgroup2"), "lipid" : ("blm_lipid1", "blm_lipid2"), 
"methyl" :("blm_methyl1", "blm_methyl1")}
ax[0].set_ylabel('Area (Å^2')
ax[0].set_xlabel('z (Å)')
ax[1].set_ylabel('SL (Å)')
ax[1].set_xlabel('z (Å)')
ax[0].plot(x, aArea)
ax[1].plot(x, anSL)
legend = ["total"]
for label in bm:
    aArea = np.zeros(dimension)
    anSL  = np.zeros(dimension)
    anSLD = np.zeros(dimension)
    for name in bm[label]:
        group = bilayer.groups[name]
        _, half_area, half_nSL = group.fnWriteProfile(np.zeros(dimension), np.zeros(dimension), dimension, stepsize, maxarea)
        anSL += half_nSL
        aArea += half_area
    ax[0].plot(x, aArea)
    ax[1].plot(x, anSL)
    legend.append(label)
ax[0].legend(legend)
ax[1].legend(legend)
plt.show() 

