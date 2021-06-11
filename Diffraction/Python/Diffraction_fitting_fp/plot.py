import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append(".")
import molgroups as mol

def graphBumpsResults(filename):
    q, F, dq, Fy= np.loadtxt(filename, skiprows=1).T
    fig, axes = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    # fig.set_size_inches(10, 5.5)
    axes[0].plot(q, Fy, zorder=2)
    axes[0].errorbar(q, F, dq, zorder=1)
    axes[0].legend(['model', 'experiment'], fontsize=8)
    axes[0].set_xlim([q[0], q[-1]])
    axes[0].set_ylabel('|F| (1/Å)', fontsize=10)
    axes[0].tick_params(axis='both', which='major', labelsize=10)
    axes[1].plot(q, (Fy-F)/dq)
    axes[1].plot(q, q*0, color='black')
    axes[1].set_ylabel('(f(q) - F) / dq', fontsize=10)
    axes[1].set_xlabel('q (1/Å)', fontsize=10)
    output_filename = filename.split(".")[0]+".pdf"
    fig.savefig(output_filename)
    plt.draw()

def getnSLD(anSL, aArea, dimension, stepsize):
        anSLD = np.zeros(dimension)
        for i in range(len(anSL)):
            if anSL[i] == 0 or aArea[i] == 0:
                anSLD[i] = 0
            else:
                anSLD[i] = anSL[i]/(aArea[i]*stepsize)
        return anSLD

def graphBilayer(bilayer, dimension, stepsize, maxarea, bulknsld, show=False, savefile=None):
    dd, aArea, anSL = bilayer.fnWriteProfile(np.zeros(dimension), np.zeros(dimension), dimension, stepsize, maxarea)
    anSLD  = np.zeros(dimension)

    for i in range(dimension):
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
    axes[2].legend(['SLD'])
    axes[2].set_ylabel('SLD (1/Å^2)')
    axes[2].set_xlabel('z (Å)')
    axes[2].tick_params(axis='both', which='major')
    if savefile: 
        fig.savefig(savefile + ".png")
    if show: plt.show() 
    else: plt.close()


def graphGroups(obj, dimension, stepsize, maxarea, show=False, savefile=None):
    fig, ax = plt.subplots(1, 3)
    __, aArea, anSL = obj.fnWriteProfile(np.zeros(dimension), np.zeros(dimension), dimension, stepsize, maxarea) 
    
    if isinstance(obj, mol.BLM_quaternary):   
        members = {"headgroup" : ("blm_headgroup1", "blm_headgroup2"), "lipid" : ("blm_lipid1", "blm_lipid2"), 
        "methyl" :("blm_methyl1", "blm_methyl1")}
    if isinstance(obj, mol.PC):
        members = {"carbonyl gylcerol": ["pc_cg"], "phosphate": ["pc_ph"], "choline": ["pc_ch"] }
    else: members = {}
    x = [stepsize * i for i in range(dimension)]

    ax[0].set_ylabel('Area (Å^2')
    ax[0].set_xlabel('z (Å)')
    ax[1].set_ylabel('SL (Å)')
    ax[1].set_xlabel('z (Å)')
    ax[2].set_ylabel('SLD (Å)')
    ax[2].set_xlabel('z (Å)')
    ax[0].plot(x, aArea)
    ax[1].plot(x, anSL)
    ax[2].plot(x, getnSLD(anSL, aArea, dimension, stepsize))
    legend = ["total"]
    
    for label in members:
        aArea = np.zeros(dimension)
        anSL  = np.zeros(dimension)
        anSLD = np.zeros(dimension)
        for name in members[label]:
            group = obj.groups[name]
            _, half_area, half_nSL = group.fnWriteProfile(np.zeros(dimension), np.zeros(dimension), dimension, stepsize, maxarea)
            anSL += half_nSL
            aArea += half_area
        anSLD = getnSLD(anSL, aArea, dimension, stepsize)
        ax[0].plot(x, aArea)
        ax[1].plot(x, anSL)
        ax[2].plot(x, anSLD)

        legend.append(label)
    fig.legend(legend)
    # ax[1].legend(legend)
    if savefile: 
        fig.savefig(savefile + ".png")
    if show: plt.show()
    else: plt.close()

def main():
    graphBumpsResults("molgroups/Diffraction/Python/Diffraction_fitting_fp/T1/run-ff.dat")
    graphBumpsResults("molgroups/Diffraction/Python/Diffraction_fitting_fp/T1_headgroups/run-ff.dat")
    graphBumpsResults("molgroups/Diffraction/Python/Diffraction_fitting_fp/T2_headgroups/run-ff.dat")
    graphBumpsResults("molgroups/Diffraction/Python/Diffraction_fitting_fp/T2/run-ff.dat")
    graphBumpsResults("molgroups/Diffraction/Python/Diffraction_fitting_fp/T3/run-ff.dat")

if __name__ == "__main__":
    main()