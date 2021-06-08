import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append(".")

def graphBumpsResults(filename):
    data = np.loadtxt(filename, delimiter=None)
    q, F, dq, Fy= data[:, 0], data[:, 1], data[:, 2], data[:, 3]
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
    # axes[0].set_ylabel('(f(q) - F) / dq', fontsize=10)
    axes[0].set_xlabel('q (1/Å)', fontsize=10)
    fig.savefig(savefile(filename))
    plt.draw()

def savefile(filename):
    return filename.split(".")[0]+".pdf"

def graphGroups(bilayer):
    for group in bilayer.groups:
        bilayer.groups[group].fnWriteProfile(bilayer.aArea, bilayer.anSL, stepsize, dMaxArea)

def main():
    graphBumpsResults("molgroups/Diffraction/Python/Diffraction_fitting_fp/T1/run-ff.dat")
    graphBumpsResults("molgroups/Diffraction/Python/Diffraction_fitting_fp/T1_headgroups/run-ff.dat")
    graphBumpsResults("molgroups/Diffraction/Python/Diffraction_fitting_fp/T2_headgroups/run-ff.dat")
    graphBumpsResults("molgroups/Diffraction/Python/Diffraction_fitting_fp/T2/run-ff.dat")
    graphBumpsResults("molgroups/Diffraction/Python/Diffraction_fitting_fp/T3/run-ff.dat")

if __name__ == "__main__":
    main()