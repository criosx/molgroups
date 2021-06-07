import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append(".")

def graphFromFile(filename):
    data = np.loadtxt(filename, delimiter=None)
    q, F, dq, Fy= data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    fig, (ax, ax1) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    # ax.figure.set_size_inches(10, 5.5)
    ax.plot(q, Fy, zorder=2)
    ax.errorbar(q, F, dq, zorder=1)
    ax.legend(['model', 'experiment'], fontsize=8)
    ax.set_xlim([q[0], q[-1]])
    ax.set_ylabel('|F| (1/Å)', fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax1.plot(q, (Fy-F)/dq)
    ax1.plot(q, q*0, color='black')
    # ax.set_ylabel('(f(q) - F) / dq', fontsize=10)
    ax.set_xlabel('q (1/Å)', fontsize=10)
    fig.savefig(savefile(filename))
    plt.draw()

def savefile(filename):
    return filename.split(".")[0]+".pdf"

def main():
    graphFromFile("molgroups/Diffraction/Python/Diffraction_fitting_fp/T1/run-ff.dat")
    graphFromFile("molgroups/Diffraction/Python/Diffraction_fitting_fp/T1_headgroups/run-ff.dat")
    graphFromFile("molgroups/Diffraction/Python/Diffraction_fitting_fp/T2_headgroups/run-ff.dat")
    graphFromFile("molgroups/Diffraction/Python/Diffraction_fitting_fp/T2/run-ff.dat")
    graphFromFile("molgroups/Diffraction/Python/Diffraction_fitting_fp/T3/run-ff.dat")

if __name__ == "__main__":
    main()