import matplotlib.pyplot as plt
import numpy as np
import mol as mol


def graphBumpsResults(filename):
    q, F, dq, Fy = np.loadtxt(filename, skiprows=1).T
    fig, axes = plt.subplots(2, sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    axes[0].plot(q, Fy, zorder=2, label='model')
    axes[0].errorbar(q, F, dq, zorder=1, label='experiment')
    axes[0].legend(fontsize=8)
    axes[0].set_xlim([q[0], q[-1]])
    axes[0].set_ylabel("|F| (1/Å)", fontsize=10)
    axes[0].tick_params(axis="both", which="major", direction="in", labelsize=10)
    axes[1].plot(q, (Fy - F) / dq)
    axes[1].plot(q, q * 0, color="black")
    axes[1].set_ylabel("(f(q) - F) / dq", fontsize=10)
    axes[1].set_xlabel("q (1/Å)", fontsize=10)
    output_filename = filename.split(".")[0] + "-model.png"
    fig.savefig(output_filename)
    plt.draw()


def getSLD(aSL, aArea, dimension, stepsize, normarea=100, bulknsld=0):
    aSLD = np.zeros(dimension)
    for i in range(len(aSL)):
        if aSL[i] == 0 or aArea[i] == 0:
            aSLD[i] = bulknsld
        else:
            aSLD[i] = aSL[i] / (aArea[i] * stepsize) * aArea[
                i
            ] / normarea + bulknsld * (1 - aArea[i] / normarea)
    return aSLD


def getED(aSL, aArea, dimension, stepsize, normarea=100, bulknsld=0):
    return getSLD(aSL, aArea, dimension, stepsize, normarea, bulknsld) * 2.814e-5


def graphBilayer(bilayer, dimension, stepsize, bulknsld, show=False, savefile=None):
    z = np.linspace(0, dimension * stepsize, dimension, endpoint=False)
    dd, aArea, aSL = bilayer.fnWriteProfile(z)
    aSLD = np.zeros(dimension)

    for i in range(dimension):
        if aArea[i] != 0:
            aSLD[i] = aSL[i] / (aArea[i] * stepsize) * aArea[i] / dd + bulknsld * (
                1 - aArea[i] / dd
            )
        else:
            aSLD[i] = bulknsld

    fig, axes = plt.subplots(3)

    axes[0].plot(z, aArea)
    axes[0].legend(["Area"])
    axes[0].set_ylabel("Area (Å^2)")
    axes[0].tick_params(axis="both", which="major")

    axes[1].plot(z, aSL)
    axes[1].legend(["SL"])
    axes[1].set_ylabel("SL (Å)")
    axes[1].tick_params(axis="both", which="major")

    axes[2].plot(z, aSLD)
    axes[2].legend(["SLD"])
    axes[2].set_ylabel("SLD (1/$Å^2$)")
    axes[2].set_xlabel("z (Å)")
    axes[2].tick_params(axis="both", which="major")
    if savefile:
        fig.savefig(savefile + ".png")
    if show:
        plt.show()
    else:
        plt.close()


def graphProfiles(
    obj, dimension, stepsize, maxarea, bulknsld, show=False, savefile=None
):
    z = np.linspace(0, dimension * stepsize, dimension, endpoint=False)
    fig, ax = plt.subplots(3)
    __, aArea, aSL = obj.fnWriteProfile(z)
    if isinstance(obj, mol.BLM_quaternary):
        members = {
            "headgroup1": ["headgroup1", "headgroup2"],
            "lipid": ["lipid1", "lipid2"],
            "methyl": ["methyl1", "methyl2"],
        }
        if obj.nf_lipid_2 > 0:
            members.update(headgroup2=["headgroup1_2", "headgroup2_2"])
        if obj.nf_lipid_3 > 0:
            members.update(headgroup3=["headgroup1_3", "headgroup2_3"])
    elif isinstance(obj, mol.PC):
        members = {
            "carbonyl gylcerol": ["cg"],
            "phosphate": ["phosphate"],
            "choline": ["choline"],
        }
    else:
        members = {}

    ax[0].set_ylabel("Area (Å^2)")
    ax[1].set_ylabel("SL (Å)")
    ax[2].set_ylabel("SLD (1/Å^2)")
    ax[2].set_xlabel("z (Å)")
    ax[0].plot(z, aArea)
    ax[1].plot(z, aSL)
    ax[2].plot(z, getSLD(aSL, aArea, dimension, stepsize, maxarea, bulknsld))
    legend = ["total"]

    for label in members:
        aArea = np.zeros(dimension)
        aSL = np.zeros(dimension)
        aSLD = np.zeros(dimension)
        for name in members[label]:
            group = obj.groups[name]
            __, partial_area, partial_SL = group.fnWriteProfile(z)
            aSL += partial_SL
            aArea += partial_area
        aSLD = getSLD(aSL, aArea, dimension, stepsize, maxarea, bulknsld)
        ax[0].plot(z, aArea)
        ax[1].plot(z, aSL)
        ax[2].plot(z, aSLD)
        legend.append(label)
    fig.legend(legend)

    if savefile:
        fig.savefig(savefile + ".png")
    if show:
        plt.show()
    else:
        plt.close()
