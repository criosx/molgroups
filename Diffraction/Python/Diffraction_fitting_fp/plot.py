import matplotlib as plt
import numpy as np
import sys
sys.path.append(".")

def graphFromFile(filename):
    data = np.loadtxt(filename, delimiter=None)
    F = data[:,1]
    q = data[:,2]
    dq = data[:,3]
    fig, ax = plt.subplots()
    ax.figure.set_size_inches(10, 5.5)
    ax.plot(q, F)
    ax.set_ylabel('|F| (1/Å)', fontsize=20)
    ax.set_xlabel('q (1/Å)', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)
    plt.draw()

def main():
    graphFromFile("T1/run-ff.dat")
    graphFromFile("T2/run-ff.dat")

if __name__ == "__main__":
    main()