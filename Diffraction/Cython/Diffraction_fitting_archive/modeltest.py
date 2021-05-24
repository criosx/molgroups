# from bumps.names import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import scipy.fft
import mol

# initialize canvas, two numpy arrays: area, nSLD; [200] (100Å, 0.5 Å stepsize)
# initialize with zeros

aArea = np.zeros(200).tolist()
anSLD = np.zeros(200).tolist()

# initialize a BLM_quaternary object (standard constructor)
# set values for lipids using the init function (we do not do this here, we use the standard from constructor)
# set values using fnSet
#     sigma = 2.5, bulknsld = 0., startz = 10, l_lipid1 = l_lipid2 = 13.6, vf_bilayer = 1., everything else preset)
# write the area and nSLD out using fnWriteProfile onto the Canvas (Maxarea=100Å2)

BLM_quaternary_1 = mol.PyBLM_quaternary()

sigma, bulknsld, startz, l_lipid1, l_lipid2, vf_bilayer = 2.5, 0, 10, 13.6, 13.6, 1
BLM_quaternary_1.Set(sigma, bulknsld, startz, l_lipid1, l_lipid2, vf_bilayer)

dimension, stepsize, maxarea = 200, .5, 100
dd, aArea, anSLD = BLM_quaternary_1.WriteProfile(aArea, anSLD, dimension, stepsize, maxarea)

# using matplotlib, plot the canvas (save the graphs for your presentation)
x = [ .5*i for i in range(200)]
aArea = np.asarray(aArea)
anSLD = np.asarray(anSLD)
plt.plot(x, 100 * aArea)
plt.plot(x, 100 * anSLD)
plt.legend(['Area', 'nSLD'])
plt.show()
# print(aArea)
# not for now:

# identify bilayer center using the bilayer object (methyl1_z+0.5*methyl1_l)
# center bilayer on canvas
# |l-----------------------x-------------------------h|
center = BLM_quaternary_1.GetCenter()
center = round(center/.5) * .5
# print(BLM_quaternary_1.GetLowerLimit())
# print(round(BLM_quaternary_1.GetUpperLimit()))
# print(center)
n = int(2*(50 - center))
centered_bilayer = np.roll(anSLD, n)
# symmetrize the bilayer around this center
# |l-----------------------x-------------------------h|          +
# |h-----------------------x-------------------------l|
# ------
# |lh----------------------x------------------------hl|
symmetrized_bilayer = np.add(centered_bilayer,centered_bilayer[::-1])

# cosine tranform of the nSLD(z) canvas to get F(q)
# (just take scipy.fft.dct, type II)
F = scipy.fft.fft(anSLD)
# plot F(q)
x = [ -(50 - center) + .5 * i for i in range(200)]
plt.plot(x, F)
plt.show()
# not not for now
# implement F(q) in bumps script with fit parameters such as bilayer thickness, completeness






#############################################################################
'''ar1=numpy.loadtxt('500nM.txt')
x1=ar1[:,0]
y1=ar1[:,1]
dy1=[]
for i in range(len(y1)):
    dy1.append(0.1)

ar2=numpy.loadtxt('1uM.txt')
x2=ar2[:,0]
y2=ar2[:,1]
dy2=[]
for i in range(len(y2)):
    dy2.append(0.1)


def fnkinetic_langmuir_dt(lx, t=0, ka=1, kd=0.1):
    #returning dX/dt=[A, B, AB]/dt
    return numpy.array([0, -1*(ka*lx[0]*lx[1]-kd*lx[2]), (ka*lx[0]*lx[1]-kd*lx[2])])


def fnkinetic_langmuir(lt, c0, bmax, ab0, toffset, yoffset, ka, kd):
    # X0=[A0, B0, AB0]
    lt[:] += toffset
    x0 = [c0, bmax-ab0, ab0]
    result = scipy.integrate.odeint(fnkinetic_langmuir_dt, x0, lt, (ka, kd))
    #select only the solution for ab
    result = result[:, 2]
    result[:] += yoffset
    return result


c0_1=0.5
bmax=0
ab0_1=0
toffset_1=0
yoffset_1=0
ka=1
kd=1
M1 = Curve(fnkinetic_langmuir,x1,y1,dy1,c0=c0_1, bmax=bmax, ab0=ab0_1, toffset=toffset_1, yoffset=yoffset_1, ka=ka, kd=kd)
M1.ka.range(0.001,10)
M1.kd.range(0.001,10)
M1.bmax.range(0,80)
M1.ab0.range(0,40)
M1.yoffset.range(440,444)
M1.toffset.range(0,10)

c0_2=1
ab0_2=0
toffset_2=0
yoffset_2=0
M2 = Curve(fnkinetic_langmuir,x2,y2,dy2,c0=c0_2, bmax=bmax, ab0=ab0_2, toffset=toffset_2, yoffset=yoffset_2, ka=ka, kd=kd)
M2.ab0.range(0,40)
M2.yoffset.range(444,449)
M2.toffset.range(0,10)
M2.bmax = M1.bmax
M2.ka = M1.ka
M2.kd = M1.kd

def constraints():
    return (0 if M1.bmax.value >= M1.ab0.value + M2.ab0.value
    and M2.ab0.value > M1.ab0.value
    else 1000+(M1.bmax.value-(M1.ab0.value + M2.ab0.value))**6
    +1000+(M2.ab0.value-M1.ab0.value)**6)


models = M1, M2
problem = MultiFitProblem(models, constraints=constraints)'''

#############################################################################
