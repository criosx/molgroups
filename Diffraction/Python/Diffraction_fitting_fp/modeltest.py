from bumps.names import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import scipy.fft 
import mol

BLM_quaternary_1 = mol.PyBLM_quaternary()
    na1, nh1, nm1, va1, vm1, vh1, lh1 = 0.00760, 0.00461, 0.000468, 972.00, 98, 331.00, 9.56 
    na2, nh2, nm2, va2, vm2, vh2, lh2 = 0, 0, 0, 0, 0, 0, 0 
    na3, nh3, nm3, va3, vm3, vh3, lh3 = 0, 0, 0, 0, 0, 0, 0
    vc, nc = 0, 0
    BLM_quaternary_1.Init(va1, na1, vm1, nm1, vh1,nh1, lh1, va2, na2, vm2, nm2, 
                            vh2, nh2, lh2, va3, na3, vm3, nm3, vh3, nh3,lh3, vc, nc)

    sigma, bulknsld, startz, l_lipid1, l_lipid2, vf_bilayer = 2.0, 9.4114E-06, 50, 11.6, 11.6, 1
    


def form_factor(bilayer, l_lipid1=11.6, l_lipid2=11.6, vf=1, sigma, bulknsld):
    # initialize canvas, two numpy arrays: area, nSLD; [200] (100Å, 0.5 Å stepsize)
    # initialize with zeros

    BLM_quaternary_1.Set(sigma, bulknsld, startz, l_lipid1, l_lipid2, vf_bilayer)
    
    maxarea = 100
    stepsize = 0.05
    dimension = 3000
    aArea = np.zeros(dimension).tolist()
    anSL  = np.zeros(dimension).tolist()
    anSLD = np.zeros(dimension).tolist()

    # initialize a BLM_quaternary object (standard constructor)
    # set values for lipids using the init function (we do not do this here, we use the standard from constructor)
    # set values using fnSet
    #     sigma = 2.5, bulknsld = 0., startz = 10, l_lipid1 = l_lipid2 = 13.6, vf_bilayer = 1., everything else preset)
    # write the area and nSLD out using fnWriteProfile onto the Canvas (Maxarea=100Å2)

    dd, aArea, anSL = BLM_quaternary_1.WriteProfile(aArea, anSLD, dimension, stepsize, maxarea)

    #fill with bulk
    for i in range(len(aArea)):
        if aArea[i] != 0:
            anSLD[i] = anSL[i] / (aArea[i]*stepsize) * aArea[i]/dd + bulknsld * (1 - aArea[i]/dd)
        else:
            anSLD[i] = bulknsld
    # print(dd)

    aArea = np.asarray(aArea)
    anSLD = np.asarray(anSLD)
    # not for now:

    # identify bilayer center using the bilayer object (methyl1_z+0.5*methyl1_l)
    # center bilayer on canvas
    # |l-----------------------x-------------------------h|
    center = BLM_quaternary_1.GetCenter()
    center = center//stepsize
    canvas_center = dimension//2
    n = int(canvas_center - center)
    centered_bilayer = np.roll(anSLD, n)
    # symmetrize the bilayer around this center
    # |l-----------------------x-------------------------h|          +
    # |h-----------------------x-------------------------l|
    # ------
    # |lh----------------------x------------------------hl|
    symmetrized_bilayer = np.add(centered_bilayer,centered_bilayer[::-1])*0.5
    symmetrized_bilayer -= bulknsld
    half_bilayer = symmetrized_bilayer[int(dimension/2):]

    # cosine tranform of the nSLD(z) canvas to get F(q)
    # (just take scipy.fft.dct, type II)
    dct_dimension = 50000
    F = scipy.fft.dct(half_bilayer, n=dct_dimension)
    F = np.abs(F)
    x = np.array([np.pi / (2*dct_dimension*stepsize)*(2*i+1)  for i in range(int(dct_dimension))])
    # plot F(q)
    x = [ -(50 - center) + .5 * i for i in range(200)]
    plt.plot(x, F)
    plt.show()
    return F
F2 = np.loadtxt("dopc.dat", delimiter=None, skiprows=1)
F2 = np.abs(F2)
# not not for now
# implement F(q) in bumps script with fit parameters such as bilayer thickness, completeness

M1 = Curve(form_factor, )
M2 = 


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
