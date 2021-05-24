from bumps.names import *
import numpy
import scipy.integrate



ar1=numpy.loadtxt('500nM.txt')
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
problem = MultiFitProblem(models, constraints=constraints)


