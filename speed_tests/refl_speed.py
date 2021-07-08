import numpy as np
import time
import ssblm_tiox_both as m
import matplotlib.pyplot as plt

npars = len(m.problem._parameters)
ncalcs = 1000

def gen_randpar(npars, ncalcs, parlist):
    randpar = np.empty((npars, ncalcs))

    for i, par in enumerate(parlist):
        lower, upper = par.to_dict()['bounds']['limits']
        randpar[i,:] = np.random.rand(ncalcs)*(upper - lower) + lower

    return randpar

if 0:
    randpar = gen_randpar(npars, ncalcs, m.problem._parameters)

    starttime = time.time()
    for i in range(ncalcs):
        m.problem.setp(randpar[:,i])
        val = m.problem.chisq_str()

    elapsed = time.time() - starttime

    print('Time elapsed to calculate reflectivity of one model: ', elapsed)

if 0:
    starttime = time.time()
    for _ in range(ncalcs):
        m.bilayer(np.arange(m.dimension)*m.stepsize, m.sigma.value, m.d2o.rho.value, m.global_rough.value, m.tiox.rho.value, m.l_submembrane.value, m.l_lipid1.value, m.l_lipid2.value, m.vf_bilayer.value)
    print('Time elapsed for bilayer calculations: ', time.time()-starttime)

if 0:
    model_step = m.Experiment(sample=m.sample, probe=m.probe, dz=m.zed, step_interfaces = True)
    problem2 = m.FitProblem([model_step])
    
    chisq1 = list()
    chisq2 = list()
    randpar = gen_randpar(npars, ncalcs, m.problem._parameters)
    for i in range(ncalcs):
        m.problem.setp(randpar[:,i])
        problem2.setp(randpar[:,i])
        chisq1.append(float(m.problem.chisq_str().split('(')[0]))
        chisq2.append(float(problem2.chisq_str().split('(')[0]))
        #print(randpar[:,i], chisq1[-1])
    
    plt.hist(chisq1, alpha=0.5)
    plt.hist(chisq2, alpha=0.5)
    plt.show()

        
