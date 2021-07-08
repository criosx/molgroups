import numpy
import math
from scipy.special import erf

class nSLDObj():

    def __init__(self, name=None):
        self.bWrapping=True
        self.bConvolution=False
        self.bProtonExchange=False
        self.dSigmaConvolution=1
        self.iNumberOfConvPoints=7
        self.absorb=0
        
        self.z = 0
        self.l = 0
        self.vol = 0
        self.nSL = 0
        self.nf = 0
        self.sigma = 0

        if name is not None:
            self.name = name
        
    def fnGetAbsorb(self, z):
        return self.absorb
        
    # returns a n-point gaussian interpolation of the area within 4 sigma
    # all area calculatioins are routed through this function, whether they use convolution or not
    # convolution works only for objects with fixed nSLD. Broadening an nSLD profile is not as direct as
    # broadening a nSL profile. For reason, however, objects report a nSLD(z) and not a nSL(z)
    # if it becomes necessary to broaden profiles with variable nSLD, structural changes to the code
    # have to be implemented.
    
    def fnGetArea(self, z):
        raise NotImplementedError()

    def fnGetProfiles(self, z):
        raise NotImplementedError()

    def fnGetConvolutedArea(self, dz):
        # TODO: use scipy image filters or numpy convolution to do this.
        if self.bConvolution:
            dnormsum = 0
            dsum = 0
            for  i in range(self.iNumberOfConvPoints):
                dd = 8 / float(self.iNumberOfConvPoints*i) - 4
                dgauss = numpy.exp((-0.5)*dd*dd)       # (sigma_convolution)^2/(sigma_convolution)^2 cancels
                dnormsum += dgauss
                dsum += self.fnGetArea(dz + dd * self.dSigmaConvolution) * dgauss
            if dnormsum != 0:
                return dsum/dnormsum
            else:
                return 0
        else:
            return self.fnGetArea(dz)
            
    def fnGetLowerLimit(self):
        raise NotImplementedError()
        
    def fnGetnSLD(self, z):
        raise NotImplementedError()
        
    def fnGetUpperLimit(self):
        raise NotImplementedError()

    def fnGetZ(self):
        return self.z
        
    def fnSetConvolution(self, sigma_convolution, iNumberOfConvPoints):
        self.bConvolution = True
        self.dSigmaConvolution = sigma_convolution
        self.iNumberOfConvPoints = iNumberOfConvPoints
        
    def fnWriteData2File(self, f, cName, z):

        header = "z"+cName+" a"+cName+" nsl"+cName
        area, nsl, _ = self.fnGetProfiles(z)
        A = numpy.vstack((z, area, nsl)).T
        numpy.savetxt(f, A, fmt='%0.6e', delimiter=' ', comments='', header=header)
        f.write('\n')
        # TODO: implement wrapping

    def fnWritePar2File(self, fp, cName, z):
        raise NotImplementedError()

    @staticmethod
    def fnWriteConstant(fp, name, darea, dSLD, z):
        header = "Constant " + name + " area " + str(darea) + " \n" + \
                 "z_" + name + " a_" + name + " nsl_" + name
        area = darea * numpy.ones_like(z)
        nsl = dSLD * darea * numpy.gradient(z)
        A = numpy.vstack((z, area, nsl)).T
        numpy.savetxt(fp, A, fmt='%0.6e', delimiter=' ', comments='', header=header)
        fp.write('\n')

    # does a Catmull-Rom Interpolation on an equal distance grid
    # 0<t<=1 is the relative position on the interval between p0 and p1
    # p-1 and p2 are needed for derivative calculation
    def CatmullInterpolate(self, t, pm1, p0, p1, p2):
        m0 = (p1-pm1)/2
        m1 = (p2-p0) /2
        t_2 = t*t
        t_3 = t_2*t
        h00 = 2*t_3-3*t_2+1
        h10 = t_3-2*t_2+t
        h01 = (-2)*t_3+3*t_2
        h11 = t_3-t_2
        return h00*p0+h10*m0+h01*p1+h11*m1
        
    # p is a [4][4][4] array, t is a [3] array
    def fnTriCubicCatmullInterpolate(self, p, t):
        dFirstStage=numpy.zeros(4,4)
        dSecondStage=numpy.zeros(4)

        for i in range(4):
            for j in range(4):
                dFirstStage[i][j] = self.CatmullInterpolate(t[0],p[0][i][j],p[1][i][j],p[2][i][j],p[3][i][j])
        
        for i in range(4):
            dSecondStage[i]=self.CatmullInterpolate(t[1],dFirstStage[0][i],dFirstStage[1][i],dFirstStage[2][i],dFirstStage[3][i])
        
        return self.CatmullInterpolate(t[2],dSecondStage[0],dSecondStage[1],dSecondStage[2],dSecondStage[3])
        
    # p is a [4][4][4][4] array, t is a [4] array
    def fnQuadCubicCatmullInterpolate(self, p, t):
        dFirstStage = numpy.zeros(4,4,4)
        dSecondStage = numpy.zeros(4,4)
        dThirdStage = numpy.zeros(4)

        for i in range(4):
            for j in range(4):
                for k in range(4):
                    dFirstStage[i][j][k] = self.CatmullInterpolate(t[0],p[0][i][j][k],p[1][i][j][k],p[2][i][j][k],p[3][i][j][k])
        
        for i in range(4):
            for j in range(4):
                dSecondStage[i][j] = self.CatmullInterpolate(t[1],dFirstStage[0][i][j],dFirstStage[1][i][j],dFirstStage[2][i][j],dFirstStage[3][i][j])
        
        for i in range(4):
            dThirdStage[i] = self.CatmullInterpolate(t[2],dSecondStage[0][i],dSecondStage[1][i],dSecondStage[2][i],dSecondStage[3][i])
        
        return self.CatmullInterpolate(t[3],dThirdStage[0],dThirdStage[1],dThirdStage[2],dThirdStage[3])

    # Philosophy for this first method: You simply add more and more volume and nSLD to the
    # volume and nSLD array. After all objects have filled up those arrays the maximal area is
    # determined which is the area per molecule and unfilled volume is filled with bulk solvent.
    # Hopefully the fit algorithm finds a physically meaningful solution. There has to be a global
    # hydration paramter for the bilayer.
    # Returns maximum area

    def fnWriteProfile(self, z, aArea=None, anSL=None):
    
        # do we want a 0.5 * stepsize shift? I believe refl1d FunctionalLayer uses 
        #z = numpy.linspace(0, dimension * stepsize, dimension, endpoint=False)
        # z, aArea, anSL must be numpy arrays with the same shape
       
        # TODO implement wrapping
        # TODO implement absorption

        aArea = numpy.zeros_like(z) if aArea is None else aArea
        anSL = numpy.zeros_like(z) if anSL is None else anSL

        assert(aArea.shape == z.shape)
        assert(anSL.shape == z.shape)

        area, nsl, _ = self.fnGetProfiles(z)
        dMaxArea = area.max()
        aArea += area
        anSL += nsl

        return dMaxArea, aArea, anSL

    def fnOverlayProfile(self, z, aArea, anSL, dMaxArea):
        # z, aArea, anSL must be numpy arrays with the same shape
        assert (aArea.shape == z.shape)
        assert (anSL.shape == z.shape)
        
        # TODO implement wrapping
        # TODO implement absorption

        area, nsl, _ = self.fnGetProfiles(z)
        temparea = aArea + area

        # find any area for which the final area will be greater than dMaxArea
        overmax = temparea > dMaxArea
        # note: unphysical overfill will trigger the following assertion error
        # TODO: implement a correction instead of throwing an error
        assert (not numpy.any((temparea - dMaxArea)>aArea))
        anSL[overmax] = anSL[overmax] * (1 - ((temparea[overmax] - dMaxArea[overmax])/aArea[overmax])) + nsl[overmax]
        aArea[overmax] = dMaxArea

        # deal with area for which the final area is not greater than dMaxArea
        aArea[~overmax] += area[~overmax]
        anSL[~overmax] += nsl[~overmax]

        return aArea, anSL

class CompositenSLDObj(nSLDObj):
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def fnFindSubgroups(self):
        # this should be run at the end of init. Could also store keys if memory is an issue and then use self.__dict__[key]
        # Could use a "subgroup adder" function that adds to subgroups
        self.subgroups = [getattr(self, attr) for attr in dir(self) if isinstance(getattr(self, attr), nSLDObj)]

    def fnGetProfiles(self, z):
        area = numpy.zeros_like(z)
        nsl = numpy.zeros_like(z)
        nsld = numpy.zeros_like(z)
        
        # make sure FindSubGroups was called. TODO: speed tests whether this is too expensive to do every time
        if not hasattr(self, 'subgroups'):
            self.fnFindSubgroups()

        for g in self.subgroups:
            newarea, newnsl, _ = g.fnGetProfiles(z)
            area += newarea
            nsl += newnsl
        
        nsld = numpy.zeros_like(z)
        pos = (area > 0)
        nsld[pos] = nsl[pos] / (area[pos] * numpy.gradient(z)[pos])

        return area * self.nf, nsl * self.nf, nsld

    def fnWritePar2File(self, f, cName, z):
        if not hasattr(self, 'subgroups'):
            self.fnFindSubgroups()

        for g in self.subgroups:
            #g.fnWritePar2File(f, cName + '.' + g.name, z) # this allows objects with the same names from multiple bilayers
            g.fnWritePar2File(f, g.name, z) # current behavior

# ------------------------------------------------------------------------------------------------------
# Function Object Implementation
# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Simple Objects
# ------------------------------------------------------------------------------------------------------
class Box2Err(nSLDObj):

    def __init__(self, dz=20, dsigma1=2, dsigma2=2, dlength=10, dvolume=10, dnSL=0, dnumberfraction=1, name=None):
        super().__init__(name=name)
        self.z = dz
        self.sigma1 = dsigma1
        self.sigma2 = dsigma2
        self.l = dlength
        self.vol = dvolume
        self.nSL = dnSL
        self.nf = dnumberfraction
        self.nsldbulk_store = 0.
        self.nSL2 = 0.

    def fnGetProfiles(self, z, bulknsld=0):

        if bulknsld != 0:
            self.nsldbulk_store = bulknsld

        # calculate area
        # Gaussian function definition, integral is volume, return value is area at positions z
        if (self.l != 0) and (self.sigma1 != 0) and (self.sigma2 != 0):
            area = erf((z - self.z + 0.5 * self.l) / (numpy.sqrt(2) * self.sigma1))
            #result = math.erf((dz - self.z + 0.5 * self.l) / math.sqrt(2) / self.sigma1)
            #result -= math.erf((dz - self.z - 0.5 * self.l) / math.sqrt(2) / self.sigma2)
            area -= erf((z - self.z - 0.5 * self.l) / (numpy.sqrt(2) * self.sigma2))            
            area *= (self.vol / self.l) * 0.5
            area *= self.nf
        else:
            area = numpy.zeros_like(z)
        
        # calculate nSLD
        nsld = self.fnGetnSL(bulknsld) / self.vol * numpy.ones_like(z) if self.vol != 0 else numpy.zeros_like(z)

        # calculate nSL.
        nsl = area * nsld * numpy.gradient(z)
        
        return area, nsl, nsld

    def fnGetnSL(self, bulknsld):
        if self.bProtonExchange:
            if self.vol != 0:
                return ((bulknsld + 0.56e-6) * self.nSL2 + (6.36e-6 - bulknsld) * self.nSL) / (6.36e-6 + 0.56e-6)
            else:
                return 0.
        else:
            return self.nSL

    # Gaussians are cut off below and above 3 sigma double
    def fnGetLowerLimit(self):
        return self.z - 0.5 * self.l - 3 * self.sigma1

    def fnGetUpperLimit(self):
        return self.z + 0.5 * self.l + 3 * self.sigma2

    # 7/6/2021 new feature: only use proton exchange if nSL2 is explicitly set!
    def fnSetnSL(self, _nSL, _nSL2=None):
        self.nSL = _nSL
        if _nSL2 is not None:
            self.nSL2 = _nSL2
            self.bProtonExchange = True
        else:
            self.bProtonExchange = False

    def fnSetSigma(self, sigma1, sigma2=None):
        self.sigma1 = sigma1
        self.sigma2 = sigma1 if sigma2 is None else sigma2

    def fnSetZ(self, dz):
        self.z = dz

    def fnWritePar2File(self, fp, cName, z):
        fp.write("Box2Err "+cName+" z "+str(self.z)+" sigma1 "+str(self.sigma1)+" sigma2 "+str(self.sigma2)+" l "
                 + str(self.l)+" vol "+str(self.vol)+" nSL "+str(self.nSL)+" nSL2 "+str(self.nSL2)+" nf "
                 + str(self.nf)+" \n")
        self.fnWriteData2File(fp, cName, z)

# ------------------------------------------------------------------------------------------------------
# Headgroups
# ------------------------------------------------------------------------------------------------------

class PC(CompositenSLDObj):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.cg = Box2Err(name='cg')
        self.phosphate = Box2Err(name='phosphate')
        self.choline = Box2Err(name='choline')

        self.groups = {"cg": self.cg, "phosphate": self.phosphate, "choline": self.choline}

        self.cg.l = 4.21 
        self.cg.sigma1, self.cg.sigma2 = 2.53, 2.29
        self.cg.vol=147 
        self.cg.nSL=3.7755e-4

        self.phosphate.l = 3.86
        self.phosphate.sigma1, self.phosphate.sigma2 = 2.29, 2.02
        self.phosphate.vol=54 
        self.phosphate.nSL=2.8350e-4 

        self.choline.l = 6.34
        self.choline.sigma1, self.choline.sigma2 = 2.02, 2.26
        self.choline.vol=120
        self.choline.nSL=-6.0930e-5

        self.cg.nf=1 
        self.phosphate.nf=1 
        self.choline.nf=1

        self.l = 9.575
        self.vol=self.cg.vol+self.phosphate.vol+self.choline.vol
        self.nSL=self.cg.nSL+self.phosphate.nSL+self.choline.nSL
        self.ph_relative_pos = .5
        self.nf = 1.
        self.fnAdjustParameters()
        self.fnFindSubgroups()

    def fnAdjustParameters(self):
        self.cg.z = self.z - 0.5*self.l + 0.5*self.cg.l
        self.choline.z = self.z + 0.5*self.l - 0.5*self.choline.l
        z0 = self.cg.z + 0.5*self.cg.l
        z1 = self.choline.z - 0.5*self.choline.l
        self.phosphate.z = z0 + (z1 - z0) * self.ph_relative_pos
    
    def fnSet(self, l=9.575, ph_relative_pos=.5, cg_nSL=1.885e-3, ch_nSL=1.407e-3, ph_nSL=1.323e-3):
        self.cg.nSL = cg_nSL
        self.choline.nSL = ch_nSL
        self.phosphate.nSL = ph_nSL
        self.l = l
        self.ph_relative_pos=ph_relative_pos
        self.fnAdjustParameters()
    
    def fnGetLowerLimit(self):
        return self.cg.fnGetLowerLimit()
    
    def fnGetUpperLimit(self):
        return self.choline.fnGetUpperLimit()
    
    def fnGetnSL(self):
        return self.cg.nSL + self.phosphate.nSL + self.choline.nSL
    
    def fnGetZ(self): 
        return self.z
    
    def fnSetSigma(self, sigma):
        self.cg.sigma1=sigma
        self.cg.sigma2=sigma
        self.phosphate.sigma1=sigma
        self.phosphate.sigma2=sigma
        self.choline.sigma1=sigma
        self.choline.sigma2=sigma
    
    def fnSetZ(self, dz):
        self.z = dz
        self.fnAdjustParameters()
    
    def fnWritePar2File (self, fp, cName, z):
        fp.write("PC "+cName+" z "+str(self.z)+" l "+str(self.l)+" vol " +
                 str(self.cg.vol + self.phosphate.vol + self.choline.vol)+" nf " + str(self.nf)+" \n")
        self.fnWriteData2File(fp, cName, z)

class PCm(PC):
    def __init__ (self, **kwargs):
        super().__init__(**kwargs)
        self.cg.sigma2=2.53
        self.cg.sigma1=2.29
        self.phosphate.sigma2=2.29
        self.phosphate.sigma1=2.02
        self.choline.sigma2=2.02
        self.choline.sigma1=2.26
        self.fnAdjustParameters()
        
    def fnAdjustParameters(self):
        self.cg.z = self.z + 0.5*self.l-0.5*self.cg.l
        self.choline.z = self.z - 0.5*self.l+0.5*self.choline.l
        z0 = self.choline.z + 0.5*self.choline.l
        z1 = self.cg.z - 0.5*self.cg.l
        self.phosphate.z = z0 + (z1 - z0) * (1-self.ph_relative_pos)

    def fnGetLowerLimit(self):
        return self.choline.fnGetLowerLimit()

    def fnGetUpperLimit(self): 
        return self.cg.fnGetUpperLimit()

    def fnWritePar2File (self,fp, cName, z):
        fp.write("PCm " + cName + " z " + str(self.z) + " l " + str(self.l) + " vol " +
                 str(self.cg.vol + self.phosphate.vol + self.choline.vol) + " nf " + str(self.nf) + " \n")
        self.fnWriteData2File(fp, cName, z)


# ------------------------------------------------------------------------------------------------------
# Lipid Bilayer
# ------------------------------------------------------------------------------------------------------

class BLM_quaternary(CompositenSLDObj):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.headgroup1 = PCm(name='headgroup1') #Box2Err() 
        self.lipid1 = Box2Err(name='lipid1')
        self.methyl1 = Box2Err(name='methyl1')
        self.methyl2 = Box2Err(name='methyl2')
        self.lipid2 = Box2Err(name='lipid2')
        self.headgroup2 = PC(name='headgroup2') #Box2Err() #PC()                          # PC head group
        self.headgroup1_2 = Box2Err(name='headgroup1_2')                        # second headgroups
        self.headgroup2_2 = Box2Err(name='headgroup2_2')
        self.headgroup1_3 = Box2Err(name='headgroup1_3')
        self.headgroup2_3 = Box2Err(name='headgroup2_3')
       
        self.defect_hydrocarbon = Box2Err(name='defect_hc')
        self.defect_headgroup = Box2Err(name='defect_hg')
        
        self.groups = {"headgroup1": self.headgroup1, "headgroup1_2": self.headgroup1_2,
                        "headgroup1_3": self.headgroup1_3, "lipid1": self.lipid1,
                        "methyl1": self.methyl1,"methyl2": self.methyl2,
                        "lipid2": self.lipid2, "headgroup2": self.headgroup2,
                        "headgroup2_2": self.headgroup2_2, "headgroup2_3": self.headgroup2_3,
                        "defect_hc": self.defect_hydrocarbon, "defect_hg": self.defect_headgroup}

        self.volacyllipid = 925
        self.nslacyllipid = -2.67e-4
        self.volmethyllipid = 98.8
        self.nslmethyllipid = -9.15e-5
        self.volacyllipid_2 = 925
        self.nslacyllipid_2 = -2.67e-4
        self.volmethyllipid_2 = 98.8
        self.nslmethyllipid_2 = -9.15e-5
        self.volacyllipid_3 = 925
        self.nslacyllipid_3 = -2.67e-4
        self.volmethyllipid_3 = 98.8
        self.nslmethyllipid_3 = -9.15e-5
        
        self.volchol = 630
        self.nslchol = 1.3215e-4

        self.headgroup1.vol = 330
        self.headgroup1_2.vol = 330
        self.headgroup1_3.vol = 330
        self.headgroup2.vol = 330
        self.headgroup2_2.vol = 330
        self.headgroup2_3.vol = 330
        self.headgroup1.nSL = 6.0012e-4
        self.headgroup1_2.nSL = 6.0012e-4
        self.headgroup1_3.nSL = 6.0012e-4
        self.headgroup2.nSL = 6.0012e-4
        self.headgroup2_2.nSL = 6.0012e-4
        self.headgroup2_3.nSL = 6.0012e-4
        self.headgroup1.l = 9.5
        self.headgroup1_2.l = 9.5
        self.headgroup1_3.l = 9.5
        self.headgroup2.l = 9.5
        self.headgroup2_2.l = 9.5
        self.headgroup2_3.l = 9.5
        
        self.nf_lipid_2 = 0.
        self.nf_lipid_3 = 0.
        self.nf_chol = 0
        self.nf = 1.

        self.vf_bilayer = 1.0
        self.absorb = 0.
        self.l_lipid1 = 11.
        self.l_lipid2 = 11.
        self.bulknsld = -0.56e-6
        self.normarea = 60.
        self.startz = 50.
        self.sigma = 2.
        self.methyl_sigma = 2.
        self.radius_defect = 100.
        self.hc_substitution_1 = 0
        self.hc_substitution_2 = 0

        self.fnAdjustParameters()
        self.fnFindSubgroups()
        
    def fnInit(self, va1, na1, vm1, nm1, vh1, nh1, lh1, va2=0, na2=0, vm2=0,
               nm2=0, vh2=0, nh2=0, lh2=0, va3=0, na3=0, vm3=0, nm3=0, vh3=0,
               nh3=0, lh3=0, vc=0, nc=0):
        self.volacyllipid = va1
        self.nslacyllipid = na1
        self.volmethyllipid = vm1
        self.nslmethyllipid = nm1
        self.volacyllipid_2 = va2
        self.nslacyllipid_2 = na2
        self.volmethyllipid_2 = vm2
        self.nslmethyllipid_2 = nm2
        self.volacyllipid_3 = va3
        self.nslacyllipid_3 = na3
        self.volmethyllipid_3 = vm3
        self.nslmethyllipid_3 = nm3
        
        self.volchol = vc
        self.nslchol = nc

        self.headgroup1.vol = vh1
        self.headgroup1_2.vol = vh2
        self.headgroup1_3.vol = vh3
        self.headgroup2.vol = vh1
        self.headgroup2_2.vol = vh2
        self.headgroup2_3.vol = vh3
        self.headgroup1.nSL = nh1
        self.headgroup1_2.nSL = nh2
        self.headgroup1_3.nSL = nh3
        self.headgroup2.nSL = nh1
        self.headgroup2_2.nSL = nh2
        self.headgroup2_3.nSL = nh3
        self.headgroup1.l = lh1
        self.headgroup1_2.l = lh2
        self.headgroup1_3.l = lh3
        self.headgroup2.l = lh1
        self.headgroup2_2.l = lh2
        self.headgroup2_3.l = lh3
        
        self.fnAdjustParameters()

    def fnAdjustParameters(self):
        # printf("Enter AdjustParameters \n")
        self.fnSetSigma(self.sigma)
        
        self.l_lipid1 = max(self.l_lipid1, 0.01)
        self.l_lipid2 = max(self.l_lipid2, 0.01)
        self.nf_lipid_2 = max(self.nf_lipid_2, 0)
        self.nf_lipid_3 = max(self.nf_lipid_3, 0)
        self.nf_chol = max(self.nf_chol, 0)
        sum = self.nf_lipid_2 + self.nf_lipid_3 + self.nf_chol
        if sum > 1: #modified so these are all divided by the same sum
            self.nf_lipid_2 = self.nf_lipid_2 / sum
            self.nf_lipid_3 = self.nf_lipid_3 / sum
            self.nf_chol = self.nf_chol / sum
        self.vf_bilayer = max(self.vf_bilayer, 1E-5)
        self.vf_bilayer = min(self.vf_bilayer, 1)
        
        # outer hydrocarbons
        l_ohc = self.l_lipid2
        nf_ohc_lipid  = 1 - self.nf_lipid_2 - self.nf_lipid_3 - self.nf_chol
        nf_ohc_lipid_2 = self.nf_lipid_2
        nf_ohc_lipid_3 = self.nf_lipid_3
        nf_ohc_chol = self.nf_chol
        V_ohc = nf_ohc_lipid * (self.volacyllipid - self.volmethyllipid) + nf_ohc_lipid_2 * (self.volacyllipid_2 - self.volmethyllipid_2) + nf_ohc_lipid_3 * (self.volacyllipid_3 - self.volmethyllipid_3) + nf_ohc_chol * self.volchol
        nSL_ohc = nf_ohc_lipid * (self.nslacyllipid - self.nslmethyllipid) + nf_ohc_lipid_2 * (self.nslacyllipid_2 - self.nslmethyllipid_2) + nf_ohc_lipid_3 * (self.nslacyllipid_3 - self.nslmethyllipid_3) + nf_ohc_chol * self.nslchol
        
        self.normarea = V_ohc / l_ohc
        c_s_ohc = self.vf_bilayer
        c_A_ohc = 1
        c_V_ohc = 1
        
        # printf("ssBLM: normarea %lf \n",normarea)
        self.lipid2.l = l_ohc
        self.lipid2.vol = V_ohc
        self.lipid2.nSL = nSL_ohc
        self.lipid2.nf = c_s_ohc * c_A_ohc * c_V_ohc
        # printf("c: c_s_ohc %lf c_A_ohc %lf c_V_ohc %lf \n", c_s_ohc, c_A_ohc, c_V_ohc)
        
        # outer methyl
        nf_om_lipid = nf_ohc_lipid
        nf_om_lipid_2 = nf_ohc_lipid_2
        nf_om_lipid_3 = nf_ohc_lipid_3
        V_om = nf_om_lipid * self.volmethyllipid + nf_om_lipid_2 * self.volmethyllipid_2 + nf_om_lipid_3 * self.volmethyllipid_3
        l_om = l_ohc * V_om / V_ohc
        nSL_om = nf_om_lipid * self.nslmethyllipid + nf_om_lipid_2 * self.nslmethyllipid_2 + nf_om_lipid_3 * self.nslmethyllipid_3
        
        c_s_om = c_s_ohc
        c_A_om = 1
        c_V_om = 1
        
        self.methyl2.l = l_om
        self.methyl2.vol = V_om
        self.methyl2.nSL = nSL_om
        self.methyl2.nf = c_s_om * c_A_om * c_V_om
        
        # inner hydrocarbons
        l_ihc = self.l_lipid1
        
        nf_ihc_lipid = nf_ohc_lipid
        nf_ihc_lipid_2 = nf_ohc_lipid_2
        nf_ihc_lipid_3 = nf_ohc_lipid_3
        nf_ihc_chol = nf_ohc_chol
        V_ihc = nf_ihc_lipid * (self.volacyllipid - self.volmethyllipid) + nf_ihc_lipid_2 * (self.volacyllipid_2 - self.volmethyllipid_2) + nf_ihc_lipid_3 * (self.volacyllipid_3 - self.volmethyllipid_3) + nf_ihc_chol * self.volchol
        nSL_ihc = nf_ihc_lipid * (self.nslacyllipid - self.nslmethyllipid) + nf_ihc_lipid_2 * (self.nslacyllipid_2 - self.nslmethyllipid_2) + nf_ihc_lipid_3 * (self.nslacyllipid_3 - self.nslmethyllipid_3) + nf_ihc_chol * self.nslchol
        
        c_s_ihc = self.vf_bilayer
        c_A_ihc = self.normarea * l_ihc / V_ihc
        c_V_ihc = 1
        
        self.lipid1.l = l_ihc
        self.lipid1.vol = V_ihc
        self.lipid1.nSL = nSL_ihc
        self.lipid1.nf = c_s_ihc * c_A_ihc * c_V_ihc
        
        # inner methyl
        nf_im_lipid = nf_ihc_lipid
        nf_im_lipid_2 = nf_ihc_lipid_2
        nf_im_lipid_3 = nf_ihc_lipid_3
        V_im = nf_im_lipid * self.volmethyllipid + nf_im_lipid_2 * self.volmethyllipid_2 + nf_im_lipid_3 * self.volmethyllipid_3
        l_im = l_ihc * V_im / V_ihc
        nSL_im = nf_im_lipid * self.nslmethyllipid + nf_im_lipid_2 * self.nslmethyllipid_2 + nf_im_lipid_3 * self.nslmethyllipid_3
        
        c_s_im = c_s_ihc
        c_A_im = c_A_ihc
        c_V_im = 1
        
        self.methyl1.l = l_im
        self.methyl1.vol = V_im
        self.methyl1.nSL = nSL_im
        self.methyl1.nf = c_s_im * c_A_im * c_V_im
        
        # outer headgroups
        self.headgroup2.nf = c_s_ohc * c_A_ohc * nf_ohc_lipid * (1 - self.hc_substitution_2)
        self.headgroup2_2.nf = c_s_ohc * c_A_ohc * nf_ohc_lipid_2 * (1 - self.hc_substitution_2)
        self.headgroup2_3.nf = c_s_ohc * c_A_ohc * nf_ohc_lipid_3 * (1 - self.hc_substitution_2)
        
        # inner headgroups
        self.headgroup1.nf   = c_s_ihc * c_A_ihc * nf_ihc_lipid * (1 - self.hc_substitution_2)
        self.headgroup1_2.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid_2 * (1 - self.hc_substitution_2)
        self.headgroup1_3.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid_3 * (1 - self.hc_substitution_2)
        # printf("c: c_s_ihc %lf c_A_ihc %lf nf_ihc_lipid %lf hc_substitution_1 %lf \n", c_s_ihc, c_A_ihc, nf_ihc_lipid, hc_substitution_1)
        
        self.lipid1.z= self.startz + self.headgroup1.l + 0.5 * self.lipid1.l
        self.headgroup1.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1.l)
        self.headgroup1.fnAdjustParameters()
        self.headgroup1_2.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1_2.l)
        self.headgroup1_3.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1_3.l)
        self.methyl1.z = self.lipid1.z + 0.5 * (self.lipid1.l + self.methyl1.l)
        self.methyl2.z = self.methyl1.z + 0.5 * (self.methyl1.l + self.methyl2.l)
        self.lipid2.z = self.methyl2.z + 0.5 * (self.methyl2.l + self.lipid2.l)
        self.headgroup2.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2.l)
        self.headgroup2.fnAdjustParameters()
        self.headgroup2_2.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2_2.l)
        self.headgroup2_3.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2_3.l)
        # printf("nf bme %lf tether %lf tetherg %lf lipid1 %lf headgroup1 %lf headgroup1_2 %lf headgroup1_3 %lf methyl1 %lf methyl2 %lf lipid2 %lf headgroup2 %lf headgroup2_2 %lf headgroup2_3 %lf \n", bME.nf, tether.nf, tetherg.nf, lipid1.nf, headgroup1.nf, headgroup1_2.nf, headgroup1_3.nf, methyl1.nf, methyl2.nf, lipid2.nf, headgroup2.nf, headgroup2_2.nf, headgroup2_3.nf)
        
        # defects
        hclength = self.lipid1.l + self.methyl1.l + self.methyl2.l + self.lipid2.l
        hglength = self.headgroup1.l + self.headgroup2.l
        
        if self.radius_defect<(0.5*(hclength+hglength)):
            self.radius_defect = 0.5 * (hclength+hglength)
        # printf("defect_radius %lf hclength %lf \n",radius_defect, hclength)
        
        volhalftorus = 3.14159265359 * 3.14159265359 * (self.radius_defect - (2 * hclength / 3 / 3.14159265359)) * hclength * hclength / 4
        volcylinder = 3.14159265359 * self.radius_defect * self.radius_defect * hclength
        # printf("volhalftorus %lf volcylinder %lf \n", volhalftorus, volcylinder)
        defectarea = volhalftorus / volcylinder * (1 - self.vf_bilayer) * self.normarea
        # printf("defectarea %lf \n", defectarea)
        
        self.defect_hydrocarbon.vol = defectarea * hclength
        self.defect_hydrocarbon.l = hclength
        self.defect_hydrocarbon.z = self.lipid1.z - 0.5 * self.lipid1.l + 0.5 * hclength
        self.defect_hydrocarbon.nSL = self.lipid2.nSL / self.lipid2.vol * self.defect_hydrocarbon.vol
        self.defect_hydrocarbon.fnSetSigma(self.sigma)
        self.defect_hydrocarbon.nf = 1
        
        defectratio = self.defect_hydrocarbon.vol / self.lipid2.vol
        self.defect_headgroup.vol = defectratio * (self.headgroup2.vol * self.headgroup2.nf + self.headgroup2_2.vol * self.headgroup2_2.nf + self.headgroup2_3.vol * self.headgroup2_3.nf)
        self.defect_headgroup.l = hclength + hglength
        self.defect_headgroup.z = self.headgroup1.fnGetZ() - 0.5 * self.headgroup1.l + 0.5 * (hclength + hglength)
        self.defect_headgroup.nSL = defectratio * (self.headgroup2.fnGetnSL()*self.headgroup2.nf * self.headgroup2.nf + self.headgroup2_2.fnGetnSL(self.bulknsld) * self.headgroup2_2.nf + self.headgroup2_3.fnGetnSL(self.bulknsld) * self.headgroup2_3.nf)
        self.defect_headgroup.fnSetSigma(self.sigma)
        self.defect_headgroup.nf = 1

    # return value of center of the membrane
    def fnGetCenter(self):
        return self.methyl1.z + 0.5 * self.methyl1.l

    # Use limits of molecular subgroups
    def fnGetLowerLimit(self):
        a = self.headgroup1.fnGetLowerLimit()
        b = self.headgroup1_2.fnGetLowerLimit()
        c = self.headgroup1_3.fnGetLowerLimit()
        return min([a, b, c])
        
    def fnGetUpperLimit(self):
        a = self.headgroup2.fnGetUpperLimit()
        b = self.headgroup2_2.fnGetUpperLimit()
        c = self.headgroup2_3.fnGetUpperLimit()
        return max([a, b, c])

    def fnSet(self, sigma, bulknsld, startz, l_lipid1, l_lipid2, vf_bilayer,
              nf_lipid_2=0., nf_lipid_3=0., nf_chol=0., hc_substitution_1=0.,
              hc_substitution_2=0., radius_defect=100.):
        self.sigma = sigma
        self.bulknsld = bulknsld
        self.startz = startz
        self.l_lipid1 = l_lipid1
        self.l_lipid2 = l_lipid2
        self.vf_bilayer = vf_bilayer
        self.nf_lipid_2 = nf_lipid_2
        self.nf_lipid_3 = nf_lipid_3
        self.nf_chol = nf_chol
        self.hc_substitution_1 = hc_substitution_1
        self.hc_substitution_2 = hc_substitution_2
        self.radius_defect = radius_defect
        
        self.fnAdjustParameters()

    def fnSetSigma(self, sigma):
        self.headgroup1.fnSetSigma(sigma)
        self.headgroup1_2.fnSetSigma(sigma)
        self.headgroup1_3.fnSetSigma(sigma)
        self.lipid1.fnSetSigma(sigma,sigma+self.methyl_sigma)
        self.methyl1.fnSetSigma(sigma+self.methyl_sigma, sigma+self.methyl_sigma)
        self.methyl2.fnSetSigma(sigma+self.methyl_sigma, sigma+self.methyl_sigma)
        self.lipid2.fnSetSigma(sigma+self.methyl_sigma,sigma)
        self.headgroup2.fnSetSigma(sigma)
        self.headgroup2_2.fnSetSigma(sigma)
        self.headgroup2_3.fnSetSigma(sigma)
        self.defect_hydrocarbon.fnSetSigma(sigma)
        self.defect_headgroup.fnSetSigma(sigma)

    def fnWritePar2File(self, fp, cName, z):
        super().fnWritePar2File(fp, cName, z)
        self.fnWriteConstant(fp, "normarea", self.normarea, 0, z)

class child_ssBLM_quaternary(BLM_quaternary):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.substrate  = Box2Err(name='substrate')
        self.siox       = Box2Err(name='siox')
        self.groups["substrate"].append(self.substrate)
        self.groups["siox"].append(self.siox)
        self.substrate.l=20
        self.substrate.z=10
        self.substrate.nf=1
        self.substrate.sigma1=2.0

        self.siox.l=20
        self.siox.z=30
        self.siox.nf=1
        self.siox.fnSetSigma(2.0)
        
        self.bulknsld=-0.56e-6

    def fnAdjustParameters(self):
        super().fnAdjustParameters()
        self.substrate.vol=self.normarea*self.substrate.l
        self.substrate.nSL=self.rho_substrate*self.substrate.vol
        self.siox.l=self.l_siox
        self.siox.vol=self.normarea*self.siox.l
        self.siox.nSL=self.rho_siox*self.siox.vol
        self.siox.z=self.substrate.l+0.5*self.siox.l
    
    def fnGetLowerLimit(self):
        return self.substrate.fnGetLowerLimit()

    def fnSet(self, sigma, global_rough, rho_substrate, bulknsld, rho_siox, l_siox,
         l_submembrane, l_lipid1, l_lipid2, vf_bilayer, nf_lipid_2=0, nf_lipid3=0,
              nf_chol=0, hc_substitution_1=0, hc_substitution_2=0, radius_defect=100):
        self.global_rough= global_rough
        self.rho_substrate= rho_substrate
        self.bulknsld= bulknsld
        self.rho_siox= rho_siox
        self.l_siox= l_siox
        self.l_submembrane= l_submembrane
        super().fnSet(sigma, bulknsld, self.startz, l_lipid1, l_lipid2, vf_bilayer)

    def fnSetSigma(self, sigma):
        super().fnSetSigma(sigma)
        self.substrate.sigma2 = self.global_rough
        self.siox.sigma1 = self.global_rough
        self.siox.sigma2 = self.global_rough

    def fnWritePar2File (self, fp, cName, dimension, stepsize):
        pass
    
    def fnWriteProfile(self, aArea, anSLD, dimension, stepsize, dMaxArea):
        pass


class ssBLM_quaternary(CompositenSLDObj):
    def __init__(self, **kwargs):

        super().__init__(**kwargs)
        self.substrate = Box2Err(name='substrate')
        self.siox = Box2Err(name='siox')
        self.headgroup1 = PCm(name='headgroup1')
        self.lipid1 = Box2Err(name='lipid1')
        self.methyl1 = Box2Err(name='methyl1')
        self.methyl2 = Box2Err(name='methyl2')
        self.lipid2 = Box2Err(name='lipid2')
        self.headgroup2 = PC(name='headgroup2')  # PC head group
        self.headgroup1_2 = Box2Err(name='headgroup1_2')  # second headgroups
        self.headgroup2_2 = Box2Err(name='headgrou2_2')
        self.headgroup1_3 = Box2Err(name='headgroup1_3')
        self.headgroup2_3 = Box2Err(name='headgroup2_3')

        self.defect_hydrocarbon = Box2Err(name='defect_hc')
        self.defect_headgroup = Box2Err(name='defect_hg')

        self.substrate.l = 40
        self.substrate.z = 0
        self.substrate.nf = 1
        self.substrate.sigma1 = 2.0
        self.rho_substrate = 1
        self.l_siox = 1
        self.rho_siox = 1

        self.siox.l = 20
        self.siox.z = 30
        self.siox.nf = 1
        self.siox.fnSetSigma(2.0)

        self.volacyllipid = 925
        self.nslacyllipid = -2.67e-4
        self.volmethyllipid = 98.8
        self.nslmethyllipid = -9.15e-5

        self.volacyllipid_2 = 925
        self.nslacyllipid_2 = -2.67e-4
        self.volmethyllipid_2 = 98.8
        self.nslmethyllipid_2 = -9.15e-5
        self.volacyllipid_3 = 925
        self.nslacyllipid_3 = -2.67e-4
        self.volmethyllipid_3 = 98.8
        self.nslmethyllipid_3 = -9.15e-5

        self.volchol = 630
        self.nslchol = 1.3215e-4

        self.headgroup1.vol = 330
        self.headgroup1_2.vol = 330
        self.headgroup1_3.vol = 330
        self.headgroup2.vol = 330
        self.headgroup2_2.vol = 330
        self.headgroup2_3.vol = 330
        self.headgroup1.nSL = 6.0012e-4
        self.headgroup1_2.nSL = 6.0012e-4
        self.headgroup1_3.nSL = 6.0012e-4
        self.headgroup2.nSL = 6.0012e-4
        self.headgroup2_2.nSL = 6.0012e-4
        self.headgroup2_3.nSL = 6.0012e-4
        self.headgroup1.l = 9.5
        self.headgroup1_2.l = 9.5
        self.headgroup1_3.l = 9.5
        self.headgroup2.l = 9.5
        self.headgroup2_2.l = 9.5
        self.headgroup2_3.l = 9.5
        self.l_submembrane = 10.

        self.hc_substitution_1 = 0
        self.hc_substitution_2 = 0

        self.headgroup1_2.vol = 330  # was 330
        self.headgroup2_2.vol = 330  # was 330
        self.headgroup1_3.vol = 330  # was 330
        self.headgroup2_3.vol = 330  # was 330
        self.headgroup1_2.nSL = 6.0012e-4  # was 6.0122e-4
        self.headgroup2_2.nSL = 6.0012e-4  # was 6.0122e-4
        self.headgroup1_3.nSL = 6.0012e-4  # was 6.0122e-4
        self.headgroup2_3.nSL = 6.0012e-4  # was 6.0122e-4
        self.headgroup1_2.l = 9.5
        self.headgroup2_2.l = 9.5
        self.headgroup1_3.l = 9.5
        self.headgroup2_3.l = 9.5

        self.volacyllipid_2 = 925
        self.nslacyllipid_2 = -2.67e-4
        self.volmethyllipid_2 = 98.8
        self.nslmethyllipid_2 = -9.15e-5

        self.volacyllipid_3 = 925
        self.nslacyllipid_3 = -2.67e-4
        self.volmethyllipid_3 = 98.8
        self.nslmethyllipid_3 = -9.15e-5

        self.volchol = 630
        self.nslchol = 1.3215e-4
        self.nf_lipid_2 = 0  # for preparing towards a general bilayer class
        self.nf_lipid_3 = 0
        self.nf_chol = 0
        self.bulknsld = -0.56e-6

        self.nf_lipid_2 = 0.
        self.nf_lipid_3 = 0.
        self.nf_chol = 0.
        self.nf = 1.

        self.vf_bilayer = 1.0
        self.absorb = 0.
        self.l_lipid1 = 11.
        self.l_lipid2 = 11.
        self.bulknsld = -0.56e-6
        self.normarea = 60.
        self.startz = 50.
        self.sigma = 2.
        self.methyl_sigma = 2.
        self.radius_defect = 100.
        self.hc_substitution_1 = 0
        self.hc_substitution_2 = 0
        self.global_rough = 2.0

        self.fnAdjustParameters()
        self.fnFindSubgroups()

    def fnInit(self, va1, na1, vm1, nm1, vh1, nh1, lh1, va2, na2, vm2, nm2, vh2, nh2, lh2, va3, na3, vm3, nm3, vh3,
               nh3, lh3, vc, nc):
        self.volacyllipid = va1
        self.nslacyllipid = na1
        self.volmethyllipid = vm1
        self.nslmethyllipid = nm1
        self.volacyllipid_2 = va2
        self.nslacyllipid_2 = na2
        self.volmethyllipid_2 = vm2
        self.nslmethyllipid_2 = nm2
        self.volacyllipid_3 = va3
        self.nslacyllipid_3 = na3
        self.volmethyllipid_3 = vm3
        self.nslmethyllipid_3 = nm3

        self.volchol = vc
        self.nslchol = nc

        self.headgroup1.vol = vh1
        self.headgroup1_2.vol = vh2
        self.headgroup1_3.vol = vh3
        self.headgroup2.vol = vh1
        self.headgroup2_2.vol = vh2
        self.headgroup2_3.vol = vh3
        self.headgroup1.nSL = nh1
        self.headgroup1_2.nSL = nh2
        self.headgroup1_3.nSL = nh3
        self.headgroup2.nSL = nh1
        self.headgroup2_2.nSL = nh2
        self.headgroup2_3.nSL = nh3
        self.headgroup1.l = lh1
        self.headgroup1_2.l = lh2
        self.headgroup1_3.l = lh3
        self.headgroup2.l = lh1
        self.headgroup2_2.l = lh2
        self.headgroup2_3.l = lh3

        self.fnAdjustParameters()

    def fnAdjustParameters(self):
        # printf("Enter AdjustParameters \n")
        self.fnSetSigma(self.sigma)

        if self.l_lipid1 <= 0:
            self.l_lipid1 = 0.01
        if self.l_lipid2 <= 0:
            self.l_lipid2 = 0.01
        if self.nf_lipid_2 < 0:
            self.nf_lipid_2 = 0
        if self.nf_lipid_3 < 0:
            self.nf_lipid_3 = 0
        if self.nf_chol < 0:
            self.nf_chol = 0
        if self.nf_lipid_2 + self.nf_lipid_3 + self.nf_chol > 1:
            self.nf_lipid_2 = self.nf_lipid_2 / (self.nf_lipid_2 + self.nf_lipid_3 + self.nf_chol)
            self.nf_lipid_3 = self.nf_lipid_3 / (self.nf_lipid_2 + self.nf_lipid_3 + self.nf_chol)
            self.nf_chol = self.nf_chol / (self.nf_lipid_2 + self.nf_lipid_3 + self.nf_chol)
        if self.vf_bilayer <= 0:
            self.vf_bilayer = 1e-5
        if self.vf_bilayer > 1:
            self.vf_bilayer = 1

        # def ssBLM_quaternary::fnAdjustParameters()
        # def fnAdjustParameters(self):
        # printf("Enter AdjustParameters \n")

        self.substrate.fnSetSigma(self.global_rough)
        self.siox.fnSetSigma(self.global_rough)
        self.fnSetSigma(self.sigma)

        # outer hydrocarbons
        l_ohc = self.l_lipid2
        nf_ohc_lipid = 1 - self.nf_lipid_2 - self.nf_lipid_3 - self.nf_chol
        nf_ohc_lipid_2 = self.nf_lipid_2
        nf_ohc_lipid_3 = self.nf_lipid_3
        nf_ohc_chol = self.nf_chol
        V_ohc = nf_ohc_lipid * (self.volacyllipid - self.volmethyllipid) + nf_ohc_lipid_2 * (
                    self.volacyllipid_2 - self.volmethyllipid_2) + nf_ohc_lipid_3 * (
                            self.volacyllipid_3 - self.volmethyllipid_3) + nf_ohc_chol * self.volchol
        nSL_ohc = nf_ohc_lipid * (self.nslacyllipid - self.nslmethyllipid) + nf_ohc_lipid_2 * (
                    self.nslacyllipid_2 - self.nslmethyllipid_2) + nf_ohc_lipid_3 * (
                              self.nslacyllipid_3 - self.nslmethyllipid_3) + nf_ohc_chol * self.nslchol

        self.normarea = V_ohc / l_ohc
        c_s_ohc = self.vf_bilayer
        c_A_ohc = 1
        c_V_ohc = 1

        # printf("ssBLM: normarea %lf \n",normarea)
        self.lipid2.l = l_ohc
        self.lipid2.vol = V_ohc
        self.lipid2.nSL = nSL_ohc
        self.lipid2.nf = c_s_ohc * c_A_ohc * c_V_ohc
        # printf("c: c_s_ohc %lf c_A_ohc %lf c_V_ohc %lf \n", c_s_ohc, c_A_ohc, c_V_ohc)

        # outer methyl
        nf_om_lipid = nf_ohc_lipid
        nf_om_lipid_2 = nf_ohc_lipid_2
        nf_om_lipid_3 = nf_ohc_lipid_3
        V_om = nf_om_lipid * self.volmethyllipid + nf_om_lipid_2 * self.volmethyllipid_2 + nf_om_lipid_3 * self.volmethyllipid_3
        l_om = l_ohc * V_om / V_ohc
        nSL_om = nf_om_lipid * self.nslmethyllipid + nf_om_lipid_2 * self.nslmethyllipid_2 + nf_om_lipid_3 * self.nslmethyllipid_3

        c_s_om = c_s_ohc
        c_A_om = 1
        c_V_om = 1

        self.methyl2.l = l_om
        self.methyl2.vol = V_om
        self.methyl2.nSL = nSL_om
        self.methyl2.nf = c_s_om * c_A_om * c_V_om
        # inner hydrocarbons
        l_ihc = self.l_lipid1

        nf_ihc_lipid = nf_ohc_lipid
        nf_ihc_lipid_2 = nf_ohc_lipid_2
        nf_ihc_lipid_3 = nf_ohc_lipid_3
        nf_ihc_chol = nf_ohc_chol
        V_ihc = nf_ihc_lipid * (self.volacyllipid - self.volmethyllipid) + nf_ihc_lipid_2 * (
                    self.volacyllipid_2 - self.volmethyllipid_2) + nf_ihc_lipid_3 * (
                            self.volacyllipid_3 - self.volmethyllipid_3) + nf_ihc_chol * self.volchol
        nSL_ihc = nf_ihc_lipid * (self.nslacyllipid - self.nslmethyllipid) + nf_ihc_lipid_2 * (
                    self.nslacyllipid_2 - self.nslmethyllipid_2) + nf_ihc_lipid_3 * (
                              self.nslacyllipid_3 - self.nslmethyllipid_3) + nf_ihc_chol * self.nslchol

        c_s_ihc = self.vf_bilayer
        c_A_ihc = self.normarea * l_ihc / V_ihc
        c_V_ihc = 1

        self.lipid1.l = l_ihc
        self.lipid1.vol = V_ihc
        self.lipid1.nSL = nSL_ihc
        self.lipid1.nf = c_s_ihc * c_A_ihc * c_V_ihc

        # inner methyl
        nf_im_lipid = nf_ihc_lipid
        nf_im_lipid_2 = nf_ihc_lipid_2
        nf_im_lipid_3 = nf_ihc_lipid_3
        V_im = nf_im_lipid * self.volmethyllipid + nf_im_lipid_2 * self.volmethyllipid_2 + nf_im_lipid_3 * self.volmethyllipid_3
        l_im = l_ihc * V_im / V_ihc
        nSL_im = nf_im_lipid * self.nslmethyllipid + nf_im_lipid_2 * self.nslmethyllipid_2 + nf_im_lipid_3 * self.nslmethyllipid_3

        c_s_im = c_s_ihc
        c_A_im = c_A_ihc
        c_V_im = 1

        self.methyl1.l = l_im
        self.methyl1.vol = V_im
        self.methyl1.nSL = nSL_im
        self.methyl1.nf = c_s_im * c_A_im * c_V_im

        # outer headgroups
        self.headgroup2.nf = c_s_ohc * c_A_ohc * nf_ohc_lipid * (1 - self.hc_substitution_2)
        self.headgroup2_2.nf = c_s_ohc * c_A_ohc * nf_ohc_lipid_2 * (1 - self.hc_substitution_2)
        self.headgroup2_3.nf = c_s_ohc * c_A_ohc * nf_ohc_lipid_3 * (1 - self.hc_substitution_2)

        # inner headgroups
        self.headgroup1.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid * (1 - self.hc_substitution_2)
        self.headgroup1_2.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid_2 * (1 - self.hc_substitution_2)
        self.headgroup1_3.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid_3 * (1 - self.hc_substitution_2)
        # printf("c: c_s_ihc %lf c_A_ihc %lf nf_ihc_lipid %lf hc_substitution_1 %lf \n", c_s_ihc, c_A_ihc, nf_ihc_lipid, hc_substitution_1)

        # substrate
        self.substrate.vol = self.normarea * self.substrate.l
        self.substrate.nSL = self.rho_substrate * self.substrate.vol
        self.siox.l = self.l_siox
        self.siox.vol = self.normarea * self.siox.l
        self.siox.nSL = self.rho_siox * self.siox.vol

        # set all lengths
        self.siox.z=self.substrate.l / 2  +0.5*self.siox.l
        self.lipid1.z= self.substrate.l / 2 + self.siox.l + self.l_submembrane + self.headgroup1.l + 0.5 * self.lipid1.l
        self.headgroup1.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1.l)
        self.headgroup1_2.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1_2.l)
        self.headgroup1_3.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1_3.l)
        self.methyl1.z = self.lipid1.z + 0.5 * (self.lipid1.l + self.methyl1.l)
        self.methyl2.z = self.methyl1.z + 0.5 * (self.methyl1.l + self.methyl2.l)
        self.lipid2.z = self.methyl2.z + 0.5 * (self.methyl2.l + self.lipid2.l)
        self.headgroup2.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2.l)
        self.headgroup2_2.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2_2.l)
        self.headgroup2_3.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2_3.l)
        # printf("nf bme %lf tether %lf tetherg %lf lipid1 %lf headgroup1 %lf headgroup1_2 %lf headgroup1_3 %lf methyl1 %lf methyl2 %lf lipid2 %lf headgroup2 %lf headgroup2_2 %lf headgroup2_3 %lf \n", bME.nf, tether.nf, tetherg.nf, lipid1.nf, headgroup1.nf, headgroup1_2.nf, headgroup1_3.nf, methyl1.nf, methyl2.nf, lipid2.nf, headgroup2.nf, headgroup2_2.nf, headgroup2_3.nf)

        # defects
        hclength = self.lipid1.l + self.methyl1.l + self.methyl2.l + self.lipid2.l
        hglength = self.headgroup1.l + self.headgroup2.l

        if self.radius_defect < (0.5 * (hclength + hglength)):
            self.radius_defect = 0.5 * (hclength + hglength)
        # printf("defect_radius %lf hclength %lf \n",radius_defect, hclength)

        volhalftorus = 3.14159265359 * 3.14159265359 * (
                    self.radius_defect - (2 * hclength / 3 / 3.14159265359)) * hclength * hclength / 4
        volcylinder = 3.14159265359 * self.radius_defect * self.radius_defect * hclength
        # printf("volhalftorus %lf volcylinder %lf \n", volhalftorus, volcylinder)
        defectarea = volhalftorus / volcylinder * (1 - self.vf_bilayer) * self.normarea
        # printf("defectarea %lf \n", defectarea)

        self.defect_hydrocarbon.vol = defectarea * hclength
        self.defect_hydrocarbon.l = hclength
        self.defect_hydrocarbon.z = self.lipid1.z - 0.5 * self.lipid1.l + 0.5 * hclength
        self.defect_hydrocarbon.nSL = self.lipid2.nSL / self.lipid2.vol * self.defect_hydrocarbon.vol
        self.defect_hydrocarbon.fnSetSigma(self.sigma)
        self.defect_hydrocarbon.nf = 1

        defectratio = self.defect_hydrocarbon.vol / self.lipid2.vol
        self.defect_headgroup.vol = defectratio * (
                    self.headgroup2.vol * self.headgroup2.nf + self.headgroup2_2.vol * self.headgroup2_2.nf + self.headgroup2_3.vol * self.headgroup2_3.nf)
        self.defect_headgroup.l = hclength + hglength
        self.defect_headgroup.z = self.headgroup1.fnGetZ() - 0.5 * self.headgroup1.l + 0.5 * (hclength + hglength)
        self.defect_headgroup.nSL = defectratio * (
                    self.headgroup2.fnGetnSL() * self.headgroup2.nf + self.headgroup2_2.fnGetnSL( self.bulknsld) *
                    self.headgroup2_2.nf + self.headgroup2_3.fnGetnSL(self.bulknsld) * self.headgroup2_3.nf)
        self.defect_headgroup.fnSetSigma(self.sigma)
        self.defect_headgroup.nf = 1

    # printf("Exit AdjustParameters \n")

    # Use limits of molecular subgroups
    def fnGetLowerLimit(self):
        return self.substrate.fnGetLowerLimit()

    def fnGetUpperLimit(self):
        a = self.headgroup2.fnGetUpperLimit()
        b = self.headgroup2_2.fnGetUpperLimit()
        c = self.headgroup2_3.fnGetUpperLimit()
        return max([a, b, c])

    # void ssBLM_quaternary::fnSet(double _sigma, double _global_rough, double _rho_substrate, double _bulknsld, double _rho_siox, double _l_siox, double _l_submembrane,  double _l_lipid1, double _l_lipid2, double _vf_bilayer, double _nf_lipid_2, double _nf_lipid_3, double _nf_chol, double _hc_substitution_1, double _hc_substitution_2, double _radius_defect){

    # printf("Enter fnSet \n")

    def fnSet(self, _sigma, _bulknsld, _global_rough, _rho_substrate, _rho_siox, _l_siox, _l_submembrane, _l_lipid1,
              _l_lipid2, _vf_bilayer, _nf_lipid_2=0., _nf_lipid_3=0., _nf_chol=0., _hc_substitution_1=0.,
              _hc_substitution_2=0., _radius_defect=100.):
        self.sigma = _sigma
        self.bulknsld = _bulknsld
        self.global_rough = _global_rough
        self.rho_substrate = _rho_substrate
        self.rho_siox = _rho_siox
        self.l_siox = _l_siox  # error undefined variable
        self.l_submembrane = _l_submembrane  # error undefined variabe
        self.l_lipid1 = _l_lipid1
        self.l_lipid2 = _l_lipid2
        self.vf_bilayer = _vf_bilayer
        self.nf_lipid_2 = _nf_lipid_2
        self.nf_lipid_3 = _nf_lipid_3
        self.nf_chol = _nf_chol
        self.hc_substitution_1 = _hc_substitution_1
        self.hc_substitution_2 = _hc_substitution_2
        self.radius_defect = _radius_defect

        self.fnAdjustParameters()

    # printf("Exit fnSet \n")

    def fnSetSigma(self, sigma):
        self.substrate.sigma2 = self.global_rough
        self.siox.sigma1 = self.global_rough
        self.siox.sigma2 = self.global_rough
        self.headgroup1.fnSetSigma(sigma)
        self.headgroup1_2.fnSetSigma(sigma)
        self.headgroup1_3.fnSetSigma(sigma)
        self.lipid1.fnSetSigma(sigma, sigma + self.methyl_sigma)
        self.methyl1.fnSetSigma(sigma + self.methyl_sigma, sigma + self.methyl_sigma)
        self.methyl2.fnSetSigma(sigma + self.methyl_sigma, sigma + self.methyl_sigma)
        self.lipid2.fnSetSigma(sigma + self.methyl_sigma, sigma)
        self.headgroup2.fnSetSigma(sigma)
        self.headgroup2_2.fnSetSigma(sigma)
        self.headgroup2_3.fnSetSigma(sigma)
        self.defect_hydrocarbon.fnSetSigma(sigma)
        self.defect_headgroup.fnSetSigma(sigma)

    def fnWritePar2File(self, fp, cName, z):
        super().fnWritePar2File(fp, cName, z)
        self.fnWriteConstant(fp, "normarea", self.normarea, 0, z)


# ------------------------------------------------------------------------------------------------------
# Tethered Lipid bilayer - binary system
# ------------------------------------------------------------------------------------------------------
class tBLM_quaternary(CompositenSLDObj):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.substrate = Box2Err(name='substrate')
        self.bME = Box2Err(name='bME')
        self.tether = Box2Err(name='tether')
        self.tetherg = Box2Err(name='tetherg')
        self.headgroup1 = PCm(name='headgroup1')  # mirrored PC head group
        self.lipid1 = Box2Err(name='lipid1')
        self.methyl1 = Box2Err(name='methyl1')
        self.methyl2 = Box2Err(name='methyl2')
        self.lipid2 = Box2Err(name='lipid2')
        self.headgroup2 = PC(name='headgroup2')  # PC head group
        self.headgroup1_2 = Box2Err(name='headgroup1_2')  # second headgroups
        self.headgroup2_2 = Box2Err(name='headgroup2_2')
        self.headgroup1_3 = Box2Err(name='headgroup1_3')
        self.headgroup2_3 = Box2Err(name='headgroup2_3')

        self.defect_hydrocarbon = Box2Err(name='defect_hc')
        self.defect_headgroup = Box2Err(name='defect_hg')

        self.substrate.l = 40
        self.substrate.z = 0
        self.substrate.nf = 1
        self.substrate.sigma1 = 2.0

        self.bME.vol = 110
        self.bME.nSL = 3.243e-5
        self.bME.l = 5.2
        self.tether.vol = 380
        self.tether.nSL = 2.1864e-4
        self.tetherg.vol = 110
        self.tetherg.nSL = 1.8654e-4

        self.volacyllipid = 925
        self.nslacyllipid = -2.67e-4
        self.volmethyllipid = 98.8
        self.nslmethyllipid = -9.15e-5
        self.volmethyltether = 98.8
        self.nslmethyltether = -9.15e-5
        self.volacyltether = 982
        self.nslacyltether = -2.85e-4

        self.hc_substitution_1 = 0
        self.hc_substitution_2 = 0
        self.headgroup1_2.vol = 330  # was 330
        self.headgroup2_2.vol = 330  # was 330
        self.headgroup1_3.vol = 330  # was 330
        self.headgroup2_3.vol = 330  # was 330
        self.headgroup1_2.nSL = 6.0012e-4  # was 6.0122e-4
        self.headgroup2_2.nSL = 6.0012e-4  # was 6.0122e-4
        self.headgroup1_3.nSL = 6.0012e-4  # was 6.0122e-4
        self.headgroup2_3.nSL = 6.0012e-4  # was 6.0122e-4
        self.headgroup1_2.l = 9.5
        self.headgroup2_2.l = 9.5
        self.headgroup1_3.l = 9.5
        self.headgroup2_3.l = 9.5

        self.volacyllipid_2 = 925
        self.nslacyllipid_2 = -2.67e-4
        self.volmethyllipid_2 = 98.8
        self.nslmethyllipid_2 = -9.15e-5

        self.volacyllipid_3 = 925
        self.nslacyllipid_3 = -2.67e-4
        self.volmethyllipid_3 = 98.8
        self.nslmethyllipid_3 = -9.15e-5

        self.volchol = 630
        self.nslchol = 1.3215e-4

        self.nf_lipid_2 = 0.
        self.nf_lipid_3 = 0.
        self.nf_chol = 0.
        self.nf = 1.

        self.vf_bilayer = 1.0
        self.absorb = 0.
        self.l_lipid1 = 11.
        self.l_lipid2 = 11.
        self.bulknsld = -0.56e-6
        self.normarea = 60.
        self.sigma = 2.
        self.methyl_sigma = 2.
        self.radius_defect = 100.
        self.hc_substitution_1 = 0
        self.hc_substitution_2 = 0
        self.global_rough = 2.0
        self.nf_tether = 0.6
        # self.nf_ihc_tether = 0.6
        self.l_tether = 3
        self.mult_tether = 5
        # self.nf_ihc_tether = 0.3
        self.substrate.nf = 1
        self.substrate.sigma1 = 2.0
        self.rho_substrate = 1

        self.fnAdjustParameters()
        self.fnFindSubgroups()

    def fnInit(self, va1, na1, vm1, nm1, vh1, nh1, lh1, va2, na2, vm2, nm2, vh2, nh2, lh2, va3, na3, vm3, nm3, vh3, nh3,
               lh3, vc, nc):
        self.volacyllipid = va1
        self.nslacyllipid = na1
        self.volmethyllipid = vm1
        self.nslmethyllipid = nm1
        self.volacyllipid_2 = va2
        self.nslacyllipid_2 = na2
        self.volmethyllipid_2 = vm2
        self.nslmethyllipid_2 = nm2
        self.volacyllipid_3 = va3
        self.nslacyllipid_3 = na3
        self.volmethyllipid_3 = vm3
        self.nslmethyllipid_3 = nm3

        self.volchol = vc
        self.nslchol = nc

        self.headgroup1.vol = vh1
        self.headgroup1_2.vol = vh2
        self.headgroup1_3.vol = vh3
        self.headgroup2.vol = vh1
        self.headgroup2_2.vol = vh2
        self.headgroup2_3.vol = vh3
        self.headgroup1.nSL = nh1
        self.headgroup1_2.nSL = nh2
        self.headgroup1_3.nSL = nh3
        self.headgroup2.nSL = nh1
        self.headgroup2_2.nSL = nh2
        self.headgroup2_3.nSL = nh3
        self.headgroup1.l = lh1
        self.headgroup1_2.l = lh2
        self.headgroup1_3.l = lh3
        self.headgroup2.l = lh1
        self.headgroup2_2.l = lh2
        self.headgroup2_3.l = lh3

        self.fnAdjustParameters()

    def fnAdjustParameters(self):
        # printf("Enter AdjustParameters \n")
        self.fnSetSigma(self.sigma)

        if self.l_lipid1 <= 0:
            self.l_lipid1 = 0.01
        if self.l_lipid2 <= 0:
            self.l_lipid2 = 0.01
        if self.l_tether <= 0:
            self.l_tether = 0.01
        if self.nf_lipid_2 < 0:
            self.nf_lipid_2 = 0
        if self.nf_lipid_3 < 0:
            self.nf_lipid_3 = 0
        if self.nf_chol < 0:
            self.nf_chol = 0
        if self.nf_lipid_2 + self.nf_lipid_3 + self.nf_chol > 1:
            self.nf_lipid_2 = self.nf_lipid_2 / (self.nf_lipid_2 + self.nf_lipid_3 + self.nf_chol)
            self.nf_lipid_3 = self.nf_lipid_3 / (self.nf_lipid_2 + self.nf_lipid_3 + self.nf_chol)
            self.nf_chol = self.nf_chol / (self.nf_lipid_2 + self.nf_lipid_3 + self.nf_chol)
        if self.vf_bilayer <= 0:
            self.vf_bilayer = 1e-5
        if self.vf_bilayer > 1:
            self.vf_bilayer = 1

        # outer hydrocarbons
        l_ohc = self.l_lipid2
        nf_ohc_lipid = 1 - self.nf_lipid_2 - self.nf_lipid_3 - self.nf_chol
        nf_ohc_lipid_2 = self.nf_lipid_2
        nf_ohc_lipid_3 = self.nf_lipid_3
        nf_ohc_chol = self.nf_chol
        V_ohc = nf_ohc_lipid * (self.volacyllipid - self.volmethyllipid) + nf_ohc_lipid_2 * (
                    self.volacyllipid_2 - self.volmethyllipid_2) + nf_ohc_lipid_3 * (
                            self.volacyllipid_3 - self.volmethyllipid_3) + nf_ohc_chol * self.volchol
        nSL_ohc = nf_ohc_lipid * (self.nslacyllipid - self.nslmethyllipid) + nf_ohc_lipid_2 * (
                    self.nslacyllipid_2 - self.nslmethyllipid_2) + nf_ohc_lipid_3 * (
                              self.nslacyllipid_3 - self.nslmethyllipid_3) + nf_ohc_chol * self.nslchol

        self.normarea = V_ohc / l_ohc
        c_s_ohc = self.vf_bilayer
        c_A_ohc = 1
        c_V_ohc = 1

        # printf("ssBLM: normarea %lf \n",normarea)

        self.lipid2.l = l_ohc
        self.lipid2.vol = V_ohc
        self.lipid2.nSL = nSL_ohc
        self.lipid2.nf = c_s_ohc * c_A_ohc * c_V_ohc
        # printf("c: c_s_ohc %lf c_A_ohc %lf c_V_ohc %lf \n", c_s_ohc, c_A_ohc, c_V_ohc)

        # outer methyl
        nf_om_lipid = nf_ohc_lipid
        nf_om_lipid_2 = nf_ohc_lipid_2
        nf_om_lipid_3 = nf_ohc_lipid_3
        V_om = nf_om_lipid * self.volmethyllipid + nf_om_lipid_2 * self.volmethyllipid_2 + nf_om_lipid_3 * self.volmethyllipid_3
        l_om = l_ohc * V_om / V_ohc
        nSL_om = nf_om_lipid * self.nslmethyllipid + nf_om_lipid_2 * self.nslmethyllipid_2 + nf_om_lipid_3 * self.nslmethyllipid_3

        c_s_om = c_s_ohc
        c_A_om = 1
        c_V_om = 1

        self.methyl2.l = l_om
        self.methyl2.vol = V_om
        self.methyl2.nSL = nSL_om
        self.methyl2.nf = c_s_om * c_A_om * c_V_om

        # inner hydrocarbons
        l_ihc = self.l_lipid1
        nf_ihc_tether = self.nf_tether
        nf_ihc_lipid = (1 - nf_ihc_tether) * nf_ohc_lipid
        nf_ihc_lipid_2 = nf_ohc_lipid_2
        nf_ihc_lipid_3 = nf_ohc_lipid_3
        nf_ihc_chol = nf_ohc_chol

        nf_ihc_lipid_2 = (1 - nf_ihc_tether) * nf_ohc_lipid_2
        nf_ihc_lipid_3 = (1 - nf_ihc_tether) * nf_ohc_lipid_3
        nf_ihc_chol = (1 - nf_ihc_tether) * nf_ohc_chol

        V_ihc = nf_ihc_lipid * (self.volacyllipid - self.volmethyllipid) + nf_ihc_lipid_2 * (
                    self.volacyllipid_2 - self.volmethyllipid_2) + nf_ihc_lipid_3 * (
                            self.volacyllipid_3 - self.volmethyllipid_3) + nf_ihc_chol * self.volchol + nf_ihc_tether * (
                            self.volacyltether - self.volmethyltether)
        nSL_ihc = nf_ihc_lipid * (self.nslacyllipid - self.nslmethyllipid) + nf_ihc_lipid_2 * (
                    self.nslacyllipid_2 - self.nslmethyllipid_2) + nf_ihc_lipid_3 * (
                              self.nslacyllipid_3 - self.nslmethyllipid_3) + nf_ihc_chol * self.nslchol + nf_ihc_tether * (
                              self.nslacyltether - self.nslmethyltether)

        c_s_ihc = self.vf_bilayer
        c_A_ihc = self.normarea * l_ihc / V_ihc
        c_V_ihc = 1

        self.lipid1.l = l_ihc
        self.lipid1.vol = V_ihc
        self.lipid1.nSL = nSL_ihc
        self.lipid1.nf = c_s_ihc * c_A_ihc * c_V_ihc

        # inner methyl -- add tether
        nf_im_lipid = nf_ihc_lipid
        nf_im_lipid_2 = nf_ihc_lipid_2
        nf_im_lipid_3 = nf_ihc_lipid_3
        nf_im_tether = nf_ihc_tether
        V_im = nf_im_lipid * self.volmethyllipid + nf_im_lipid_2 * self.volmethyllipid_2 + nf_im_lipid_3 * self.volmethyllipid_3 + nf_im_tether * self.volmethyltether
        l_im = l_ihc * V_im / V_ihc
        nSL_im = nf_im_lipid * self.nslmethyllipid + nf_im_lipid_2 * self.nslmethyllipid_2 + nf_im_lipid_3 * self.nslmethyllipid_3 + nf_im_tether * self.nslmethyltether

        c_s_im = c_s_ihc
        c_A_im = c_A_ihc
        c_V_im = 1

        self.methyl1.l = l_im
        self.methyl1.vol = V_im
        self.methyl1.nSL = nSL_im
        self.methyl1.nf = c_s_im * c_A_im * c_V_im

        # outer headgroups
        self.headgroup2.nf = c_s_ohc * c_A_ohc * nf_ohc_lipid * (1 - self.hc_substitution_2)
        self.headgroup2_2.nf = c_s_ohc * c_A_ohc * nf_ohc_lipid_2 * (1 - self.hc_substitution_2)
        self.headgroup2_3.nf = c_s_ohc * c_A_ohc * nf_ohc_lipid_3 * (1 - self.hc_substitution_2)

        # inner headgroups
        self.headgroup1.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid * (1 - self.hc_substitution_2)
        self.headgroup1_2.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid_2 * (1 - self.hc_substitution_2)
        self.headgroup1_3.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid_3 * (1 - self.hc_substitution_2)

        # tether glycerol part -- transliteration here
        V_tg = self.tetherg.vol

        c_s_tg = c_s_ihc
        c_A_tg = c_A_ihc
        c_V_tg = nf_ihc_tether * (1 - self.hc_substitution_2)

        self.tetherg.l = self.tetherg.vol / ((self.volacyltether - self.volmethyltether) / self.lipid1.l) / 0.9
        self.tetherg.nf = c_s_tg * c_A_tg * c_V_tg

        # tether EO part
        l_EO = self.l_tether
        V_EO = self.tether.vol

        c_s_EO = c_s_ihc
        c_A_EO = c_A_ihc
        c_V_EO = nf_ihc_tether * (1 - self.hc_substitution_2)
        self.tether.nf = c_s_EO * c_A_EO * c_V_EO
        self.tether.l = l_EO

        if (self.tether.nf * self.tether.vol / self.tether.l) > self.normarea:
            self.tether.l = (self.tether.nf * self.tether.vol) / self.normarea
        self.l_tether = self.tether.l

        # bME
        self.bME.l = 5.2
        l_bME = self.bME.l
        self.headgroup1.l = 9.575
        V_bME = self.bME.vol

        d1 = self.headgroup1.l + self.bME.l - self.tether.l - self.tetherg.l
        if (d1 > 0):
            self.bME.l = self.bME.l - d1 / 2
            self.headgroup1.l = self.headgroup1.l - d1 / 2

        if ((
                self.tether.nf * self.tether.vol / self.tether.l + self.mult_tether * self.tether.nf * self.bME.vol / self.bME.l) > self.normarea):
            # print(self.tether.nf, self.tether.vol, self.tether.l, self.mult_tether, self.bME.vol, self.bME.l)
            self.mult_tether = ((self.normarea - self.tether.nf * self.tether.vol / self.tether.l) / (
                        self.bME.vol / self.bME.l)) / self.tether.nf
            if (self.mult_tether < 0):
                self.mult_tether = 0

        self.bME.nf = self.tether.nf * self.mult_tether  # 2.333

        # substrate
        self.substrate.vol = self.normarea * self.substrate.l
        self.substrate.nSL = self.rho_substrate * self.substrate.vol

        # set all lengths
        self.bME.z = 0.5 * self.bME.l + self.substrate.l * 0.5
        self.tether.z = 0.5 * self.tether.l + self.substrate.l * 0.5
        self.tetherg.z = self.tether.z + 0.5 * self.tether.l + 0.5 * self.tetherg.l
        self.lipid1.z = self.tetherg.z + 0.5 * (self.tetherg.l + self.lipid1.l)
        self.headgroup1.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1.l)
        self.headgroup1_2.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1_2.l)
        self.headgroup1_3.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1_3.l)
        self.methyl1.z = self.lipid1.z + 0.5 * (self.lipid1.l + self.methyl1.l)
        self.methyl2.z = self.methyl1.z + 0.5 * (self.methyl1.l + self.methyl2.l)
        self.lipid2.z = self.methyl2.z + 0.5 * (self.methyl2.l + self.lipid2.l)
        self.headgroup2.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2.l)
        self.headgroup2_2.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2_2.l)
        self.headgroup2_3.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2_3.l)

        # printf("nf bme %lf tether %lf tetherg %lf lipid1 %lf headgroup1 %lf headgroup1_2 %lf headgroup1_3 %lf methyl1 %lf methyl2 %lf lipid2 %lf headgroup2 %lf headgroup2_2 %lf headgroup2_3 %lf \n", bME.nf, tether.nf, tetherg.nf, lipid1.nf, headgroup1.nf, headgroup1_2.nf, headgroup1_3.nf, methyl1.nf, methyl2.nf, lipid2.nf, headgroup2.nf, headgroup2_2.nf, headgroup2_3.nf)
        # defects
        hclength = self.lipid1.l + self.methyl1.l + self.methyl2.l + self.lipid2.l
        hglength = self.headgroup1.l + self.headgroup2.l

        if self.radius_defect < (0.5 * (hclength + hglength)):
            self.radius_defect = 0.5 * (hclength + hglength)
        # printf("defect_radius %lf hclength %lf \n",radius_defect, hclength)

        volhalftorus = 3.14159265359 * 3.14159265359 * (
                    self.radius_defect - (2 * hclength / 3 / 3.14159265359)) * hclength * hclength / 4
        volcylinder = 3.14159265359 * self.radius_defect * self.radius_defect * hclength
        # printf("volhalftorus %lf volcylinder %lf \n", volhalftorus, volcylinder)
        defectarea = volhalftorus / volcylinder * (1 - self.vf_bilayer) * self.normarea
        # printf("defectarea %lf \n", defectarea)

        self.defect_hydrocarbon.vol = defectarea * hclength
        self.defect_hydrocarbon.l = hclength
        self.defect_hydrocarbon.z = self.lipid1.z - 0.5 * self.lipid1.l + 0.5 * hclength
        self.defect_hydrocarbon.nSL = self.lipid2.nSL / self.lipid2.vol * self.defect_hydrocarbon.vol
        self.defect_hydrocarbon.fnSetSigma(self.sigma)
        self.defect_hydrocarbon.nf = 1

        defectratio = self.defect_hydrocarbon.vol / self.lipid2.vol
        self.defect_headgroup.vol = defectratio * (
                    self.headgroup2.vol * self.headgroup2.nf + self.headgroup2_2.vol * self.headgroup2_2.nf + self.headgroup2_3.vol * self.headgroup2_3.nf)
        self.defect_headgroup.l = hclength + hglength
        self.defect_headgroup.z = self.headgroup1.fnGetZ() - 0.5 * self.headgroup1.l + 0.5 * (hclength + hglength)
        self.defect_headgroup.nSL = defectratio * (
                    self.headgroup2.fnGetnSL() * self.headgroup2.nf + self.headgroup2_2.fnGetnSL(
                self.bulknsld) * self.headgroup2_2.nf + self.headgroup2_3.fnGetnSL(
                self.bulknsld) * self.headgroup2_3.nf)
        self.defect_headgroup.fnSetSigma(self.sigma)
        self.defect_headgroup.nf = 1

    # Use limits of molecular subgroups
    def fnGetLowerLimit(self):
        return self.substrate.fnGetLowerLimit()

    def fnGetUpperLimit(self):
        a = self.headgroup2.fnGetUpperLimit()
        b = self.headgroup2_2.fnGetUpperLimit()
        c = self.headgroup2_3.fnGetUpperLimit()
        return max([a, b, c])

    def fnSet(self, _sigma, _bulknsld, _global_rough, _rho_substrate, _nf_tether, _mult_tether, _l_tether, _l_lipid1,
              _l_lipid2, _vf_bilayer, _nf_lipid_2=0., _nf_lipid_3=0., _nf_chol=0., _hc_substitution_1=0.,
              _hc_substitution_2=0., _radius_defect=100.):
        self.sigma = _sigma
        self.bulknsld = _bulknsld
        self.global_rough = _global_rough
        self.rho_substrate = _rho_substrate
        self.nf_tether = _nf_tether
        self.mult_tether = _mult_tether
        self.l_tether = _l_tether
        self.l_lipid1 = _l_lipid1
        self.l_lipid2 = _l_lipid2
        self.vf_bilayer = _vf_bilayer
        self.nf_lipid_2 = _nf_lipid_2
        self.nf_lipid_3 = _nf_lipid_3
        self.nf_chol = _nf_chol
        self.hc_substitution_1 = _hc_substitution_1
        self.hc_substitution_2 = _hc_substitution_2
        self.radius_defect = _radius_defect

        self.fnAdjustParameters()

    def fnSetSigma(self, sigma):
        self.substrate.sigma2 = self.global_rough
        self.bME.sigma1 = self.global_rough
        self.bME.sigma2 = self.global_rough
        self.tether.sigma1 = self.global_rough
        self.tether.sigma2 = sigma
        self.tetherg.fnSetSigma(sigma)
        self.headgroup1.fnSetSigma(sigma)
        self.headgroup1_2.fnSetSigma(sigma)
        self.headgroup1_3.fnSetSigma(sigma)
        self.lipid1.fnSetSigma(sigma, sigma + self.methyl_sigma)
        self.methyl1.fnSetSigma(sigma + self.methyl_sigma, sigma + self.methyl_sigma)
        self.methyl2.fnSetSigma(sigma + self.methyl_sigma, sigma + self.methyl_sigma)
        self.lipid2.fnSetSigma(sigma + self.methyl_sigma, sigma)
        self.headgroup2.fnSetSigma(sigma)
        self.headgroup2_2.fnSetSigma(sigma)
        self.headgroup2_3.fnSetSigma(sigma)
        self.defect_hydrocarbon.fnSetSigma(sigma)
        self.defect_headgroup.fnSetSigma(sigma)

    def fnWritePar2File(self, fp, cName, z):
        super().fnWritePar2File(fp, cName, z)
        self.fnWriteConstant(fp, "normarea", self.normarea, 0, z)

# ------------------------------------------------------------------------------------------------------
# Hermite Spline
# ------------------------------------------------------------------------------------------------------

class Hermite(nSLDObj):
    def __init__(self, n, dstartposition, dnSLD, dnormarea, **kwargs):
        super().__init__(**kwargs)
        self.numberofcontrolpoints = n
        self.nSLD=dnSLD
        self.normarea=dnormarea
        self.monotonic=1
        self.damping=1
        self.dampthreshold=0.001
        self.dampFWHM=0.0002
        self.damptrigger=0.04
        
        self.dp     = [0.]*n
        self.vf     = [0.]*n
        self.damp   = [0.]*n
        
    def fnGetArea(self, dz): 
        temp = self.fnGetSplineArea(dz, self.dp, self.vf, self.damping)*self.normarea*self.nf
        return max(temp, 0)

    def fnGetnSLD(self, dz): 
        return self.nSLD

    def fnGetLowerLimit(self): 
        return self.dp[0]

    def fnGetUpperLimit(self): 
        return self.dp[-1]

    def fnGetVolume(self, dz1, dz2): 
        return self.fnGetSplineIntegral(dz1, dz2, self.dp, self.vf, self.damping)*self.normarea*self.nf

    def fnSetNormarea(self, dnormarea): 
        self.normarea = dnormarea

    def fnSetnSLD(self, dnSLD): 
        self.nSLD = dnSLD

    def fnSetRelative(self, dSpacing, dStart, dDp, dVf, dnf): 
        for i in range(self.numberofcontrolpoints):
            self.vf[i] = max(0, dVf[i])
            self.dp[i] = dStart + dSpacing*i + dDp[i]
        self.nf = dnf

    def fnGetSplineAntiDerivative(self, dz, dp, dh): 
        interval, m0, m1, p0, p1 = self.fnGetSplinePars(dz, dp, dh, 0, 0, 0, 0)
        if (0 <= interval < self.numberofcontrolpoints-1):
            dd=dp[interval+1]-dp[interval] 
            t=(dz-dp[interval])/dd 
            t_2=t*t 
            t_3=t_2*t 
            t_4=t_3*t 
            h00= (1/2)*t_4-t_3+t 
            h01=(-1/2)*t_4+t_3             
            h10= (1/4)*t_4-(2/3)*t_3+(1/2)*t_2   
            h11= (1/4)*t_4-(1/3)*t_3             
            return dd*(h00*p0 + h10*dd*m0 + h01*p1 + h11*dd*m1)         
        return 0

    def fnGetSplineArea(self, dz, dp, dh, damping):
        m0, m1, p0, p1 = 0, 0, 0, 0
        if (damping == 0):
            self.damp = dh.copy()  
        else:
            dampfactor = 1 
            for i in range (self.numberofcontrolpoints):
                    self.damp[i] = dh[i] * dampfactor 
                    if (dh[i] >= self.damptrigger):
                        dampfactor=dampfactor*(1/(1 + math.exp(-2.1*(dh[i]- self.dampthreshold)/self.dampFWHM))) 

        interval, m0, m1, p0, p1 = self.fnGetSplinePars(dz, dp, self.damp, m0, m1, p0, p1) 

        if (0 <= interval < self.numberofcontrolpoints-1):
            dd=dp[interval+1]-dp[interval] 
            t=(dz-dp[interval])/dd 
            t_2=t*t 
            t_3=t_2*t 
            h00= 2*t_3 - 3*t_2 + 1 
            h10= t_3 - 2*t_2 + t 
            h01= (-2)*t_3 + 3*t_2 
            h11= t_3-t_2
            return h00*p0+h10*dd*m0+h01*p1+h11*dd*m1 
        return 0   

    def fnGetSplinePars(self, dz, dp, dh, m0, m1, p0, p1): 
        interval=-1 
        for i in range(self.numberofcontrolpoints-1):
            if ((dp[i] <= dz) and (dp[i+1] > dz)):
                # printf("Found interval %i \n", i)
                interval = i
        if (dz==dp[-1]):
            interval=self.numberofcontrolpoints-2 
        
        if (interval>=0):                                           #tangent calculation
            if (self.monotonic==1):                                     #Monotonic cubic spline, see Wikipedia
                
                if (dh[interval]==dh[interval+1]):
                    m0=0 
                    m1=0 
                else:
                    if (interval==0):
                        k0=(dh[interval+1]-dh[interval])/(dp[interval+1]-dp[interval]) 
                        k1=(dh[interval+2]-dh[interval+1])/(dp[interval+2]-dp[interval+1]) 
                        k2=(dh[interval+3]-dh[interval+2])/(dp[interval+3]-dp[interval+2]) 
                        km1=k0 

                    elif (interval==(self.numberofcontrolpoints-2)):
                        km1=(dh[interval]-dh[interval-1])/(dp[interval]-dp[interval-1]) 
                        k0=(dh[interval+1]-dh[interval])/(dp[interval+1]-dp[interval]) 
                        k1=k0 
                        k2=k0 

                    elif (interval==(self.numberofcontrolpoints-3)):
                        km1=(dh[interval]-dh[interval-1])/(dp[interval]-dp[interval-1]) 
                        k0=(dh[interval+1]-dh[interval])/(dp[interval+1]-dp[interval]) 
                        k1=(dh[interval+2]-dh[interval+1])/(dp[interval+2]-dp[interval+1]) 
                        k2=k0 

                    else:
                        km1=(dh[interval]-dh[interval-1])/(dp[interval]-dp[interval-1]) 
                        k0=(dh[interval+1]-dh[interval])/(dp[interval+1]-dp[interval]) 
                        k1=(dh[interval+2]-dh[interval+1])/(dp[interval+2]-dp[interval+1]) 
                        k2=(dh[interval+3]-dh[interval+2])/(dp[interval+3]-dp[interval+2]) 

                    m0=(k0+km1)/2 
                    m1=(k1+k0)/2 
                    m2=(k2+k1)/2 
                    
                    if (k0==0):
                        m0=0 
                        m1=0 

                    else :
                        alpha0=m0/k0  
                        beta0=m1/k0                     
                        if ((alpha0<0) or (beta0<0)):
                            m0=0 
    
                        elif ((alpha0*alpha0+beta0*beta0)>9):
                            tau=3/math.sqrt(alpha0*alpha0+beta0*beta0) 
                            m0=tau*alpha0*k0 
                            m1=tau*beta0*k0 

                    if (k1==0):
                        m1=0 
                        m2=0 

                    else :                        
                        alpha1=m1/k1  
                        beta1=m2/k1 
                        if ((alpha1<0) or (beta1<0)):
                            m1=0 
    
                        elif ((alpha1*alpha1+beta1*beta1)>9):
                            tau=3/math.sqrt(alpha1*alpha1+beta1*beta1) 
                            m1=tau*alpha1*k1 
                            m2=tau*beta1*k1 

            else:                                                    #Catmull-Rom spline, see Wikipedia
                if (interval==0):
                    m0=0 
                    m1=(dh[2]-dh[0])/(dp[2]-dp[0]) 

                elif (interval==(self.numberofcontrolpoints-2)):
                    m0=(dh[interval+1]-dh[interval-1])/(dp[interval+1]-dp[interval-1]) 
                    m1=0 

                else:
                    m0=(dh[interval+1]-dh[interval-1])/(dp[interval+1]-dp[interval-1]) 
                    m1=(dh[interval+2]-dh[interval])/(dp[interval+2]-dp[interval]) 
  
            p0=dh[interval] 
            p1=dh[interval+1] 
    
        return interval, m0, m1, p0, p1
        
    def fnGetSplineIntegral(self, dz1, dz2, dp, dh, damping): 
        if (dz1 > dz2):
            dz1, dz2 = dz2, dz1
        
        #check for boundaries
        dz1 = max(dz1, dp[0])
        dz2 = min(dz2, dp[-1])
        integral=0  
        d=dz1 
        while (d <= dz2):
            integral += self.fnGetSplineArea(d, dp, dh, damping)*0.5 
            d += 0.5  
        return integral 

    def fnGetSplineProductIntegral(self, dz1, dz2, dp, dh1, dh2, damping1, damping2): 
        if (dz1>dz2):
            dz1, dz2 = dz2, dz1
        
        #check for boundaries
        dz1 = max(dp[0], dz1)
        dz2 = min(dz2, dp[-1])
        
        integral=0  
        d=dz1 
        while (d<=dz2):
            integral+=self.fnGetSplineArea(d, dp, dh1, damping1)*self.fnGetSplineArea(d,dp,dh2, damping2)*0.5 
            d+=0.5 
        
        return integral