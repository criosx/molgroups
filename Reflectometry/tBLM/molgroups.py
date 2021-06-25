import numpy
import math

class nSLDObj():

    def __init__(self):
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

    def fnGetConvolutedArea(self, dz):
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
        
    def fnWriteData2File(self, f, cName, dimension, stepsize):
        f.write("z"+cName+" a"+cName+" nsl"+cName+" \n")
        dLowerLimit = self.fnGetLowerLimit()
        dUpperLimit = self.fnGetUpperLimit()
        d = numpy.floor(dLowerLimit / stepsize + 0.5) * stepsize
        
        for i in range(dimension):
            d = float(i) * stepsize
            dmirror = d - float(2*i) * stepsize
            if self.bWrapping and dmirror >= dLowerLimit:
                dAreaInc = self.fnGetConvolutedArea(d) + self.fnGetConvolutedArea(dmirror)
                dnSLDInc = (self.fnGetnSLD(d) * self.fnGetConvolutedArea(d) + self.fnGetnSLD(dmirror) * self.fnGetConvolutedArea(dmirror))
                dnSLDInc /= (self.fnGetConvolutedArea(d) + self.fnGetConvolutedArea(dmirror))
                # printf("Bin %i Area %f nSLD %e nSL %e \n", i, dAreaInc, fnGetnSLD(d), fnGetnSLD(d)*dAreaInc*stepsize)
            else:
                dAreaInc = self.fnGetConvolutedArea(d)
                dnSLDInc = self.fnGetnSLD(d)
                # printf("Bin %i z %g Area %f nSLD %e nSL %e \n", i, d, dAreaInc, fnGetnSLD(d), fnGetnSLD(d)*dAreaInc*stepsize)
                
            f.write(str(d)+" "+str(dAreaInc)+" "+str(dAreaInc*stepsize)+"\n")
        f.write("\n")
        
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
    
    def fnWriteProfile(self, aArea, anSL, dimension, stepsize, dMaxArea):
        dLowerLimit = self.fnGetLowerLimit()
        dUpperLimit = self.fnGetUpperLimit()
        if dUpperLimit==0:
            dUpperLimit = float(dimension) * stepsize
        d = numpy.floor(dLowerLimit / stepsize + 0.5) * stepsize
        
        while d<=dUpperLimit:
            i = int(d/stepsize)
            dprefactor=1
            if (i<0) and self.bWrapping:
                i = -1 * i
            if (i==0) and self.bWrapping:
                dprefactor = 2                                            #avoid too low filling when mirroring
            if (i>=0) and (i<dimension):
                dAreaInc = self.fnGetConvolutedArea(d)
                aArea[i] = aArea[i] + dAreaInc * dprefactor
                if (aArea[i]>dMaxArea):
                    dMaxArea = aArea[i]
                anSL[i] = anSL[i] + self.fnGetnSLD(d) * dAreaInc*stepsize * dprefactor
                # printf("Bin %i AreaInc %g total %g MaxArea %g nSL %f total %f \n", i, dAreaInc, aArea[i], dMaxArea, fnGetnSLD(d)*dAreaInc*stepsize, anSL[i])
            d = d + stepsize
        return dMaxArea, aArea, anSL

    def fnWriteProfileAbsorb(self, aArea, anSL, aAbsorb, dimension, stepsize, dMaxArea):
        dLowerLimit = self.fnGetLowerLimit()
        dUpperLimit = self.fnGetUpperLimit()
        if dUpperLimit == 0:
            dUpperLimit = float(dimension) * stepsize
        d = numpy.floor(dLowerLimit / stepsize + 0.5) * stepsize
        
        while d <= dUpperLimit:
            i = int(d / stepsize)
            dprefactor = 1
            # printf("Here we are %i, dimension %i \n", i, dimension)
            if (i<0) and self.bWrapping:
                i = -1 * i
            if (i == 0) and self.bWrapping:
                dprefactor = 2                                            #avoid too low filling when mirroring
            if 0 <= i < dimension:
                # printf("Bin %i Areainc %f area now %f nSLD %g Absorbinc %g Absorb now %g nSLinc %g nSL now %g \n", i, dAreaInc, aArea[i], fnGetnSLD(d), aAbsorb[i], fnGetAbsorb(d)*dAreaInc*stepsize, fnGetnSLD(d)*dAreaInc*stepsize, anSL[i])
                dAreaInc = self.fnGetConvolutedArea(d)
                aArea[i] = self.aArea[i] + dAreaInc * dprefactor
                if (aArea[i]>dMaxArea):
                    dMaxArea = aArea[i]
                anSL[i] = anSL[i] + self.fnGetnSLD(d) * dAreaInc * stepsize * dprefactor
                aAbsorb[i] = aAbsorb[i] + self.fnGetAbsorb(d) * dAreaInc * stepsize * dprefactor
                # printf("Bin %i Area %f total %f nSL %f total %f \n", i, dAreaInc, aArea[i], fnGetnSLD(d)*dAreaInc*stepsize, anSL[i])
            d = d + stepsize
        return dMaxArea
        
    def fnOverlayProfile(self, aArea, anSL, dimension, stepsize, dMaxArea):
        dLowerLimit = self.fnGetLowerLimit()
        dUpperLimit = self.fnGetUpperLimit()
        if dUpperLimit == 0:
            dUpperLimit = float(dimension) * stepsize
        d = numpy.floor(dLowerLimit / stepsize + 0.5) * stepsize
        
        while d<=dUpperLimit:
            i = int(d / stepsize)
            dprefactor = 1
            if (i < 0) and self.bWrapping:
                i=-1*i
            if (i == 0) and self.bWrapping:
                dprefactor = 2                                           #avoid too low filling when mirroring
            if (i >= 0) and i<dimension:
                dAreaInc = self.fnGetConvolutedArea(d)
                temparea = dAreaInc * dprefactor + aArea[i]
                if temparea<=dMaxArea:
                    aArea[i] = aArea[i] + dAreaInc * dprefactor
                    anSL[i] = anSL[i] + self.fnGetnSLD(d) * dAreaInc * stepsize * dprefactor
                else:
                    if (temparea-dMaxArea) <= aArea[i]:                           #overfill is not larger than existing area
                        anSL[i] = anSL[i] * (1-((temparea-dMaxArea)/aArea[i]))    #eliminate the overfilled portion using original content
                        anSL[i] = anSL[i] + self.fnGetnSLD(d) * dAreaInc * stepsize * dprefactor
                        aArea[i] = dMaxArea
                        # printf("Replace: Bin %i temparea %g Areainc %g area now %g dMaxArea %g nSLD %g nSLinc %g nSL now %g \n", i, temparea, dAreaInc, aArea[i], dMaxArea, fnGetnSLD(d), fnGetnSLD(d)*dAreaInc*stepsize, anSL[i])
                    else:                                                         #overfill is larger!!, this is non-physical
                        anSL[i] = self.fnGetnSLD(d) * dMaxArea*stepsize
                        aArea[i] = dMaxArea
            d = d + stepsize

    def fnOverlayProfileAbsorb(self, aArea, anSL, aAbsorb, dimension, stepsize, dMaxArea):
        dLowerLimit = self.fnGetLowerLimit()
        dUpperLimit = self.fnGetUpperLimit()
        if dUpperLimit==0:
            dUpperLimit = float(dimension)*stepsize
        d = math.floor(dLowerLimit / stepsize + 0.5) * stepsize
        
        while d<=dUpperLimit:
            i = int(d/stepsize)
            dprefactor = 1
            # printf("Here we are %i, dimension %i, maxarea %f \n", i, dimension, dMaxArea)
            if i < 0 and self.bWrapping:
                i = -1 * i
            if i == 0 and self.bWrapping:
                dprefactor = 2                                                   #avoid too low filling when mirroring
            if 0 <= i < dimension:
                dAreaInc = self.fnGetConvolutedArea(d)
                temparea = dAreaInc * dprefactor+aArea[i]
                if temparea>dMaxArea:
                    # printf("Bin %i Areainc %f area now %f nSLD %g Absorbinc %g Absorb now %g nSLinc %g nSL now %g \n", i, dAreaInc, aArea[i], fnGetnSLD(d), aAbsorb[i], fnGetAbsorb(d)*dAreaInc*stepsize, fnGetnSLD(d)*dAreaInc*stepsize, anSL[i])
                    # eliminate the overfilled portion using original content
                    anSL[i] = anSL[i] * (1-((temparea-dMaxArea)/aArea[i]))
                    anSL[i] = anSL[i] + self.fnGetnSLD(d) * dAreaInc * stepsize * dprefactor
                    # eliminate the overfilled portion using original content
                    aAbsorb[i] = aAbsorb[i] * (1-((temparea-dMaxArea)/aArea[i]))
                    aAbsorb[i] = aAbsorb[i] + self.fnGetAbsorb(d) * dAreaInc * stepsize*dprefactor
                    aArea[i] = dMaxArea
                else:
                    # printf("Bin %i Areainc %f area now %f nSLD %g Absorbinc %g Absorb now %g nSLinc %g nSL now %g \n", i, dAreaInc, aArea[i], fnGetnSLD(d), aAbsorb[i], fnGetAbsorb(d)*dAreaInc*stepsize, fnGetnSLD(d)*dAreaInc*stepsize, anSL[i])
                    aArea[i] = aArea[i] + dAreaInc * dprefactor
                    anSL[i] = anSL[i] + self.fnGetnSLD(d) * dAreaInc * stepsize * dprefactor
                    aAbsorb[i] = aAbsorb[i] + self.fnGetAbsorb(d) * dAreaInc * stepsize * dprefactor
            d = d+stepsize


# ------------------------------------------------------------------------------------------------------
# Function Object Implementation
# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Simple Objects
# ------------------------------------------------------------------------------------------------------
class Box2Err(nSLDObj):

    def __init__(self, dz=20, dsigma1=2, dsigma2=2, dlength=10, dvolume=10, dnSL=0, dnumberfraction=1):
        super().__init__()
        self.z = dz
        self.sigma1 = dsigma1
        self.sigma2 = dsigma2
        self.l = dlength
        self.vol = dvolume
        self.nSL = dnSL
        self.nf = dnumberfraction
        self.nsldbulk_store = 0.
        self.nSL2 = 0.

    # Gaussian function definition, integral is volume, return value is area at position z
    def fnGetArea(self, dz):
        if (self.l != 0) and (self.sigma1 != 0) and (self.sigma2 != 0):
            result = math.erf((dz - self.z + 0.5 * self.l) / math.sqrt(2) / self.sigma1)
            result -= math.erf((dz - self.z - 0.5 * self.l) / math.sqrt(2) / self.sigma2)
            result *= (self.vol / self.l) * 0.5
            result *= self.nf
            return  result
        else:
            return 0.

    def fnGetnSL(self, bulknsld):
        if self.bProtonExchange:
            if self.vol != 0:
                return ((bulknsld + 0.56e-6) * self.nSL2 + (6.36e-6 - bulknsld) * self.nSL) / (6.36e-6 + 0.56e-6)
            else:
                return 0.
        else:
            return self.nSL

    # constant nSLD
    def fnGetnSLD(self, dz, bulknsld = 0.):
        if bulknsld != 0:
            self.nsldbulk_store = bulknsld
        if self.vol != 0:
            if self.bProtonExchange:
                if bulknsld == 0.:
                    bulknsld = self.nsldbulk_store
                return ((bulknsld+0.56e-6)*self.nSL2+(6.36e-6-bulknsld)*self.nSL)/(6.36e-6+0.56e-6)/self.vol
            else:
                return self.nSL/self.vol
        else:
            return 0

    # Gaussians are cut off below and above 3 sigma double
    def fnGetLowerLimit(self):
        return self.z - 0.5 * self.l - 3 * self.sigma1

    def fnGetUpperLimit(self):
        return self.z + 0.5 * self.l + 3 * self.sigma2

    def fnSetnSL(self, _nSL, _nSL2):
        self.nSL = _nSL
        self.nSL2 = _nSL2
        self.bProtonExchange = True

    def fnSetSigma(self, sigma1, sigma2=0.):
        self.sigma1 = sigma1
        if sigma2 == 0:
            sigma2 = sigma1
        self.sigma2 = sigma2

    def fnSetZ(self, dz):
        self.z = dz

    def fnWritePar2File(self, fp, cName, dimension, stepsize):
        fp.write("Box2Err "+cName+" z "+str(self.z)+" sigma1 "+str(self.sigma1)+" sigma2 "+str(self.sigma2)+" l "
                 + str(self.l)+" vol "+str(self.vol)+" nSL "+str(self.nSL)+" nSL2 "+str(self.nSL2)+" nf "
                 + str(self.nf)+" \n")
        nSLDObj.fnWriteData2File(self, fp, cName, dimension, stepsize)

# ------------------------------------------------------------------------------------------------------
# Headgroups
# ------------------------------------------------------------------------------------------------------
class PC(nSLDObj):
    def __init__(self):
        super().__init__()
        self.z
        self.cg = Box2Err()
        self.phosphate = Box2Err()
        self.choline = Box2Err()
        self.cg.l = 4.21 
        self.phosphate.l = 3.86
        self.choline.l = 6.34
        self.cg.sigma1, self.cg.sigma2 = 2.53, 2.29
        self.phosphate.sigma1, self.phosphate.sigma2 = 2.29, 2.02
        self.choline.sigma1, self.choline.sigma2 = 2.02, 2.26
        self.l = 9.575
        self.cg.vol=147 
        self.phosphate.vol=54 
        self.choline.vol=120
        self.cg.nSL=3.7755e-4 
        self.phosphate.nSL=2.8350e-4 
        self.choline.nSL=-6.0930e-5
        self.cg.nf=1 
        self.phosphate.nf=1 
        self.choline.nf=1
        self.vol=self.cg.vol+self.phosphate.vol+self.choline.vol
        self.nSL=self.cg.nSL+self.phosphate.nSL+self.choline.nSL
        self.fnAdjustParameters()

    def fnAdjustParameters(self):
        self.cg.z = self.z - 0.5*self.l + 0.5*self.cg.l
        self.choline.z = self.z + 0.5*self.l - 0.5*self.choline.l
        self.phosphate.z = (self.cg.z+0.5*self.cg.l+self.choline.z-self.choline.l*0.5)/2
    
    def fnGetLowerLimit(self):
        return self.cg.fnGetLowerLimit()
    
    def fnGetUpperLimit(self):
        return self.choline.fnGetUpperLimit()
    
    def fnGetTotalnSL(self):
        return self.cg.nSL + self.phosphate.nSL + self.choline.nSL
        
    def fnGetnSL(self, bulknsld):
        if self.bProtonExchange:
            if self.vol != 0:
                return ((bulknsld + 0.56e-6) * self.nSL2 + (6.36e-6 - bulknsld) * self.nSL) / (6.36e-6 + 0.56e-6)
            else:
                return 0.
        else:
            return self.nSL

    def fnGetArea(self, dz):
        return (self.cg.fnGetArea(dz)+self.phosphate.fnGetArea(dz)+self.choline.fnGetArea(dz))*self.nf
    
    def fnGetnSLD(self, dz, bulknsld = 0):
        cgarea=self.cg.fnGetArea(dz)
        pharea=self.phosphate.fnGetArea(dz)
        charea=self.choline.fnGetArea(dz)
        sum=cgarea+pharea+charea
        #bulknsld=-0.56e-6  
        if (sum == 0):
            return 0
        else:
            return (self.cg.fnGetnSLD(dz)*cgarea+\
                self.phosphate.fnGetnSLD(dz)*\
                pharea+self.choline.fnGetnSLD(dz)*charea)/sum
    
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
    
    def fnWritePar2File (self, fp, cName, dimension, stepsize):
        pass

class PCm(PC):
    def __init__ (self):
        super().__init__()
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
        self.phosphate.z = (self.cg.z-0.5*self.cg.l+self.choline.z+self.choline.l*0.5)/2
    def fnGetLowerLimit(self):
        return self.choline.fnGetLowerLimit()
    def fnGetUpperLimit(self): 
        return self.cg.fnGetUpperLimit()
    def fnWritePar2File (self,fp, cName, dimension, stepsize):
        pass

    def fnGetnSLD(self, dz, bulknsld = 0):
        cgarea=self.cg.fnGetArea(dz)
        pharea=self.phosphate.fnGetArea(dz)
        charea=self.choline.fnGetArea(dz)
        sum=cgarea+pharea+charea
        #bulknsld=-0.56e-6  
        if (sum == 0):
            return 0
        else:
            return (self.cg.fnGetnSLD(dz)*cgarea+\
                self.phosphate.fnGetnSLD(dz)*\
                pharea+self.choline.fnGetnSLD(dz)*charea)/sum
        
        """
        
        if (sum == 0):
            return 0
        elif sum!=0:
            if self.vol != 0:
                if self.bProtonExchange:
                    if bulknsld == 0.:
                        bulknsld = self.nsldbulk_store
                     return ((bulknsld+0.56e-6)*self.nSL2+(6.36e-6-bulknsld)*self.nSL)/(6.36e-6+0.56e-6)/self.vol
                    return (self.cg.fnGetnSLD(dz)*cgarea+\
                    self.phosphate.fnGetnSLD(dz)*\
                    pharea+self.choline.fnGetnSLD(dz)*charea)/sum
        elif bulknsld != 0:
            self.nsldbulk_store = bulknsld
        if sum !=0   
        
       """

        
            
# ------------------------------------------------------------------------------------------------------
# Lipid Bilayer
# ------------------------------------------------------------------------------------------------------
class BLM_quaternary(nSLDObj):
    def __init__(self):
        
        super().__init__()
        self.headgroup1 = PCm()
        self.lipid1 = Box2Err()
        self.methyl1 = Box2Err()
        self.methyl2 = Box2Err()
        self.lipid2 = Box2Err()
        self.headgroup2 = PC()                          # PC head group
        self.headgroup1_2 = Box2Err()                        # second headgroups
        self.headgroup2_2 = Box2Err()
        self.headgroup1_3 = Box2Err()
        self.headgroup2_3 = Box2Err()
        
        self.defect_hydrocarbon = Box2Err()
        self.defect_headgroup = Box2Err()
        
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

        self.vf_bilayer = 1.0
        self.absorb = 0.
        self.l_lipid1 = 11.
        self.l_lipid2 = 11.
        self.bulknsld = -0.56e-6
        self.normarea = 60.
        self.startz = 50.
        self.sigma = 2.
        self.radius_defect = 100.
        self.hc_substitution_1 = 0
        self.hc_substitution_2 = 0

        self.fnAdjustParameters()
        
    def fnInit(self, va1, na1, vm1, nm1, vh1, nh1, lh1, va2, na2, vm2, nm2, vh2, nh2, lh2, va3, na3, vm3, nm3, vh3, nh3, lh3, vc, nc):
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
        # printf("Enter AdjustParameters \n");
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
        
        # printf("ssBLM: normarea %lf \n",normarea);
        self.lipid2.l = l_ohc
        self.lipid2.vol = V_ohc
        self.lipid2.nSL = nSL_ohc
        self.lipid2.nf = c_s_ohc * c_A_ohc * c_V_ohc
        # printf("c: c_s_ohc %lf c_A_ohc %lf c_V_ohc %lf \n", c_s_ohc, c_A_ohc, c_V_ohc);
        
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
        # printf("c: c_s_ihc %lf c_A_ihc %lf nf_ihc_lipid %lf hc_substitution_1 %lf \n", c_s_ihc, c_A_ihc, nf_ihc_lipid, hc_substitution_1);
        
        self.lipid1.z= self.startz + self.headgroup1.l + 0.5 * self.lipid1.l
        self.headgroup1.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1.l)
        self.headgroup1_2.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1_2.l)
        self.headgroup1_3.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1_3.l)
        self.methyl1.z = self.lipid1.z + 0.5 * (self.lipid1.l + self.methyl1.l)
        self.methyl2.z = self.methyl1.z + 0.5 * (self.methyl1.l + self.methyl2.l)
        self.lipid2.z = self.methyl2.z + 0.5 * (self.methyl2.l + self.lipid2.l)
        self.headgroup2.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2.l)
        self.headgroup2_2.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2_2.l)
        self.headgroup2_3.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2_3.l)
        # printf("nf bme %lf tether %lf tetherg %lf lipid1 %lf headgroup1 %lf headgroup1_2 %lf headgroup1_3 %lf methyl1 %lf methyl2 %lf lipid2 %lf headgroup2 %lf headgroup2_2 %lf headgroup2_3 %lf \n", bME.nf, tether.nf, tetherg.nf, lipid1.nf, headgroup1.nf, headgroup1_2.nf, headgroup1_3.nf, methyl1.nf, methyl2.nf, lipid2.nf, headgroup2.nf, headgroup2_2.nf, headgroup2_3.nf);
        
        # defects
        hclength = self.lipid1.l + self.methyl1.l + self.methyl2.l + self.lipid2.l
        hglength = self.headgroup1.l + self.headgroup2.l
        
        if self.radius_defect<(0.5*(hclength+hglength)):
            self.radius_defect = 0.5 * (hclength+hglength)
        # printf("defect_radius %lf hclength %lf \n",radius_defect, hclength);
        
        volhalftorus = 3.14159265359 * 3.14159265359 * (self.radius_defect - (2 * hclength / 3 / 3.14159265359)) * hclength * hclength / 4
        volcylinder = 3.14159265359 * self.radius_defect * self.radius_defect * hclength
        # printf("volhalftorus %lf volcylinder %lf \n", volhalftorus, volcylinder);
        defectarea = volhalftorus / volcylinder * (1 - self.vf_bilayer) * self.normarea
        # printf("defectarea %lf \n", defectarea);
        
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
        self.defect_headgroup.nSL = defectratio * (self.headgroup2.fnGetnSL(self.bulknsld) * self.headgroup2.nf + self.headgroup2_2.fnGetnSL(self.bulknsld) * self.headgroup2_2.nf + self.headgroup2_3.fnGetnSL(self.bulknsld) * self.headgroup2_3.nf)
        self.defect_headgroup.fnSetSigma(self.sigma)
        self.defect_headgroup.nf = 1

    # Return value is area at position z
    def fnGetArea(self, dz):
        result  = self.lipid1.fnGetArea(dz) + self.headgroup1.fnGetArea(dz) + self.methyl1.fnGetArea(dz)
        result += self.methyl2.fnGetArea(dz) + self.lipid2.fnGetArea(dz) + self.headgroup2.fnGetArea(dz)
        result += self.headgroup1_2.fnGetArea(dz) + self.headgroup2_2.fnGetArea(dz) + self.headgroup1_3.fnGetArea(dz)
        result += self.headgroup2_3.fnGetArea(dz) + self.defect_hydrocarbon.fnGetArea(dz) + self.defect_headgroup.fnGetArea(dz)
        return result

    # return value of center of the membrane
    def fnGetCenter(self):
        return self.methyl1.z + 0.5 * self.methyl1.l

    # get nSLD from molecular subgroups
    def fnGetnSLD(self, dz):
        # printf("Enter fnGetnSLD \n");
        lipid1area = self.lipid1.fnGetArea(dz)
        headgroup1area = self.headgroup1.fnGetArea(dz)
        headgroup1_2_area = self.headgroup1_2.fnGetArea(dz)
        headgroup1_3_area = self.headgroup1_3.fnGetArea(dz)
        methyl1area = self.methyl1.fnGetArea(dz)
        methyl2area = self.methyl2.fnGetArea(dz)
        lipid2area = self.lipid2.fnGetArea(dz)
        headgroup2area = self.headgroup2.fnGetArea(dz)
        headgroup2_2_area = self.headgroup2_2.fnGetArea(dz)
        headgroup2_3_area = self.headgroup2_3.fnGetArea(dz)
        defect_hydrocarbon_area = self.defect_hydrocarbon.fnGetArea(dz)
        defect_headgroup_area = self.defect_headgroup.fnGetArea(dz)
        
        sum  = lipid1area + headgroup1area + methyl1area + methyl2area + lipid2area + headgroup2area + headgroup1_2_area
        sum += headgroup2_2_area + headgroup1_3_area + headgroup2_3_area + defect_headgroup_area + defect_hydrocarbon_area
        
        # printf("%e \n", defect_headgroup.fnGetnSLD(dz));
        if sum==0:
            return 0
        else:
            result  = self.headgroup1.fnGetnSLD(dz, self.bulknsld) * headgroup1area
            result += self.headgroup1_2.fnGetnSLD(dz, self.bulknsld) * headgroup1_2_area
            result += self.headgroup1_3.fnGetnSLD(dz, self.bulknsld) * headgroup1_3_area
            result += self.lipid1.fnGetnSLD(dz) * lipid1area
            result += self.methyl1.fnGetnSLD(dz) * methyl1area
            result += self.methyl2.fnGetnSLD(dz) * methyl2area
            result += self.lipid2.fnGetnSLD(dz) * lipid2area
            result += self.headgroup2.fnGetnSLD(dz, self.bulknsld) * headgroup2area
            result += self.headgroup2_2.fnGetnSLD(dz, self.bulknsld) * headgroup2_2_area
            result += self.headgroup2_3.fnGetnSLD(dz, self.bulknsld) * headgroup2_3_area
            result += self.defect_hydrocarbon.fnGetnSLD(dz) * defect_hydrocarbon_area
            result += self.defect_headgroup.fnGetnSLD(dz) * defect_headgroup_area
            result /= sum
            return result

    # Use limits of molecular subgroups
    def fnGetLowerLimit(self):
        return self.headgroup1.fnGetLowerLimit()

    def fnGetUpperLimit(self):
        a = self.headgroup2.fnGetUpperLimit()
        b = self.headgroup2_2.fnGetUpperLimit()
        c = self.headgroup2_3.fnGetUpperLimit()
        return max([a, b, c])

    def fnSet(self, _sigma, _bulknsld, _startz, _l_lipid1, _l_lipid2, _vf_bilayer, _nf_lipid_2=0., _nf_lipid_3=0., _nf_chol=0., _hc_substitution_1=0., _hc_substitution_2=0., _radius_defect=100.):
        self.sigma = _sigma
        self.bulknsld = _bulknsld
        self.startz = _startz
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
        self.headgroup1.fnSetSigma(sigma)
        self.headgroup1_2.fnSetSigma(sigma)
        self.headgroup1_3.fnSetSigma(sigma)
        self.lipid1.fnSetSigma(sigma,sigma+2)
        self.methyl1.fnSetSigma(sigma+2,sigma+2)
        self.methyl2.fnSetSigma(sigma+2,sigma+2)
        self.lipid2.fnSetSigma(sigma+2,sigma)
        self.headgroup2.fnSetSigma(sigma)
        self.headgroup2_2.fnSetSigma(sigma)
        self.headgroup2_3.fnSetSigma(sigma)
        self.defect_hydrocarbon.fnSetSigma(sigma)
        self.defect_headgroup.fnSetSigma(sigma)

    def fnWriteProfile(self, aArea, anSL, dimension, stepsize, dMaxArea):
        _, aArea, anSL = nSLDObj.fnWriteProfile(self, aArea, anSL, dimension, stepsize, dMaxArea)
        return self.normarea, aArea, anSL

    def fnWritePar2File(self, fp, cName, dimension, stepsize):
        self.headgroup1.fnWritePar2File(fp, "blm_headgroup1", dimension, stepsize)
        self.headgroup1_2.fnWritePar2File(fp, "blm_headgroup1_2", dimension, stepsize)
        self.headgroup1_3.fnWritePar2File(fp, "blm_headgroup1_3", dimension, stepsize)
        self.lipid1.fnWritePar2File(fp, "blm_lipid1", dimension, stepsize)
        self.methyl1.fnWritePar2File(fp, "blm_methyl1", dimension, stepsize)
        self.methyl2.fnWritePar2File(fp, "blm_methyl2", dimension, stepsize)
        self.lipid2.fnWritePar2File(fp, "blm_lipid2", dimension, stepsize)
        self.headgroup2.fnWritePar2File(fp, "blm_headgroup2", dimension, stepsize)
        self.headgroup2_2.fnWritePar2File(fp, "blm_headgroup2_2", dimension, stepsize)
        self.headgroup2_3.fnWritePar2File(fp, "blm_headgroup2_3", dimension, stepsize)
        self.defect_hydrocarbon.fnWritePar2File(fp, "blm_defect_hc", dimension, stepsize)
        self.defect_headgroup.fnWritePar2File(fp, "blm_defect_hg", dimension, stepsize)
        self.fnWriteConstant(fp, "blm_normarea", self.normarea, 0, dimension, stepsize)

##---start of modfied code-----


class ssBLM_quaternary(nSLDObj):
    def __init__(self):

        super().__init__()
        self.substrate  = Box2Err()
        self.siox       = Box2Err()
        self.headgroup1 = PCm()
        self.lipid1 = Box2Err()
        self.methyl1 = Box2Err()
        self.methyl2 = Box2Err()
        self.lipid2 = Box2Err()
        self.headgroup2 = PC()                          # PC head group
        self.headgroup1_2 = Box2Err()                        # second headgroups
        self.headgroup2_2 = Box2Err()
        self.headgroup1_3 = Box2Err()
        self.headgroup2_3 = Box2Err()
        
        self.defect_hydrocarbon = Box2Err()
        self.defect_headgroup    = Box2Err()

        self.substrate.l = 20
        self.substrate.z = 10
        self.substrate.nf = 1
        self.substrate.sigma1 = 2.0
        self.rho_substrate = 1
        self.l_siox =1
        self.rho_siox = 1

        self.siox.l=20
        self.siox.z=30
        self.siox.nf=1
        self.siox.fnSetSigma(2.0)
	
        self.volacyllipid=925
        self.nslacyllipid=-2.67e-4
        self.volmethyllipid=98.8
        self.nslmethyllipid=-9.15e-5



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
        

        self.hc_substitution_1=0
        self.hc_substitution_2=0
       
       
        self.headgroup1_2.vol=330       #was 330
        self.headgroup2_2.vol=330      #was 330
        self.headgroup1_3.vol=330       #was 330
        self.headgroup2_3.vol=330       #was 330
        self.headgroup1_2.nSL=6.0012e-4 #was 6.0122e-4
        self.headgroup2_2.nSL=6.0012e-4 #was 6.0122e-4
        self.headgroup1_3.nSL=6.0012e-4 #was 6.0122e-4
        self.headgroup2_3.nSL=6.0012e-4 #was 6.0122e-4
        self.headgroup1_2.l=9.5
        self.headgroup2_2.l=9.5
        self.headgroup1_3.l=9.5
        self.headgroup2_3.l=9.5
	
        self.volacyllipid_2=925
        self.nslacyllipid_2=-2.67e-4
        self.volmethyllipid_2=98.8
        self.nslmethyllipid_2=-9.15e-5


        self.volacyllipid_3=925
        self.nslacyllipid_3=-2.67e-4
        self.volmethyllipid_3=98.8
        self.nslmethyllipid_3=-9.15e-5

        self.volchol=630
        self.nslchol=1.3215e-4    
        self.nf_lipid_2=0                                                                #for preparing towards a general bilayer class
        self.nf_lipid_3=0
        self.nf_chol=0
        self.bulknsld=-0.56e-6     

        self.nf_lipid_2 = 0.
        self.nf_lipid_3 = 0.
        self.nf_chol = 0.

        self.vf_bilayer = 1.0
        self.absorb = 0.
        self.l_lipid1 = 11.
        self.l_lipid2 = 11.
        self.bulknsld = -0.56e-6
        self.normarea = 60.
        self.startz = 50.
        self.sigma = 2.
        self.radius_defect = 100.
        self.hc_substitution_1 = 0
        self.hc_substitution_2 = 0 
        self.global_rough = 2.0
        
        
        self.fnAdjustParameters()


# is this necessary?
    def fnInit(self, va1, na1, vm1, nm1, vh1, nh1, lh1, va2, na2, vm2, nm2, vh2, nh2, lh2, va3, na3, vm3, nm3, vh3, nh3, lh3, vc, nc):
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

        #def ssBLM_quaternary::fnAdjustParameters()
        #def fnAdjustParameters(self):  
        #printf("Enter AdjustParameters \n")

        self.substrate.fnSetSigma(self.global_rough)
        self.siox.fnSetSigma(self.global_rough)
        self.fnSetSigma(self.sigma)

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
    
#substrate 
        self.substrate.vol=self.normarea*self.substrate.l
        self.substrate.nSL=self.rho_substrate*self.substrate.vol
        self.siox.l= self.l_siox
        self.siox.vol=self.normarea*self.siox.l
        self.siox.nSL=self.rho_siox*self.siox.vol
    
    #set all lengths
        self.siox.z=self.substrate.l+0.5*self.siox.l
        self.lipid1.z= self.startz + self.headgroup1.l + 0.5 * self.lipid1.l
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
    
        if self.radius_defect<(0.5*(hclength+hglength)):
            self.radius_defect = 0.5 * (hclength+hglength)
        # printf("defect_radius %lf hclength %lf \n",radius_defect, hclength)
    
        volhalftorus=3.14159265359*3.14159265359*(self.radius_defect-(2*hclength/3/3.14159265359))*hclength*hclength/4
        volcylinder=3.14159265359*self.radius_defect*self.radius_defect*hclength
    #printf("volhalftorus %lf volcylinder %lf \n", volhalftorus, volcylinder)
        defectarea=volhalftorus/volcylinder*(1- self.vf_bilayer)* self.normarea
    #printf("defectarea %lf \n", defectarea)
    

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
        self.defect_headgroup.nSL = defectratio * (self.headgroup2.fnGetnSL(self.bulknsld) * self.headgroup2.nf + self.headgroup2_2.fnGetnSL(self.bulknsld) * self.headgroup2_2.nf + self.headgroup2_3.fnGetnSL(self.bulknsld) * self.headgroup2_3.nf)
        self.defect_headgroup.fnSetSigma(self.sigma)
        self.defect_headgroup.nf = 1
    #printf("Exit AdjustParameters \n")
    

#Return value is area at position z
    def fnGetArea(self, dz):
        result  = self.substrate.fnGetArea(dz)+ self.siox.fnGetArea(dz) + self.lipid1.fnGetArea(dz) + self.headgroup1.fnGetArea(dz) + self.methyl1.fnGetArea(dz)
        result += self.methyl2.fnGetArea(dz) + self.lipid2.fnGetArea(dz) + self.headgroup2.fnGetArea(dz)
        result += self.headgroup1_2.fnGetArea(dz) + self.headgroup2_2.fnGetArea(dz) + self.headgroup1_3.fnGetArea(dz)
        result += self.headgroup2_3.fnGetArea(dz) + self.defect_hydrocarbon.fnGetArea(dz) + self.defect_headgroup.fnGetArea(dz)
        return result


#get nSLD from molecular subgroups
    def fnGetnSLD(self, dz):
        # printf("Enter fnGetnSLD \n")
        substratearea=self.substrate.fnGetArea(dz)
        sioxarea=self.siox.fnGetArea(dz)
        lipid1area = self.lipid1.fnGetArea(dz)
        headgroup1area = self.headgroup1.fnGetArea(dz)
        headgroup1_2_area = self.headgroup1_2.fnGetArea(dz)
        headgroup1_3_area = self.headgroup1_3.fnGetArea(dz)
        methyl1area = self.methyl1.fnGetArea(dz)
        methyl2area = self.methyl2.fnGetArea(dz)
        lipid2area = self.lipid2.fnGetArea(dz)
        headgroup2area = self.headgroup2.fnGetArea(dz)
        headgroup2_2_area = self.headgroup2_2.fnGetArea(dz)
        headgroup2_3_area = self.headgroup2_3.fnGetArea(dz)
        defect_hydrocarbon_area = self.defect_hydrocarbon.fnGetArea(dz)
        defect_headgroup_area = self.defect_headgroup.fnGetArea(dz)
        
        sum  = substratearea + sioxarea + lipid1area + headgroup1area + methyl1area + methyl2area + lipid2area + headgroup2area + headgroup1_2_area
        sum += headgroup2_2_area + headgroup1_3_area + headgroup2_3_area + defect_headgroup_area + defect_hydrocarbon_area
        
        # printf("%e \n", defect_headgroup.fnGetnSLD(dz))
    

        if sum==0:
            return 0
        else:
            result = self.substrate.fnGetnSLD(dz, self.bulknsld) * substratearea
            result += self.siox.fnGetnSLD(dz, self.bulknsld) * sioxarea
            result += self.headgroup1.fnGetnSLD(dz, self.bulknsld) * headgroup1area
            result += self.headgroup1_2.fnGetnSLD(dz, self.bulknsld) * headgroup1_2_area
            result += self.headgroup1_3.fnGetnSLD(dz, self.bulknsld) * headgroup1_3_area
            result += self.lipid1.fnGetnSLD(dz) * lipid1area
            result += self.methyl1.fnGetnSLD(dz) * methyl1area
            result += self.methyl2.fnGetnSLD(dz) * methyl2area
            result += self.lipid2.fnGetnSLD(dz) * lipid2area
            result += self.headgroup2.fnGetnSLD(dz, self.bulknsld) * headgroup2area
            result += self.headgroup2_2.fnGetnSLD(dz, self.bulknsld) * headgroup2_2_area
            result += self.headgroup2_3.fnGetnSLD(dz, self.bulknsld) * headgroup2_3_area
            result += self.defect_hydrocarbon.fnGetnSLD(dz) * defect_hydrocarbon_area
            result += self.defect_headgroup.fnGetnSLD(dz) * defect_headgroup_area
            result /= sum
            return result

 # Use limits of molecular subgroups
    def fnGetLowerLimit(self):
        return self.substrate.fnGetLowerLimit()

    def fnGetUpperLimit(self):
        a = self.headgroup2.fnGetUpperLimit()
        b = self.headgroup2_2.fnGetUpperLimit()
        c = self.headgroup2_3.fnGetUpperLimit()
        return max([a, b, c])




#void ssBLM_quaternary::fnSet(double _sigma, double _global_rough, double _rho_substrate, double _bulknsld, double _rho_siox, double _l_siox, double _l_submembrane,  double _l_lipid1, double _l_lipid2, double _vf_bilayer, double _nf_lipid_2, double _nf_lipid_3, double _nf_chol, double _hc_substitution_1, double _hc_substitution_2, double _radius_defect){
    
#printf("Enter fnSet \n")
    

    def fnSet(self, _sigma, _bulknsld, _global_rough, _rho_substrate, _rho_siox, _l_siox, _l_submembrane, _l_lipid1, _l_lipid2, _vf_bilayer, _nf_lipid_2 = 0., _nf_lipid_3=0., _nf_chol=0., _hc_substitution_1=0., _hc_substitution_2=0., _radius_defect=100.):
        self.sigma = _sigma
        self.bulknsld = _bulknsld
        self.global_rough= _global_rough
        self.rho_substrate=_rho_substrate
        self.rho_siox= _rho_siox
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

    #printf("Exit fnSet \n")

    def fnSetSigma(self, sigma):
        self.substrate.sigma2=self.global_rough
        self.siox.sigma1=self.global_rough
        self.siox.sigma2=self.global_rough
        self.headgroup1.fnSetSigma(sigma)
        self.headgroup1_2.fnSetSigma(sigma)
        self.headgroup1_3.fnSetSigma(sigma)
        self.lipid1.fnSetSigma(sigma,sigma+2)
        self.methyl1.fnSetSigma(sigma+2,sigma+2)
        self.methyl2.fnSetSigma(sigma+2,sigma+2)
        self.lipid2.fnSetSigma(sigma+2,sigma)
        self.headgroup2.fnSetSigma(sigma)
        self.headgroup2_2.fnSetSigma(sigma)
        self.headgroup2_3.fnSetSigma(sigma)
        self.defect_hydrocarbon.fnSetSigma(sigma)
        self.defect_headgroup.fnSetSigma(sigma)



    def fnWriteProfile(self, aArea, anSL, dimension, stepsize, dMaxArea):
        _, aArea, anSL = nSLDObj.fnWriteProfile(self, aArea, anSL, dimension, stepsize, dMaxArea)
        return self.normarea, aArea, anSL


   #char *str = new char[80]
    
    #fprintf(fp, "PC %s z %lf l %lf nf %lf \n",cName, z, l, nf)
 
       #delete []str


    def fnWritePar2File(self, fp, cName, dimension, stepsize):
        self.substrate.fnWritePar2File(fp, "substrate", dimension, stepsize)
        self.siox.fnWritePar2File(fp, "siox", dimension, stepsize)
        self.headgroup1.fnWritePar2File(fp, "blm_headgroup1", dimension, stepsize)
        self.headgroup1_2.fnWritePar2File(fp, "blm_headgroup1_2", dimension, stepsize)
        self.headgroup1_3.fnWritePar2File(fp, "blm_headgroup1_3", dimension, stepsize)
        self.lipid1.fnWritePar2File(fp, "blm_lipid1", dimension, stepsize)
        self.methyl1.fnWritePar2File(fp, "blm_methyl1", dimension, stepsize)
        self.methyl2.fnWritePar2File(fp, "blm_methyl2", dimension, stepsize)
        self.lipid2.fnWritePar2File(fp, "blm_lipid2", dimension, stepsize)
        self.headgroup2.fnWritePar2File(fp, "blm_headgroup2", dimension, stepsize)
        self.headgroup2_2.fnWritePar2File(fp, "blm_headgroup2_2", dimension, stepsize)
        self.headgroup2_3.fnWritePar2File(fp, "blm_headgroup2_3", dimension, stepsize)
        self.defect_hydrocarbon.fnWritePar2File(fp, "blm_defect_hc", dimension, stepsize)
        self.defect_headgroup.fnWritePar2File(fp, "blm_defect_hg", dimension, stepsize)
        self.fnWriteConstant(fp, "blm_normarea", self.normarea, 0, dimension, stepsize)

#------------------------------------------------------------------------------------------------------
# Tethered Lipid bilayer - binary system
#------------------------------------------------------------------------------------------------------
class tBLM_quaternary(nSLDObj):
    def __init__(self):

        super().__init__() 
        self.substrate  = Box2Err()
        self.bME        = Box2Err()
        self.tether	   = Box2Err()
        self.tetherg    = Box2Err()
        self.headgroup1 = PCm()                                                    #mirrored PC head group
        self.lipid1     = Box2Err()
        self.methyl1    = Box2Err()
        self.methyl2	   = Box2Err()
        self.lipid2	   = Box2Err()
        self.headgroup2 = PC()                                                          #PC head group
        self.headgroup1_2 = Box2Err()                                                  #second headgroups
        self.headgroup2_2 = Box2Err()
        self.headgroup1_3 = Box2Err()
        self.headgroup2_3 = Box2Err()
        
        self.defect_hydrocarbon = Box2Err()
        self.defect_headgroup    = Box2Err()

        self.substrate.l = 20
        self.substrate.z = 10
        self.substrate.nf = 1
        self.substrate.sigma1 = 2.0
	
        self.bME.vol=110
        self.bME.nSL=3.243e-5
        self.bME.l=5.2
        self.tether.vol=380
        self.tether.nSL=2.1864e-4
        self.tetherg.vol=110
        self.tetherg.nSL=1.8654e-4

        self.volacyllipid=925
        self.nslacyllipid=-2.67e-4
        self.volmethyllipid=98.8
        self.nslmethyllipid=-9.15e-5
        self.volmethyltether=98.8
        self.nslmethyltether=-9.15e-5
        self.volacyltether=982
        self.nslacyltether=-2.85e-4
	
        self.hc_substitution_1=0
        self.hc_substitution_2=0
        self.headgroup1_2.vol=330       #was 330
        self.headgroup2_2.vol=330       #was 330
        self.headgroup1_3.vol=330       #was 330
        self.headgroup2_3.vol=330       #was 330
        self.headgroup1_2.nSL=6.0012e-4 # was 6.0122e-4
        self.headgroup2_2.nSL=6.0012e-4 # was 6.0122e-4
        self.headgroup1_3.nSL=6.0012e-4 # was 6.0122e-4
        self.headgroup2_3.nSL=6.0012e-4 # was 6.0122e-4
        self.headgroup1_2.l=9.5
        self.headgroup2_2.l=9.5
        self.headgroup1_3.l=9.5
        self.headgroup2_3.l=9.5
	
        self.volacyllipid_2=925
        self.nslacyllipid_2=-2.67e-4
        self.volmethyllipid_2=98.8
        self.nslmethyllipid_2=-9.15e-5
	
        self.volacyllipid_3=925
        self.nslacyllipid_3=-2.67e-4
        self.volmethyllipid_3=98.8
        self.nslmethyllipid_3=-9.15e-5
    
        self.volchol=630
        self.nslchol=1.3215e-4    
    
        self.nf_lipid_2 = 0.
        self.nf_lipid_3 = 0.
        self.nf_chol = 0.

        self.vf_bilayer = 1.0
        self.absorb = 0.
        self.l_lipid1 = 11.
        self.l_lipid2 = 11.
        self.bulknsld = -0.56e-6
        self.normarea = 60.
        self.sigma = 2.
        self.radius_defect = 100.
        self.hc_substitution_1 = 0
        self.hc_substitution_2 = 0 
        self.global_rough = 2.0
        self.nf_tether = 0.6
        #self.nf_ihc_tether = 0.6
        self.l_tether = 3
        self.mult_tether = 5
        #self.nf_ihc_tether = 0.3
        self.substrate.l = 20
        self.substrate.z = 10
        self.substrate.nf = 1
        self.substrate.sigma1 = 2.0
        self.rho_substrate = 1

        self.fnAdjustParameters()

#define tether, bme, etc.
       
    
    #fnAdjustParameters()


# is this necessary?
    def fnInit(self, va1, na1, vm1, nm1, vh1, nh1, lh1, va2, na2, vm2, nm2, vh2, nh2, lh2, va3, na3, vm3, nm3, vh3, nh3, lh3, vc, nc):
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

            

    #printf("Enter AdjustParameters \n")
    
    def fnAdjustParameters(self):
        # printf("Enter AdjustParameters \n")
        self.fnSetSigma(self.sigma)
        
        if self.l_lipid1 <= 0:
            self.l_lipid1 = 0.01
        if self.l_lipid2 <= 0:
            self.l_lipid2 = 0.01
        if self.l_tether<=0:
            self.l_tether=0.01
        if self.nf_lipid_2<0:
            self.nf_lipid_2=0
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
        nf_ihc_tether= self.nf_tether
        nf_ihc_lipid=(1-nf_ihc_tether)*nf_ohc_lipid
        nf_ihc_lipid_2 = nf_ohc_lipid_2
        nf_ihc_lipid_3 = nf_ohc_lipid_3
        nf_ihc_chol = nf_ohc_chol

        nf_ihc_lipid_2=(1-nf_ihc_tether)*nf_ohc_lipid_2
        nf_ihc_lipid_3=(1-nf_ihc_tether)*nf_ohc_lipid_3
        nf_ihc_chol=(1-nf_ihc_tether)*nf_ohc_chol

        V_ihc = nf_ihc_lipid * (self.volacyllipid - self.volmethyllipid) + nf_ihc_lipid_2 * (self.volacyllipid_2 - self.volmethyllipid_2) + nf_ihc_lipid_3 * (self.volacyllipid_3 - self.volmethyllipid_3) + nf_ihc_chol * self.volchol + nf_ihc_tether*(self.volacyltether-self.volmethyltether)
        nSL_ihc = nf_ihc_lipid * (self.nslacyllipid - self.nslmethyllipid) + nf_ihc_lipid_2 * (self.nslacyllipid_2 - self.nslmethyllipid_2) + nf_ihc_lipid_3 * (self.nslacyllipid_3 - self.nslmethyllipid_3) + nf_ihc_chol * self.nslchol + nf_ihc_tether*(self.nslacyltether-self.nslmethyltether)
        
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
        V_im = nf_im_lipid * self.volmethyllipid + nf_im_lipid_2 * self.volmethyllipid_2 + nf_im_lipid_3 * self.volmethyllipid_3 + nf_im_tether*self.volmethyltether
        l_im = l_ihc * V_im / V_ihc
        nSL_im= nf_im_lipid*self.nslmethyllipid+nf_im_lipid_2*self.nslmethyllipid_2+nf_im_lipid_3*self.nslmethyllipid_3+ nf_im_tether*self.nslmethyltether

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
    
        
        
        #tether glycerol part -- transliteration here
        V_tg= self.tetherg.vol
    
        c_s_tg=c_s_ihc
        c_A_tg=c_A_ihc
        c_V_tg=nf_ihc_tether*(1-self.hc_substitution_2)
    
        self.tetherg.l=self.tetherg.vol/((self.volacyltether-self.volmethyltether)/self.lipid1.l)/0.9
        self.tetherg.nf=c_s_tg*c_A_tg*c_V_tg


        #tether EO part
        l_EO=self.l_tether
        V_EO=self.tether.vol
    
        c_s_EO=c_s_ihc
        c_A_EO=c_A_ihc
        c_V_EO=nf_ihc_tether*(1-self.hc_substitution_2)
        self.tether.nf=c_s_EO*c_A_EO*c_V_EO
        self.tether.l=l_EO
    
        if (self.tether.nf*self.tether.vol/self.tether.l)>self.normarea:
            self.tether.l=(self.tether.nf*self.tether.vol)/self.normarea
            #self.l_tether=self.tether.l
    
        self.l_tether= self.tether.l
        
        #bME
        self.bME.l=5.2
        l_bME= self.bME.l
        self.headgroup1.l=9.575
        V_bME= self.bME.vol
    
    
        d1= self.headgroup1.l+ self.bME.l- self.tether.l- self.tetherg.l
        if (d1>0):
            self.bME.l= self.bME.l-d1/2
            self.headgroup1.l= self.headgroup1.l-d1/2
    
    
        # look into if statements here ----
        if (( self.tether.nf* self.tether.vol/ self.tether.l+ self.mult_tether*self.tether.nf*self.bME.vol/self.bME.l)>self.normarea):
            #print(self.tether.nf, self.tether.vol, self.tether.l, self.mult_tether, self.bME.vol, self.bME.l)
            self.mult_tether=((self.normarea-self.tether.nf*self.tether.vol/self.tether.l)/(self.bME.vol/self.bME.l))/self.tether.nf
            if (self.mult_tether<0):
                self.mult_tether=0


    
        self.bME.nf=self.tether.nf*self.mult_tether #2.333
    
    
        #substrate
        self.substrate.vol=self.normarea*self.substrate.l
        self.substrate.nSL=self.rho_substrate*self.substrate.vol
    
    
        # set all lengths
        self.bME.z=0.5*self.bME.l+self.substrate.l
        self.tether.z=0.5*self.tether.l+self.substrate.l
        self.tetherg.z=self.tether.z+0.5*self.tether.l+0.5*self.tetherg.l
        self.lipid1.z= self.tetherg.z+0.5*(self.tetherg.l+self.lipid1.l)
        self.headgroup1.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1.l)
        self.headgroup1_2.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1_2.l)
        self.headgroup1_3.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * self.headgroup1_3.l)
        self.methyl1.z = self.lipid1.z + 0.5 * (self.lipid1.l + self.methyl1.l)
        self.methyl2.z = self.methyl1.z + 0.5 * (self.methyl1.l + self.methyl2.l)
        self.lipid2.z = self.methyl2.z + 0.5 * (self.methyl2.l + self.lipid2.l)
        self.headgroup2.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2.l)
        self.headgroup2_2.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2_2.l)
        self.headgroup2_3.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * self.headgroup2_3.l)
        
        
        
   
    
    #printf("nf bme %lf tether %lf tetherg %lf lipid1 %lf headgroup1 %lf headgroup1_2 %lf headgroup1_3 %lf methyl1 %lf methyl2 %lf lipid2 %lf headgroup2 %lf headgroup2_2 %lf headgroup2_3 %lf \n", bME.nf, tether.nf, tetherg.nf, lipid1.nf, headgroup1.nf, headgroup1_2.nf, headgroup1_3.nf, methyl1.nf, methyl2.nf, lipid2.nf, headgroup2.nf, headgroup2_2.nf, headgroup2_3.nf)
   # defects
        hclength = self.lipid1.l + self.methyl1.l + self.methyl2.l + self.lipid2.l
        hglength = self.headgroup1.l + self.headgroup2.l
    
        if self.radius_defect<(0.5*(hclength+hglength)):
            self.radius_defect = 0.5 * (hclength+hglength)
        # printf("defect_radius %lf hclength %lf \n",radius_defect, hclength)
    
        volhalftorus=3.14159265359*3.14159265359*(self.radius_defect-(2*hclength/3/3.14159265359))*hclength*hclength/4
        volcylinder=3.14159265359*self.radius_defect*self.radius_defect*hclength
    #printf("volhalftorus %lf volcylinder %lf \n", volhalftorus, volcylinder)
        defectarea=volhalftorus/volcylinder*(1- self.vf_bilayer)* self.normarea
    #printf("defectarea %lf \n", defectarea)
    

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
        self.defect_headgroup.nSL = defectratio * (self.headgroup2.fnGetnSL(self.bulknsld) * self.headgroup2.nf + self.headgroup2_2.fnGetnSL(self.bulknsld) * self.headgroup2_2.nf + self.headgroup2_3.fnGetnSL(self.bulknsld) * self.headgroup2_3.nf)
        self.defect_headgroup.fnSetSigma(self.sigma)
        self.defect_headgroup.nf = 1
    #printf("Exit AdjustParameters \n")


    #Return value is area at position z
    def fnGetArea(self, dz):
        result  = self.substrate.fnGetArea(dz)+ self.bME.fnGetArea(dz) + self.tether.fnGetArea(dz) + self.tetherg.fnGetArea(dz) + self.lipid1.fnGetArea(dz) + self.headgroup1.fnGetArea(dz) + self.methyl1.fnGetArea(dz)
        result += self.methyl2.fnGetArea(dz) + self.lipid2.fnGetArea(dz) + self.headgroup2.fnGetArea(dz)
        result += self.headgroup1_2.fnGetArea(dz) + self.headgroup2_2.fnGetArea(dz) + self.headgroup1_3.fnGetArea(dz)
        result += self.headgroup2_3.fnGetArea(dz) + self.defect_hydrocarbon.fnGetArea(dz) + self.defect_headgroup.fnGetArea(dz)
        return result

#get nSLD from molecular subgroups

    
   
    
    def fnGetnSLD(self, dz):
        # printf("Enter fnGetnSLD \n")
        substratearea=self.substrate.fnGetArea(dz)
        bMEarea=self.bME.fnGetArea(dz)
        tetherarea=self.tether.fnGetArea(dz)
        tethergarea=self.tetherg.fnGetArea(dz)
        lipid1area = self.lipid1.fnGetArea(dz)
        headgroup1area = self.headgroup1.fnGetArea(dz)
        headgroup1_2_area = self.headgroup1_2.fnGetArea(dz)
        headgroup1_3_area = self.headgroup1_3.fnGetArea(dz)
        methyl1area = self.methyl1.fnGetArea(dz)
        methyl2area = self.methyl2.fnGetArea(dz)
        lipid2area = self.lipid2.fnGetArea(dz)
        headgroup2area = self.headgroup2.fnGetArea(dz)
        headgroup2_2_area = self.headgroup2_2.fnGetArea(dz)
        headgroup2_3_area = self.headgroup2_3.fnGetArea(dz)
        defect_hydrocarbon_area = self.defect_hydrocarbon.fnGetArea(dz)
        defect_headgroup_area = self.defect_headgroup.fnGetArea(dz)
        
        sum  = substratearea + bMEarea + tetherarea+tethergarea + lipid1area + headgroup1area + methyl1area + methyl2area + lipid2area + headgroup2area + headgroup1_2_area
        sum += headgroup2_2_area + headgroup1_3_area + headgroup2_3_area + defect_headgroup_area + defect_hydrocarbon_area



        if sum==0:
            return 0
        else:
            result = self.substrate.fnGetnSLD(dz, self.bulknsld) * substratearea
            result += self.bME.fnGetnSLD(dz, self.bulknsld) * bMEarea
            result += self.tether.fnGetnSLD(dz, self.bulknsld) * tetherarea
            result += self.tetherg.fnGetnSLD(dz, self.bulknsld) * tethergarea
            result += self.headgroup1.fnGetnSLD(dz, self.bulknsld) * headgroup1area
            result += self.headgroup1_2.fnGetnSLD(dz, self.bulknsld) * headgroup1_2_area
            result += self.headgroup1_3.fnGetnSLD(dz, self.bulknsld) * headgroup1_3_area
            result += self.lipid1.fnGetnSLD(dz) * lipid1area
            result += self.methyl1.fnGetnSLD(dz) * methyl1area
            result += self.methyl2.fnGetnSLD(dz) * methyl2area
            result += self.lipid2.fnGetnSLD(dz) * lipid2area
            result += self.headgroup2.fnGetnSLD(dz, self.bulknsld) * headgroup2area
            result += self.headgroup2_2.fnGetnSLD(dz, self.bulknsld) * headgroup2_2_area
            result += self.headgroup2_3.fnGetnSLD(dz, self.bulknsld) * headgroup2_3_area
            result += self.defect_hydrocarbon.fnGetnSLD(dz) * defect_hydrocarbon_area
            result += self.defect_headgroup.fnGetnSLD(dz) * defect_headgroup_area
            result /= sum
            return result
# Use limits of molecular subgroups
    def fnGetLowerLimit(self):
        return self.substrate.fnGetLowerLimit()

    def fnGetUpperLimit(self):
        a = self.headgroup2.fnGetUpperLimit()
        b = self.headgroup2_2.fnGetUpperLimit()
        c = self.headgroup2_3.fnGetUpperLimit()
        return max([a, b, c])


    def fnSet(self, _sigma, _bulknsld, _global_rough, _rho_substrate, _nf_tether, _mult_tether, _l_tether, _l_lipid1, _l_lipid2, _vf_bilayer, _nf_lipid_2 = 0., _nf_lipid_3=0., _nf_chol=0., _hc_substitution_1=0., _hc_substitution_2=0., _radius_defect=100.):
        self.sigma = _sigma
        self.bulknsld = _bulknsld
        self.global_rough= _global_rough
        self.rho_substrate=_rho_substrate
        self.nf_tether=_nf_tether
        self.mult_tether=_mult_tether
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

    #printf("Exit fnSet \n")

    def fnSetSigma(self, sigma):
        self.substrate.sigma2=self.global_rough
        self.bME.sigma1=self.global_rough
        self.bME.sigma2=self.global_rough
        self.tether.sigma1=self.global_rough 
        self.tether.sigma2=sigma
        self.tetherg.fnSetSigma(sigma)
        self.headgroup1.fnSetSigma(sigma)
        self.headgroup1_2.fnSetSigma(sigma)
        self.headgroup1_3.fnSetSigma(sigma)
        self.lipid1.fnSetSigma(sigma,sigma+2)
        self.methyl1.fnSetSigma(sigma+2,sigma+2)
        self.methyl2.fnSetSigma(sigma+2,sigma+2)
        self.lipid2.fnSetSigma(sigma+2,sigma)
        self.headgroup2.fnSetSigma(sigma)
        self.headgroup2_2.fnSetSigma(sigma)
        self.headgroup2_3.fnSetSigma(sigma)
        self.defect_hydrocarbon.fnSetSigma(sigma)
        self.defect_headgroup.fnSetSigma(sigma)



    def fnWriteProfile(self, aArea, anSL, dimension, stepsize, dMaxArea):
        _, aArea, anSL = nSLDObj.fnWriteProfile(self, aArea, anSL, dimension, stepsize, dMaxArea)
        return self.normarea, aArea, anSL

    def fnWritePar2File(self, fp, cName, dimension, stepsize):
        self.substrate.fnWritePar2File(fp, "substrate", dimension, stepsize)
        self.bME.fnWritePar2File(fp, "bME", dimension, stepsize)
        self.tether.fnWritePar2File(fp, "tether", dimension, stepsize)
        self.tetherg.fnWritePar2File(fp, "tetherg", dimension, stepsize)
        self.headgroup1.fnWritePar2File(fp, "blm_headgroup1", dimension, stepsize)
        self.headgroup1_2.fnWritePar2File(fp, "blm_headgroup1_2", dimension, stepsize)
        self.headgroup1_3.fnWritePar2File(fp, "blm_headgroup1_3", dimension, stepsize)
        self.lipid1.fnWritePar2File(fp, "blm_lipid1", dimension, stepsize)
        self.methyl1.fnWritePar2File(fp, "blm_methyl1", dimension, stepsize)
        self.methyl2.fnWritePar2File(fp, "blm_methyl2", dimension, stepsize)
        self.lipid2.fnWritePar2File(fp, "blm_lipid2", dimension, stepsize)
        self.headgroup2.fnWritePar2File(fp, "blm_headgroup2", dimension, stepsize)
        self.headgroup2_2.fnWritePar2File(fp, "blm_headgroup2_2", dimension, stepsize)
        self.headgroup2_3.fnWritePar2File(fp, "blm_headgroup2_3", dimension, stepsize)
        self.defect_hydrocarbon.fnWritePar2File(fp, "blm_defect_hc", dimension, stepsize)
        self.defect_headgroup.fnWritePar2File(fp, "blm_defect_hg", dimension, stepsize)
        self.fnWriteConstant(fp, "blm_normarea", self.normarea, 0, dimension, stepsize)





class Hermite(object):
    def __init__(self, n, dstartposition, dnSLD, dnormarea):
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

    def fnSetSigma(self, sigma): 
        pass

    def fnWritePar2File(self, fp, cName, dimension, stepsize): 
        pass

    def fnGetSplineAntiDerivative(self, dz, dp, dh): 
        interval = self.fnGetSplinePars(dz, dp, dh, self.m0, self.m1, self.p0, self.p1) 
        if (0 <= interval < self.numberofcontrolpoints-1):
            dd=dp[interval+1]-dp[interval] 
            t=(dz-dp[interval])/dd 
            t_2=t*t 
            t_3=t_2*t 
            t_4=t_3*t 
            h00= (1/2)*t_4-      t_3          +t 
            h01=(-1/2)*t_4+      t_3             
            h10= (1/4)*t_4-(2/3)*t_3+(1/2)*t_2   
            h11= (1/4)*t_4-(1/3)*t_3             
            return dd*(h00*self.p0 + h10*dd*self.m0 + h01*self.p1 + h11*dd*self.m1)         
        return 0

    def fnGetSplineArea(self, dz, dp, dh, damping): 
        dampfactor=1 
        for i in range (0, self.numberofcontrolpoints):
            if (damping==0):
                self.damp[i]=dh[i] 
            else:
                self.damp[i]=dh[i]*dampfactor 
                if (dh[i] >= self.damptrigger):
                    dampfactor=dampfactor*(1/(1 + math.exp(-2.1*(dh[i]- self.dampthreshold)/self.dampFWHM))) 

        interval=self.fnGetSplinePars(self, dz, dp, self.damp, self.m0, self.m1, self.p0, self.p1) 

        if (0 <= interval < self.numberofcontrolpoints-1):
            
            dd=dp[interval+1]-dp[interval] 
            t=(dz-dp[interval])/dd 
            t_2=t*t 
            t_3=t_2*t 
            h00=   2*t_3-3*t_2+1 
            h10=     t_3-2*t_2+t 
            h01=(-2)*t_3+3*t_2 
            h11=     t_3-t_2 
            
            return h00*self.p0+h10*dd*self.m0+h01*self.p1+h11*dd*self.m1 

        return 0   

    def fnGetSplinePars(self, d, dp, dh, m0, m1, p0, p1): 
        interval=-1 
        for i in range(self.numberofcontrolpoints):
            if ((dp[i]<=self.dz) and (dp[i+1]>self.dz)):
                # printf("Found interval %i \n", i)
                interval = i
        if (self.dz==dp[-1]):
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

                    elif (interval==(self.numberofcontrolpointss-3)):
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
    
        return interval 
        
    def fnGetSplineIntegral(self, dz1, dz2, dp, dh, damping): 
        if (dz1>dz2):
            temp=dz2 
            dz2=dz1 
            dz1=temp 
        
        #check for boundaries
        if (dz1<dp[0]):
            dz1=dp[0]
        if (dz2>dp[-1]):
            dz2=dp[-1] 
        integral=0  
        d=dz1 
        while (d<=dz2):
            integral+=self.fnGetSplineArea(d,dp,dh, damping)*0.5 
            d+=0.5 
        
        return integral 

    def fnGetSplineProductIntegral(self, dz1, dz2, dp, dh1, dh2, damping1, damping2): 
        if (dz1>dz2) :
            temp=dz2 
            dz2=dz1 
            dz1=temp 
        
        #check for boundaries
        if (dz1<dp[0]):
            dz1=dp[0] 
        if (dz2>self.dp[-1]):
            dz2=self.dp[-1] 
        
        integral=0  
        d=dz1 
        while (d<=dz2):
            integral+=self.fnGetSplineArea(d,dp,dh1, damping1)*self.fnGetSplineArea(d,dp,dh2, damping2)*0.5 
             #printf("d %g damping1 %i area1 %g damping2 %i area2 %g \n", d, damping1, fnGetSplineArea(d,dp,dh1, damping1), damping2, fnGetSplineArea(d,dp,dh2, damping2)) 
            d+=0.5 
        
        return integral 
    
