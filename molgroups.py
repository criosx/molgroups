import numpy
from scipy.special import erf
from scipy.interpolate import PchipInterpolator, CubicHermiteSpline, interpn
from scipy.spatial.transform import Rotation
from scipy.ndimage.filters import gaussian_filter

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
        anSL[overmax] = anSL[overmax] * (1 - ((temparea[overmax] - dMaxArea)/aArea[overmax])) + nsl[overmax]
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
            g.fnWritePar2File(f, cName+'_'+g.name, z) # current behavior

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
# Lipids and headgroups
# ------------------------------------------------------------------------------------------------------

from periodictable.fasta import Molecule
# in formulas, use H[1] for exchangeable hydrogens, then Molecule.sld, Molecule.Dsld
# to get limiting slds, or Molecule.D2Osld(D2O_fraction=??) to get arbitrary D2O fraction.

class Component(Molecule):
    # Subclasses Molecule to automatically store a component length for use later
    # and calculate total scattering lengths
    def __init__(self, length=9.575, **kwargs):
        super().__init__(**kwargs)
        self.length = length
        self.nSL = self.sld * self.cell_volume * 1e-6 # units of Ang
        self.DnSL = self.Dsld * self.cell_volume * 1e-6 # units of Ang

def Component2Box(name=None, molecule=None):
    # Creates a new Box2Err from Component data

    box = Box2Err(name=name, dvolume=molecule.cell_volume, dlength=molecule.length)
    box.fnSetnSL(molecule.sld * molecule.cell_volume * 1e-6, molecule.Dsld * molecule.cell_volume * 1e-6)

    return box

def AddMolecules(component_list, length=None):
    # Adds together components and molecules. Note that length information is lost
    # if length is not specified
    total_formula = ' '.join([str(c.formula) for c in component_list])
    total_cell_volume = sum([c.cell_volume for c in component_list])
    total_name = ' + '.join([c.name for c in component_list])
    if length is None:
        return Molecule(name=total_name, formula=total_formula, cell_volume=total_cell_volume)
    else:
        return Component(name=total_name, formula=total_formula, cell_volume=total_cell_volume, length=length)

# PC headgroup pieces
choline = Component(name='choline', formula='C5 H13 N', cell_volume=120., length=6.34)
phosphate = Component(name='phosphate', formula='PO4', cell_volume=54., length=3.86)
carbonyl_glycerol = Component(name='carbonyl + glycerol', formula='C5 O4 H5', cell_volume=147., length=4.21)

# standard headgroups
pc = Component(name='PC', formula='C10 H18 O8 N P', cell_volume=331.00, length=9.575)
pe = Component(name='PE', formula='C7 H9 H[1]3 O8 N P', cell_volume=262., length=7.7)
pg = Component(name='PG', formula='C8 H10 H[1]2 O10 P', cell_volume=240., length=7.8)
ps = Component(name='PS', formula='C8 H8 H[1]3 N O10 P', cell_volume=280., length=8.1)
pa = Component(name='PA', formula='C5 H5 H[1] O8 P', cell_volume=174., length=5.0)
pi = Component(name='PI', formula='C11 H7 H[1]5 O13 P', cell_volume=370.7, length=10.7)
pip2 = Component(name='PI(4,5)P2', formula='C11 H7 H[1]5 O19 P3', cell_volume=500., length=12.0) # doesn't match molgroups.cc
cardiolipin = Component(name='cardiolipin', formula='C13 H15 H[1] O17 P2', cell_volume=684.4, length=9.56)

# standard acyl chains
oleoyl = Molecule(name='oleoyl', formula='C17 H33', cell_volume=972./2.0)
palmitoyl = Molecule(name='palmitoyl', formula='C15 H31', cell_volume=770./2.0)
myristoyl = Molecule(name='myristoyl', formula='C13 H27', cell_volume=770./2.0)
phytanoyl = Molecule(name='phytanoyl', formula='C19 H39', cell_volume=1095./2.0)
cholesterol = Molecule(name='cholesterol', formula='C27 H45 H[1] O', cell_volume=630.)

# methyl
methyl = Molecule(name='methyl', formula='CH3', cell_volume=98.8/2.0)
Dmethyl = Molecule(name='Dmethyl', formula='CD3', cell_volume=98.8/2.0)

# Tether components
SAc = Molecule(name='thiol acetate', formula = 'C2H3OS', cell_volume=117.0)
EO6 = Molecule(name='6x ethylene oxide', formula='(C2H4O)6', cell_volume=360.0)
tetherg_ether = Molecule(name='tether glycerol ether', formula='C5H9O2', cell_volume=125.40)
tetherg_ester = carbonyl_glycerol
ethanoyl = Molecule(name='ethanoyl', formula='C2H5O', cell_volume=(117 - 25.75))
thiol = Molecule(name='sulfur', formula='S', cell_volume=25.75)
bmeSAc = AddMolecules([SAc, ethanoyl])
bme = AddMolecules([thiol, ethanoyl])
SEO6 = AddMolecules([thiol, EO6])
SAcEO6 = AddMolecules([SAc, EO6])

# TODO: A headgroup class must contain an "innerleaflet" flag that determines whether the headgroup
# is in the inner or outer leaflet. The profiles so obtained should be flipped. This should perhaps
# be standardized, (1) by creating a CompositeHeadgroup class that has this flag already, and (2) by
# finding a reasonable way of flipping the profile after calculating it (removes lots of if statements)

class PC(CompositenSLDObj):
    def __init__(self, name='PC', innerleaflet=False, **kwargs):
        # innerleaflet flag locates it in the inner leaflet and flips the order of cg, phosphate,
        # and choline groups. If False, it's the outer leaflet
        super().__init__(name=name, **kwargs)
        self.cg = Component2Box(name='cg', molecule=carbonyl_glycerol)
        self.phosphate = Component2Box(name='phosphate', molecule=phosphate)
        self.choline = Component2Box(name='choline', molecule=choline)
        
        self.innerleaflet = innerleaflet

        self.groups = {"cg": self.cg, "phosphate": self.phosphate, "choline": self.choline}

        if innerleaflet:

            self.cg.sigma2, self.cg.sigma1 = 2.53, 2.29
            self.phosphate.sigma2, self.phosphate.sigma1 = 2.29, 2.02
            self.choline.sigma2, self.choline.sigma1 = 2.02, 2.26

        else:
            self.cg.sigma1, self.cg.sigma2 = 2.53, 2.29
            self.phosphate.sigma1, self.phosphate.sigma2 = 2.29, 2.02
            self.choline.sigma1, self.choline.sigma2 = 2.02, 2.26

        self.cg.nf=1 
        self.phosphate.nf=1 
        self.choline.nf=1

        self.l = 9.575
        self.vol=self.cg.vol+self.phosphate.vol+self.choline.vol
        self.nSL=self.cg.nSL+self.phosphate.nSL+self.choline.nSL
        self.ph_relative_pos = .58
        self.nf = 1.
        self.fnAdjustParameters()
        self.fnFindSubgroups()

    def fnAdjustParameters(self):
        if self.innerleaflet:
            self.cg.z = self.z + 0.5 * self.l - 0.5 * self.cg.l
            self.choline.z = self.z - 0.5 * self.l + 0.5 * self.choline.l
            z0 = self.z - 0.5 * self.l + 0.5 * self.phosphate.l
            z1 = self.z + 0.5 * self.l - 0.5 * self.phosphate.l
            self.phosphate.z = z0 + (z1 - z0) * (1 - self.ph_relative_pos)
        else:
            self.cg.z = self.z - 0.5 * self.l + 0.5*self.cg.l
            self.choline.z = self.z + 0.5 * self.l - 0.5*self.choline.l
            z0 = self.z - 0.5*self.l + 0.5 * self.phosphate.l
            z1 = self.z + 0.5 * self.l - 0.5*self.phosphate.l
            self.phosphate.z = z0 + (z1 - z0) * self.ph_relative_pos
    
    def fnSet(self, l=9.575, ph_relative_pos=.5, cg_nSL=1.885e-3, ch_nSL=1.407e-3, ph_nSL=1.323e-3):
        # TODO: these overwrite the known nSL parameters. Do we really want to do that here?
        self.cg.nSL = cg_nSL
        self.choline.nSL = ch_nSL
        self.phosphate.nSL = ph_nSL
        self.l = l
        self.ph_relative_pos=ph_relative_pos
        self.fnAdjustParameters()
    
    def fnGetLowerLimit(self):
        return self.choline.fnGetLowerLimit() if self.innerleaflet else self.cg.fnGetLowerLimit()
    
    def fnGetUpperLimit(self):
        return self.cg.fnGetUpperLimit() if self.innerleaflet else self.choline.fnGetUpperLimit()
    
    def fnGetnSL(self, bulknsld=None):
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
        prefix = "PCm" if self.innerleaflet else "PC"
        fp.write(prefix + " "+cName+" z "+str(self.z)+" l "+str(self.l)+" vol " +
                 str(self.cg.vol + self.phosphate.vol + self.choline.vol)+" nf " + str(self.nf)+" \n")
        self.fnWriteData2File(fp, cName, z)

        super().fnWritePar2File(fp, cName, z)

class PCm(PC):
    # deprecated. Use PC(innerleaflet=True)
    def __init__ (self, **kwargs):
        super().__init__(**kwargs)
        self.cg.sigma2 = 2.53
        self.cg.sigma1 = 2.29
        self.phosphate.sigma2 = 2.29
        self.phosphate.sigma1 = 2.02
        self.choline.sigma2 = 2.02
        self.choline.sigma1 = 2.26
        self.fnAdjustParameters()
        
    def fnAdjustParameters(self):
        self.cg.z = self.z + 0.5 * self.l - 0.5 * self.cg.l
        self.choline.z = self.z - 0.5 * self.l + 0.5 * self.choline.l
        z0 = self.z - 0.5 * self.l + 0.5 * self.phosphate.l
        z1 = self.z + 0.5 * self.l - 0.5 * self.phosphate.l
        self.phosphate.z = z0 + (z1 - z0) * (1 - self.ph_relative_pos)

    def fnGetLowerLimit(self):
        return self.choline.fnGetLowerLimit()

    def fnGetUpperLimit(self): 
        return self.cg.fnGetUpperLimit()

    def fnWritePar2File (self,fp, cName, z):
        fp.write("PCm " + cName + " z " + str(self.z) + " l " + str(self.l) + " vol " +
                 str(self.cg.vol + self.phosphate.vol + self.choline.vol) + " nf " + str(self.nf) + " \n")
        self.fnWriteData2File(fp, cName, z)

        super().fnWritePar2File(fp, cName, z)

# TODO: make a Tether class that can be passed to tethered bilayers. Perhaps make Tether class None for supported bilayers?
# TODO: think about making a generic bilayer class with "substrate=", "tether=", and "filler=" inputs. If substrate is None, bilayer floats.
#       If tether is None but substrate exists, it's a ssBLM. "filler" is some material that the substrate is decorated with like bME.

class Lipid(object):
    def __init__(self, name=None, headgroup=PC, tails=[oleoyl, oleoyl], methyls=methyl):
        # hg = Molecule object with headgroup information or PC molgroups object
        # tails = List of molecule objects containing lipid tail information
        # methyls = List of methyl groups, one for each lipid tail; OR, a single group
        #           that is copied to each tail; OR, None, which uses a null group with
        #           no volume (use with cholesterol). Use methyl=Dmethyl for CD3.
        # default is DOPC using the PC class.
        
        n_tails = len(tails)

        # tails section
        if isinstance(tails, list):
            tail_volume = sum([t.cell_volume for t in tails])
            tail_formula = ' '.join([str(t.formula) for t in tails])
            self.tails = Molecule(name='tails',formula=tail_formula, cell_volume=tail_volume)
        elif isinstance(tails, Molecule):
            self.tails = tails
        else:
            raise TypeError('Lipid.tails must be Molecule or list of Molecules')
        
        # headgroup
        self.headgroup = headgroup if headgroup is not None else Component(name=None, formula = '', cell_volume=0.0, length=9.575)

        # Create methyl groups
        if isinstance(methyls, list):
            assert n_tails==len(methyls), 'Lipid tails and lipid methyl lists must have equal length, not %i and %i' % (len(tails), len(methyls))
        elif isinstance(methyls, Molecule) | (methyls is None):
            methyls = n_tails * [methyls]
        else:
            raise TypeError('Lipid.methyls must be Molecule or list of Molecules')
        methyls = [m if m is not None else Molecule(name='', formula='', cell_volume=0.0) for m in methyls]  # Replace None with null molecule
        methylvolsum = sum([m.cell_volume for m in methyls])
        methylformula = ' '.join([str(m.formula) for m in methyls])
        self.methyls = Molecule(name='methyls', formula = methylformula, cell_volume=methylvolsum)

        if name is not None:
            self.name = name
        else:
            self.name = ' + '.join([c.name for c in (tails + [self.headgroup]) if c.name is not None])

class Tether(Lipid):
    """Subclass of Lipid for use with tether molecules.
        Uses hg field for tether glycerol.
        Tether includes both volume from """
    def __init__(self, name=None, tether=SEO6, tetherg=tetherg_ether, tails=[oleoyl, oleoyl], methyls=methyl):
        super().__init__(name=name, headgroup=tetherg, tails=tails, methyls=methyls)
        self.tether = tether
        self.tetherg = self.headgroup

        if name is not None:
            self.name = name
        else:
            self.name = ' + '.join([c.name for c in ([self.tether] + tails) if c.name is not None])



DOPC = Lipid(name='DOPC', headgroup=PC, tails=2*[oleoyl])
POPC = Lipid(name='POPC', headgroup=PC, tails=[palmitoyl, oleoyl])
DOPS = Lipid(name='DOPS', headgroup=ps, tails=2*[oleoyl])
chol = Lipid(name='chol', headgroup=None, tails=[cholesterol], methyls=None)

WC14 = Tether(name='WC14', tether=SEO6, tetherg=tetherg_ether, tails=[myristoyl, myristoyl])
HC18 = Tether(name='HC18', tether=SEO6, tetherg=tetherg_ether, tails=[oleoyl, oleoyl])
HC18SAc = Tether(name='HC18SAc', tether=SAcEO6, tetherg=tetherg_ether, tails=[oleoyl, oleoyl])

# ------------------------------------------------------------------------------------------------------
# Lipid Bilayer
# ------------------------------------------------------------------------------------------------------

def _unpack_lipids(self, lipids):
    """ Helper function for BLM classes that unpacks a lipid list into headgroup objects
        and lists of acyl chain and methyl volumes and nSLs. Creates the following attributes:
        o headgroups1: a list of inner leaflet headgroup objects
        o headgroups2: a list of outer leaflet headgroup objects
        NB: each object in headgroups1 and headgroups2 is created as a standalone attribute in "self"
            so that searching for all nSLDObj in self will return each headgroup individually as 
            headgroup1_1, headgroup1_2, ... headgroup1_n for n lipids. Similarly for headgroup2.
        o vol_acyl_lipids: a numpy array of acyl chain volumes
        o vol_methyl_lipids: a numpy array of methyl volumes
        o nsl_acyl_lipids: a numpy array of acyl chain nSL
        o nsl_methyl_lipids: a numpy array of methyl nSL
        """
    n_lipids = len(lipids)
    self.lipids = lipids # store this information

    self.headgroups1 = []   # list of inner headgroups
    self.headgroups2 = []   # list of outer headgroups
    self.vol_acyl_lipids = numpy.zeros(n_lipids)     # list of acyl chain volumes
    self.nsl_acyl_lipids = numpy.zeros(n_lipids)     # list of acyl chain nsls
    self.vol_methyl_lipids = numpy.zeros(n_lipids)     # list of acyl chain methyl volumes
    self.nsl_methyl_lipids = numpy.zeros(n_lipids)     # list of acyl chain methyl nsls

    # create objects for each lipid headgroup and a composite tail object
    # TODO: create separate lipid objects so methyl sigmas can be different for different species
    for i, lipid in enumerate(lipids):
        ihg_name = 'headgroup1_%i' % (i+1)
        ohg_name = 'headgroup2_%i' % (i+1)

        if isinstance(lipid.headgroup, Component):
            # populates nSL, nSL2, vol, and l
            ihg_obj = Component2Box(name=ihg_name, molecule=lipid.hg)
            ohg_obj = Component2Box(name=ohg_name, molecule=lipid.hg)

        elif issubclass(lipid.headgroup, CompositenSLDObj):
            ihg_obj = lipid.headgroup(name=ihg_name, innerleaflet=True)
            ohg_obj = lipid.headgroup(name=ohg_name, innerleaflet=False)

        else:
            raise TypeError('Lipid.hg must be a Headgroup object or a subclass of CompositenSLDObj')

        # Note that because there's always one lipid, the objects self.headgroups1[0] and self.headgroups2[0] always exist
        self.__setattr__(ihg_name, ihg_obj)
        self.headgroups1.append(self.__getattribute__(ihg_name))
        self.__setattr__(ohg_name, ohg_obj)
        self.headgroups2.append(self.__getattribute__(ohg_name))

        # find null headgroups to exclude from averaging over headgroup properties
        self.null_hg1 = numpy.array([hg.vol <= 0.0 for hg in self.headgroups1], dtype=bool)
        self.null_hg2 = numpy.array([hg.vol <= 0.0 for hg in self.headgroups2], dtype=bool)

        self.vol_acyl_lipids[i] = lipid.tails.cell_volume
        self.nsl_acyl_lipids[i] = lipid.tails.sld * lipid.tails.cell_volume * 1e-6

        self.vol_methyl_lipids[i] = lipid.methyls.cell_volume
        self.nsl_methyl_lipids[i] = lipid.methyls.sld * lipid.methyls.cell_volume * 1e-6

class BLM_arbitrary(CompositenSLDObj):
    def __init__(self, lipids=[DOPC], lipid_nf=[1.0], **kwargs):
        super().__init__(**kwargs)
        assert len(lipids)==len(lipid_nf), 'List of lipids and number fractions must be of equal length, not %i and %i' % (len(lipids), len(lipid_nf))
        assert len(lipids)>0, 'Must specify at least one lipid'

        # normalize number fractions. This allows ratios of lipids to be given instead of number fractions
        self.lipid_nf = numpy.array(lipid_nf) / numpy.sum(lipid_nf)

        _unpack_lipids(self, lipids)

        self.lipid1 = Box2Err(name='lipid1')
        self.methyl1 = Box2Err(name='methyl1')
        self.methyl2 = Box2Err(name='methyl2')
        self.lipid2 = Box2Err(name='lipid2')

        self.defect_hydrocarbon = Box2Err(name='defect_hc')
        self.defect_headgroup = Box2Err(name='defect_hg')
        
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
        
    def fnAdjustParameters(self):

        self._adjust_lipids()
        self._adjust_z(self.startz + self.av_hg1_l)
        self._adjust_defects()
        self.fnSetSigma(self.sigma)

    def _adjust_lipids(self):

        self.l_lipid1 = max(self.l_lipid1, 0.01)
        self.l_lipid2 = max(self.l_lipid2, 0.01)
        for nf in self.lipid_nf:
            nf = max(nf, 0)

        # this is changed normalization behavior, but it's a bit more intuitive
        self.lipid_nf /= numpy.sum(self.lipid_nf)

        self.vf_bilayer = max(self.vf_bilayer, 1E-5)
        self.vf_bilayer = min(self.vf_bilayer, 1)

        # calculate average headgroup lengths, ignore zero volume (e.g. cholesterol)
        self.av_hg1_l = numpy.sum(numpy.array([hg.l for hg, use in zip(self.headgroups1, ~self.null_hg1) if use]) * self.lipid_nf[~self.null_hg1]) / numpy.sum(self.lipid_nf[~self.null_hg1])
        self.av_hg2_l = numpy.sum(numpy.array([hg.l for hg, use in zip(self.headgroups2, ~self.null_hg2) if use]) * self.lipid_nf[~self.null_hg2]) / numpy.sum(self.lipid_nf[~self.null_hg2])
        
        # outer hydrocarbons
        l_ohc = self.l_lipid2
        nf_ohc_lipid = self.lipid_nf
        V_ohc = numpy.sum(nf_ohc_lipid * (self.vol_acyl_lipids - self.vol_methyl_lipids))
        nSL_ohc = numpy.sum(nf_ohc_lipid * (self.nsl_acyl_lipids - self.nsl_methyl_lipids))
        #V_ohc = nf_ohc_lipid * (self.volacyllipid - self.volmethyllipid) + nf_ohc_lipid_2 * (self.volacyllipid_2 - self.volmethyllipid_2) + nf_ohc_lipid_3 * (self.volacyllipid_3 - self.volmethyllipid_3) + nf_ohc_chol * self.volchol
        #nSL_ohc = nf_ohc_lipid * (self.nslacyllipid - self.nslmethyllipid) + nf_ohc_lipid_2 * (self.nslacyllipid_2 - self.nslmethyllipid_2) + nf_ohc_lipid_3 * (self.nslacyllipid_3 - self.nslmethyllipid_3) + nf_ohc_chol * self.nslchol
        
        self.normarea = V_ohc / l_ohc
        c_s_ohc = self.vf_bilayer
        c_A_ohc = 1
        c_V_ohc = 1
        
        self.lipid2.l = l_ohc
        self.lipid2.vol = V_ohc
        self.lipid2.nSL = nSL_ohc
        self.lipid2.nf = c_s_ohc * c_A_ohc * c_V_ohc
        
        # outer methyl
        #nf_om_lipid = nf_ohc_lipid
        #nf_om_lipid_2 = nf_ohc_lipid_2
        #nf_om_lipid_3 = nf_ohc_lipid_3
        nf_om_lipid = nf_ohc_lipid
        #V_om = nf_om_lipid * self.volmethyllipid + nf_om_lipid_2 * self.volmethyllipid_2 + nf_om_lipid_3 * self.volmethyllipid_3
        V_om = numpy.sum(nf_om_lipid * self.vol_methyl_lipids)   # cholesterol has a methyl volume of zero and does not contribute
        l_om = l_ohc * V_om / V_ohc
        #nSL_om = nf_om_lipid * self.nslmethyllipid + nf_om_lipid_2 * self.nslmethyllipid_2 + nf_om_lipid_3 * self.nslmethyllipid_3
        nSL_om = numpy.sum(nf_om_lipid * self.nsl_methyl_lipids)
        
        c_s_om = c_s_ohc
        c_A_om = 1
        c_V_om = 1
        
        self.methyl2.l = l_om
        self.methyl2.vol = V_om
        self.methyl2.nSL = nSL_om
        self.methyl2.nf = c_s_om * c_A_om * c_V_om
        
        # inner hydrocarbons
        l_ihc = self.l_lipid1
        
        #nf_ihc_lipid = nf_ohc_lipid
        #nf_ihc_lipid_2 = nf_ohc_lipid_2
        #nf_ihc_lipid_3 = nf_ohc_lipid_3
        #nf_ihc_chol = nf_ohc_chol
        nf_ihc_lipid = nf_ohc_lipid
        V_ihc = numpy.sum(nf_ihc_lipid * (self.vol_acyl_lipids - self.vol_methyl_lipids))
        #V_ihc = nf_ihc_lipid * (self.volacyllipid - self.volmethyllipid) + nf_ihc_lipid_2 * (self.volacyllipid_2 - self.volmethyllipid_2) + nf_ihc_lipid_3 * (self.volacyllipid_3 - self.volmethyllipid_3) + nf_ihc_chol * self.volchol
        #nSL_ihc = nf_ihc_lipid * (self.nslacyllipid - self.nslmethyllipid) + nf_ihc_lipid_2 * (self.nslacyllipid_2 - self.nslmethyllipid_2) + nf_ihc_lipid_3 * (self.nslacyllipid_3 - self.nslmethyllipid_3) + nf_ihc_chol * self.nslchol
        nSL_ihc = numpy.sum(nf_ihc_lipid * (self.nsl_acyl_lipids - self.nsl_methyl_lipids))
        
        c_s_ihc = self.vf_bilayer
        c_A_ihc = self.normarea * l_ihc / V_ihc
        c_V_ihc = 1
        
        self.lipid1.l = l_ihc
        self.lipid1.vol = V_ihc
        self.lipid1.nSL = nSL_ihc
        self.lipid1.nf = c_s_ihc * c_A_ihc * c_V_ihc
        
        # inner methyl
        #nf_im_lipid = nf_ihc_lipid
        #nf_im_lipid_2 = nf_ihc_lipid_2
        #nf_im_lipid_3 = nf_ihc_lipid_3
        nf_im_lipid = nf_ihc_lipid
        V_im = numpy.sum(nf_im_lipid * self.vol_methyl_lipids)   # cholesterol has a methyl volume of zero and does not contribute
        #V_im = nf_im_lipid * self.volmethyllipid + nf_im_lipid_2 * self.volmethyllipid_2 + nf_im_lipid_3 * self.volmethyllipid_3
        l_im = l_ihc * V_im / V_ihc
        nSL_im = numpy.sum(nf_im_lipid * self.nsl_methyl_lipids)
        #nSL_im = nf_im_lipid * self.nslmethyllipid + nf_im_lipid_2 * self.nslmethyllipid_2 + nf_im_lipid_3 * self.nslmethyllipid_3
        
        c_s_im = c_s_ihc
        c_A_im = c_A_ihc
        c_V_im = 1
        
        self.methyl1.l = l_im
        self.methyl1.vol = V_im
        self.methyl1.nSL = nSL_im
        self.methyl1.nf = c_s_im * c_A_im * c_V_im

        # headgroup size and position
        for hg1, nf_ihc, hg2, nf_ohc in zip(self.headgroups1, nf_ihc_lipid, self.headgroups2, nf_ohc_lipid):
            hg2.nf = c_s_ohc * c_A_ohc * nf_ohc * (1 - self.hc_substitution_2)
            hg1.nf = c_s_ihc * c_A_ihc * nf_ihc * (1 - self.hc_substitution_1)
            if hasattr(hg1, 'fnAdjustParameters'):
                hg1.fnAdjustParameters()
            if hasattr(hg2, 'fnAdjustParameters'):
                hg2.fnAdjustParameters()

    def _adjust_defects(self):# defects
        hclength = self.lipid1.l + self.methyl1.l + self.methyl2.l + self.lipid2.l
        hglength = self.av_hg1_l + self.av_hg2_l

        if self.radius_defect<(0.5*(hclength+hglength)):
            self.radius_defect = 0.5 * (hclength+hglength)
        
        volhalftorus = numpy.pi**2 * (self.radius_defect - (2. * hclength / 3. / numpy.pi)) * hclength * hclength / 4
        volcylinder = numpy.pi * self.radius_defect * self.radius_defect * hclength
        defectarea = volhalftorus / volcylinder * (1 - self.vf_bilayer) * self.normarea
        
        self.defect_hydrocarbon.vol = defectarea * hclength
        self.defect_hydrocarbon.l = hclength
        self.defect_hydrocarbon.z = self.lipid1.z - 0.5 * self.lipid1.l + 0.5 * hclength
        self.defect_hydrocarbon.nSL = self.lipid2.nSL / self.lipid2.vol * self.defect_hydrocarbon.vol
        self.defect_hydrocarbon.fnSetSigma(self.sigma)
        self.defect_hydrocarbon.nf = 1
        
        defectratio = self.defect_hydrocarbon.vol / self.lipid2.vol
        #self.defect_headgroup.vol = defectratio * (self.headgroup2.vol * self.headgroup2.nf + self.headgroup2_2.vol * self.headgroup2_2.nf + self.headgroup2_3.vol * self.headgroup2_3.nf)
        self.defect_headgroup.vol = defectratio * numpy.sum([hg.nf * hg.vol for hg in self.headgroups2])
        self.defect_headgroup.l = hclength + hglength
        self.defect_headgroup.z = self.headgroups1[0].fnGetZ() - 0.5 * self.headgroups1[0].l + 0.5 * (hclength + hglength)
        self.defect_headgroup.nSL = defectratio * numpy.sum([hg.nf * hg.fnGetnSL(self.bulknsld) for hg in self.headgroups2])
        #self.defect_headgroup.nSL = defectratio * (self.headgroup2.fnGetnSL()*self.headgroup2.nf * self.headgroup2.nf + self.headgroup2_2.fnGetnSL(self.bulknsld) * self.headgroup2_2.nf + self.headgroup2_3.fnGetnSL(self.bulknsld) * self.headgroup2_3.nf)
        self.defect_headgroup.fnSetSigma(self.sigma)
        self.defect_headgroup.nf = 1

    def _adjust_z(self, startz):
        # startz is the position of the hg1/lipid1 interface.
        self.lipid1.z= startz + 0.5 * self.lipid1.l # change here: use average headgroup length instead of a specific headgroup
        self.methyl1.z = self.lipid1.z + 0.5 * (self.lipid1.l + self.methyl1.l)
        self.methyl2.z = self.methyl1.z + 0.5 * (self.methyl1.l + self.methyl2.l)
        self.lipid2.z = self.methyl2.z + 0.5 * (self.methyl2.l + self.lipid2.l)

        for hg1, hg2 in zip(self.headgroups1, self.headgroups2):
            hg1.fnSetZ(self.lipid1.z - 0.5 * self.lipid1.l - 0.5 * hg1.l)
            hg2.fnSetZ(self.lipid2.z + 0.5 * self.lipid2.l + 0.5 * hg2.l)

    # return value of center of the membrane
    def fnGetCenter(self):
        return self.methyl1.z + 0.5 * self.methyl1.l

    # Use limits of molecular subgroups
    def fnGetLowerLimit(self):
        return min([hg.fnGetLowerLimit() for hg in self.headgroups1])
        
    def fnGetUpperLimit(self):
        return max([hg.fnGetUpperLimit() for hg in self.headgroups2])

    def fnSet(self, sigma, bulknsld, startz, l_lipid1, l_lipid2, vf_bilayer,
              nf_lipids=None, hc_substitution_1=0.,
              hc_substitution_2=0., radius_defect=100.):
        self.sigma = sigma
        self.bulknsld = bulknsld
        self.startz = startz
        self.l_lipid1 = l_lipid1
        self.l_lipid2 = l_lipid2
        self.vf_bilayer = vf_bilayer
        if nf_lipids is not None:   # pass None to keep lipid_nf unchanged
            assert len(nf_lipids)==len(self.lipids), 'nf_lipids must be same length as number of lipids in bilayer, not %i and %i' % (len(nf_lipids), len(self.lipids))
            self.lipid_nf = numpy.array(nf_lipids)
        self.hc_substitution_1 = hc_substitution_1
        self.hc_substitution_2 = hc_substitution_2
        self.radius_defect = radius_defect
        
        self.fnAdjustParameters()

    def fnSetSigma(self, sigma):
        for hg1, hg2 in zip(self.headgroups1, self.headgroups2):
            hg1.fnSetSigma(sigma)
            hg2.fnSetSigma(sigma)
        self.lipid1.fnSetSigma(sigma,numpy.sqrt(sigma**2+self.methyl_sigma**2))
        self.methyl1.fnSetSigma(numpy.sqrt(sigma**2+self.methyl_sigma**2), numpy.sqrt(sigma**2+self.methyl_sigma**2))
        self.methyl2.fnSetSigma(numpy.sqrt(sigma**2+self.methyl_sigma**2), numpy.sqrt(sigma**2+self.methyl_sigma**2))
        self.lipid2.fnSetSigma(numpy.sqrt(sigma**2+self.methyl_sigma**2),sigma)
        self.defect_hydrocarbon.fnSetSigma(sigma)
        self.defect_headgroup.fnSetSigma(sigma)

    def fnWritePar2File(self, fp, cName, z):
        super().fnWritePar2File(fp, cName, z)
        self.fnWriteConstant(fp, cName+"_normarea", self.normarea, 0, z)

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

        self.volacyllipid = 972
        self.nslacyllipid = -2.67e-4
        self.volmethyllipid = 98.8
        self.nslmethyllipid = -9.15e-5
        self.volacyllipid_2 = 972
        self.nslacyllipid_2 = -2.67e-4
        self.volmethyllipid_2 = 98.8
        self.nslmethyllipid_2 = -9.15e-5
        self.volacyllipid_3 = 972
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
        self.fnSetSigma(self.sigma)
        
        self.l_lipid1 = max(self.l_lipid1, 0.01)
        self.l_lipid2 = max(self.l_lipid2, 0.01)
        self.nf_lipid_2 = max(self.nf_lipid_2, 0)
        self.nf_lipid_3 = max(self.nf_lipid_3, 0)
        self.nf_chol = max(self.nf_chol, 0)
        sum = self.nf_lipid_2 + self.nf_lipid_3 + self.nf_chol
        if sum > 1:
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
        
        self.lipid2.l = l_ohc
        self.lipid2.vol = V_ohc
        self.lipid2.nSL = nSL_ohc
        self.lipid2.nf = c_s_ohc * c_A_ohc * c_V_ohc
        
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
        self.headgroup1.nf   = c_s_ihc * c_A_ihc * nf_ihc_lipid * (1 - self.hc_substitution_1)
        self.headgroup1_2.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid_2 * (1 - self.hc_substitution_1)
        self.headgroup1_3.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid_3 * (1 - self.hc_substitution_1)

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
        
        # defects
        hclength = self.lipid1.l + self.methyl1.l + self.methyl2.l + self.lipid2.l
        hglength = self.headgroup1.l + self.headgroup2.l
        
        if self.radius_defect<(0.5*(hclength+hglength)):
            self.radius_defect = 0.5 * (hclength+hglength)
        
        volhalftorus = 3.14159265359 * 3.14159265359 * (self.radius_defect - (2 * hclength / 3 / 3.14159265359)) * hclength * hclength / 4
        volcylinder = 3.14159265359 * self.radius_defect * self.radius_defect * hclength
        defectarea = volhalftorus / volcylinder * (1 - self.vf_bilayer) * self.normarea
        
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
        self.fnWriteConstant(fp, cName+"_normarea", self.normarea, 0, z)

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

        self.lipid2.l = l_ohc
        self.lipid2.vol = V_ohc
        self.lipid2.nSL = nSL_ohc
        self.lipid2.nf = c_s_ohc * c_A_ohc * c_V_ohc

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
        self.headgroup1.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid * (1 - self.hc_substitution_1)
        self.headgroup1_2.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid_2 * (1 - self.hc_substitution_1)
        self.headgroup1_3.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid_3 * (1 - self.hc_substitution_1)

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

        # defects
        hclength = self.lipid1.l + self.methyl1.l + self.methyl2.l + self.lipid2.l
        hglength = self.headgroup1.l + self.headgroup2.l

        if self.radius_defect < (0.5 * (hclength + hglength)):
            self.radius_defect = 0.5 * (hclength + hglength)

        volhalftorus = 3.14159265359 * 3.14159265359 * (
                    self.radius_defect - (2 * hclength / 3 / 3.14159265359)) * hclength * hclength / 4
        volcylinder = 3.14159265359 * self.radius_defect * self.radius_defect * hclength
        defectarea = volhalftorus / volcylinder * (1 - self.vf_bilayer) * self.normarea

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


    # Use limits of molecular subgroups
    def fnGetLowerLimit(self):
        return self.substrate.fnGetLowerLimit()

    def fnGetUpperLimit(self):
        a = self.headgroup2.fnGetUpperLimit()
        b = self.headgroup2_2.fnGetUpperLimit()
        c = self.headgroup2_3.fnGetUpperLimit()
        return max([a, b, c])

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

class ssBLM_arbitrary(BLM_arbitrary):
    def __init__(self, **kwargs):

        # add ssBLM-specific subgroups
        self.substrate = Box2Err(name='substrate')
        self.siox = Box2Err(name='siox')

        self.substrate.l = 40
        self.substrate.z = 0
        self.substrate.nf = 1
        self.rho_substrate = 2.07e-6
        self.l_siox = 1
        self.rho_siox = 1

        self.siox.l = 20
        self.siox.z = self.substrate.z + self.substrate.l / 2.0 + self.siox.l / 2.0
        self.siox.nf = 1

        self.l_submembrane = 10.
        self.global_rough = 2.0

        super().__init__(**kwargs)

    def fnAdjustParameters(self):
        
        self._adjust_lipids()

        # substrate
        self.substrate.vol = self.normarea * self.substrate.l
        self.substrate.nSL = self.rho_substrate * self.substrate.vol
        self.siox.l = self.l_siox
        self.siox.vol = self.normarea * self.siox.l
        self.siox.nSL = self.rho_siox * self.siox.vol
        self.siox.z=self.substrate.z + self.substrate.l / 2 + 0.5 * self.siox.l

        self._adjust_z(self.substrate.z + self.substrate.l / 2 + self.siox.l + self.l_submembrane + self.av_hg1_l)
        self._adjust_defects()

        self.fnSetSigma(self.sigma)

    def fnSetSigma(self, sigma):
        super().fnSetSigma(sigma)
        self.substrate.fnSetSigma(self.global_rough)
        self.siox.fnSetSigma(self.global_rough)

    def fnGetLowerLimit(self):
        return self.substrate.fnGetLowerLimit()
        # does this make sense since this goes to negative z and isn't intended to be used?

    def fnSet(self, _sigma, _bulknsld, _global_rough, _rho_substrate, _rho_siox, _l_siox, _l_submembrane, _l_lipid1,
              _l_lipid2, _vf_bilayer, _nf_lipids=None, _hc_substitution_1=0.,
              _hc_substitution_2=0., _radius_defect=100.):
        self.sigma = _sigma
        self.bulknsld = _bulknsld
        self.global_rough = _global_rough
        self.rho_substrate = _rho_substrate
        self.rho_siox = _rho_siox
        self.l_siox = _l_siox 
        self.l_submembrane = _l_submembrane 
        self.l_lipid1 = _l_lipid1
        self.l_lipid2 = _l_lipid2
        self.vf_bilayer = _vf_bilayer
        if _nf_lipids is not None:   # pass None to keep lipid_nf unchanged
            assert len(_nf_lipids)==len(self.lipids), 'nf_lipids must be same length as number of lipids in bilayer, not %i and %i' % (len(_nf_lipids), len(self.lipids))
            self.lipid_nf = numpy.array(_nf_lipids)
        self.hc_substitution_1 = _hc_substitution_1
        self.hc_substitution_2 = _hc_substitution_2
        self.radius_defect = _radius_defect

        self.fnAdjustParameters()


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
        self.headgroup1 = PC(name='headgroup1', innerleaflet=True)  # mirrored PC head group
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

        self.lipid2.l = l_ohc
        self.lipid2.vol = V_ohc
        self.lipid2.nSL = nSL_ohc
        self.lipid2.nf = c_s_ohc * c_A_ohc * c_V_ohc

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
        self.headgroup1.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid * (1 - self.hc_substitution_1)
        self.headgroup1_2.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid_2 * (1 - self.hc_substitution_1)
        self.headgroup1_3.nf = c_s_ihc * c_A_ihc * nf_ihc_lipid_3 * (1 - self.hc_substitution_1)

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
        self.bME.z = 0.5 * self.bME.l + self.substrate.z + self.substrate.l * 0.5
        self.tether.z = 0.5 * self.tether.l + self.substrate.z + self.substrate.l * 0.5
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

        # defects
        hclength = self.lipid1.l + self.methyl1.l + self.methyl2.l + self.lipid2.l
        hglength = self.headgroup1.l + self.headgroup2.l

        if self.radius_defect < (0.5 * (hclength + hglength)):
            self.radius_defect = 0.5 * (hclength + hglength)

        volhalftorus = 3.14159265359 * 3.14159265359 * (
                    self.radius_defect - (2 * hclength / 3 / 3.14159265359)) * hclength * hclength / 4
        volcylinder = 3.14159265359 * self.radius_defect * self.radius_defect * hclength
        defectarea = volhalftorus / volcylinder * (1 - self.vf_bilayer) * self.normarea

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


class tBLM_arbitrary(BLM_arbitrary):
    def __init__(self, tether=HC18, filler=bme, **kwargs):

        # add ssBLM-specific subgroups
        self.substrate = Box2Err(name='substrate')
        self.nf_tether = 0.6
        self.l_tether=3.0
        self.mult_tether = 5.0

        # change these to Headgroup2Box
        self.bME = Component2Box(name='bME', molecule=Component(name='bME', formula=filler.formula, cell_volume=filler.cell_volume, length=5.2))
        self.tether = Component2Box(name='tether', molecule=Component(name='tether', formula=tether.tether.formula, cell_volume=tether.tether.cell_volume, length=self.l_tether))
        self.tetherg = Component2Box(name='tetherg', molecule=Component(name='tetherg', formula=tether.tetherg.formula, cell_volume=tether.tetherg.cell_volume, length=10.0))

        self.substrate.l = 40
        self.substrate.z = 0
        self.substrate.nf = 1
        self.rho_substrate = 2.07e-6

        self.global_rough = 2.0

        self.vol_acyl_tether = tether.tails.cell_volume
        self.nSL_acyl_tether = tether.tails.sld * tether.tails.cell_volume * 1e-6
        self.vol_methyl_tether = tether.methyls.cell_volume
        self.nSL_methyl_tether = tether.methyls.sld * tether.methyls.cell_volume * 1e-6

        super().__init__(**kwargs)

    def fnAdjustParameters(self):
        
        self._adjust_lipids()

        # substrate
        self.substrate.vol = self.normarea * self.substrate.l
        self.substrate.nSL = self.rho_substrate * self.substrate.vol

        # set all lengths
        self.bME.z = 0.5 * self.bME.l + self.substrate.z + self.substrate.l * 0.5
        self.tether.z = 0.5 * self.tether.l + self.substrate.z + self.substrate.l * 0.5
        self.tetherg.z = self.tether.z + 0.5 * self.tether.l + 0.5 * self.tetherg.l
        
        self._adjust_z(self.tetherg.z + 0.5 * self.tetherg.l)
        
        self._adjust_defects()

        self.fnSetSigma(self.sigma)

    def _adjust_lipids(self):

        # sufficiently different from parent class that code is repeated here with small adjustments
        # TODO: split into inner and outer leaflets so this can be made more efficient
        
        self.l_lipid1 = max(self.l_lipid1, 0.01)
        self.l_lipid2 = max(self.l_lipid2, 0.01)
        for nf in self.lipid_nf:
            nf = max(nf, 0)

        # this is changed normalization behavior, but it's a bit more intuitive
        self.lipid_nf /= numpy.sum(self.lipid_nf)

        self.vf_bilayer = max(self.vf_bilayer, 1E-5)
        self.vf_bilayer = min(self.vf_bilayer, 1)

        # calculate average headgroup lengths, ignore zero volume (e.g. cholesterol)
        self.av_hg1_l = numpy.sum(numpy.array([hg.l for hg, use in zip(self.headgroups1, ~self.null_hg1) if use]) * self.lipid_nf[~self.null_hg1]) / numpy.sum(self.lipid_nf[~self.null_hg1])
        self.av_hg2_l = numpy.sum(numpy.array([hg.l for hg, use in zip(self.headgroups2, ~self.null_hg2) if use]) * self.lipid_nf[~self.null_hg2]) / numpy.sum(self.lipid_nf[~self.null_hg2])
        
        # outer hydrocarbons
        l_ohc = self.l_lipid2
        nf_ohc_lipid = self.lipid_nf
        V_ohc = numpy.sum(nf_ohc_lipid * (self.vol_acyl_lipids - self.vol_methyl_lipids))
        nSL_ohc = numpy.sum(nf_ohc_lipid * (self.nsl_acyl_lipids - self.nsl_methyl_lipids))
        #V_ohc = nf_ohc_lipid * (self.volacyllipid - self.volmethyllipid) + nf_ohc_lipid_2 * (self.volacyllipid_2 - self.volmethyllipid_2) + nf_ohc_lipid_3 * (self.volacyllipid_3 - self.volmethyllipid_3) + nf_ohc_chol * self.volchol
        #nSL_ohc = nf_ohc_lipid * (self.nslacyllipid - self.nslmethyllipid) + nf_ohc_lipid_2 * (self.nslacyllipid_2 - self.nslmethyllipid_2) + nf_ohc_lipid_3 * (self.nslacyllipid_3 - self.nslmethyllipid_3) + nf_ohc_chol * self.nslchol
        
        self.normarea = V_ohc / l_ohc
        c_s_ohc = self.vf_bilayer
        c_A_ohc = 1
        c_V_ohc = 1
        
        self.lipid2.l = l_ohc
        self.lipid2.vol = V_ohc
        self.lipid2.nSL = nSL_ohc
        self.lipid2.nf = c_s_ohc * c_A_ohc * c_V_ohc
        
        # outer methyl
        #nf_om_lipid = nf_ohc_lipid
        #nf_om_lipid_2 = nf_ohc_lipid_2
        #nf_om_lipid_3 = nf_ohc_lipid_3
        nf_om_lipid = nf_ohc_lipid
        #V_om = nf_om_lipid * self.volmethyllipid + nf_om_lipid_2 * self.volmethyllipid_2 + nf_om_lipid_3 * self.volmethyllipid_3
        V_om = numpy.sum(nf_om_lipid * self.vol_methyl_lipids)   # cholesterol has a methyl volume of zero and does not contribute
        l_om = l_ohc * V_om / V_ohc
        #nSL_om = nf_om_lipid * self.nslmethyllipid + nf_om_lipid_2 * self.nslmethyllipid_2 + nf_om_lipid_3 * self.nslmethyllipid_3
        nSL_om = numpy.sum(nf_om_lipid * self.nsl_methyl_lipids)
        
        c_s_om = c_s_ohc
        c_A_om = 1
        c_V_om = 1
        
        self.methyl2.l = l_om
        self.methyl2.vol = V_om
        self.methyl2.nSL = nSL_om
        self.methyl2.nf = c_s_om * c_A_om * c_V_om
        
        # inner hydrocarbons
        l_ihc = self.l_lipid1
        
        nf_ihc_lipid = nf_ohc_lipid * (1 - self.nf_tether)
        nf_ihc_tether = self.nf_tether
        V_ihc = numpy.sum(nf_ihc_lipid * (self.vol_acyl_lipids - self.vol_methyl_lipids)) + nf_ihc_tether * (self.vol_acyl_tether - self.vol_methyl_tether)
        nSL_ihc = numpy.sum(nf_ihc_lipid * (self.nsl_acyl_lipids - self.nsl_methyl_lipids)) + nf_ihc_tether * (self.nSL_acyl_tether - self.nSL_methyl_tether)
        
        c_s_ihc = self.vf_bilayer
        c_A_ihc = self.normarea * l_ihc / V_ihc
        c_V_ihc = 1
        
        self.lipid1.l = l_ihc
        self.lipid1.vol = V_ihc
        self.lipid1.nSL = nSL_ihc
        self.lipid1.nf = c_s_ihc * c_A_ihc * c_V_ihc
        
        # inner methyl
        nf_im_lipid = nf_ihc_lipid
        nf_im_tether = nf_ihc_tether
        V_im = numpy.sum(nf_im_lipid * self.vol_methyl_lipids) + nf_im_tether * self.vol_methyl_tether
        l_im = l_ihc * V_im / V_ihc
        nSL_im = numpy.sum(nf_im_lipid * self.nsl_methyl_lipids) + nf_im_tether * self.nSL_methyl_tether
        
        c_s_im = c_s_ihc
        c_A_im = c_A_ihc
        c_V_im = 1
        
        self.methyl1.l = l_im
        self.methyl1.vol = V_im
        self.methyl1.nSL = nSL_im
        self.methyl1.nf = c_s_im * c_A_im * c_V_im

        # headgroup size and position
        for hg1, nf_ihc, hg2, nf_ohc in zip(self.headgroups1, nf_ihc_lipid, self.headgroups2, nf_ohc_lipid):
            hg2.nf = c_s_ohc * c_A_ohc * nf_ohc * (1 - self.hc_substitution_2)
            hg1.nf = c_s_ihc * c_A_ihc * nf_ihc * (1 - self.hc_substitution_1)
            if hasattr(hg1, 'fnAdjustParameters'):
                hg1.fnAdjustParameters()
            if hasattr(hg2, 'fnAdjustParameters'):
                hg2.fnAdjustParameters()

        # tether glycerol part -- transliteration here
        V_tg = self.tetherg.vol

        c_s_tg = c_s_ihc
        c_A_tg = c_A_ihc
        c_V_tg = nf_ihc_tether * (1 - self.hc_substitution_2)

        self.tetherg.l = self.tetherg.vol / ((self.vol_acyl_tether - self.vol_methyl_tether) / self.lipid1.l) / 0.9
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
        #self.headgroup1.l = 9.575
        V_bME = self.bME.vol

        d1 = self.av_hg1_l + self.bME.l - self.tether.l - self.tetherg.l
        if (d1 > 0):
            self.bME.l = self.bME.l - d1 / 2
            self.av_hg1_l = self.av_hg1_l - d1 / 2

        if ((
                self.tether.nf * self.tether.vol / self.tether.l + self.mult_tether * self.tether.nf * self.bME.vol / self.bME.l) > self.normarea):
            # print(self.tether.nf, self.tether.vol, self.tether.l, self.mult_tether, self.bME.vol, self.bME.l)
            self.mult_tether = ((self.normarea - self.tether.nf * self.tether.vol / self.tether.l) / (
                        self.bME.vol / self.bME.l)) / self.tether.nf
            if (self.mult_tether < 0):
                self.mult_tether = 0

        self.bME.nf = self.tether.nf * self.mult_tether  # 2.333

    def fnSetSigma(self, sigma):
        super().fnSetSigma(sigma)
        self.substrate.fnSetSigma(self.global_rough)
        self.bME.fnSetSigma(self.global_rough)
        self.tether.fnSetSigma(self.global_rough, sigma)
        self.tetherg.fnSetSigma(sigma)

    def fnGetLowerLimit(self):
        return self.substrate.fnGetLowerLimit()
        # does this make sense since this goes to negative z and isn't intended to be used?

    def fnSet(self, _sigma, _bulknsld, _global_rough, _rho_substrate, _nf_tether, _mult_tether, _l_tether, _l_lipid1,
              _l_lipid2, _vf_bilayer, _nf_lipids=None, _hc_substitution_1=0.,
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
        if _nf_lipids is not None:   # pass None to keep lipid_nf unchanged
            assert len(_nf_lipids)==len(self.lipids), 'nf_lipids must be same length as number of lipids in bilayer, not %i and %i' % (len(_nf_lipids), len(self.lipids))
            self.lipid_nf = numpy.array(_nf_lipids)
        self.hc_substitution_1 = _hc_substitution_1
        self.hc_substitution_2 = _hc_substitution_2
        self.radius_defect = _radius_defect

        self.fnAdjustParameters()

# ------------------------------------------------------------------------------------------------------
# Hermite Spline
# ------------------------------------------------------------------------------------------------------

"""
Notes on usage:
1. Instantiate the spline object, e.g. h = SLDHermite()
    NB: the only initialization variable that isn't overwritten by fnSetRelative is dnormarea, so the others have been removed
        and dnSLD has been moved to fnSetRelative instead of the initialization. This is to keep the fnSetRelative function calls
        similar between Hermite and SLDHermite (except in Hermite nSLD is a constant and in SLDHermite it's a list of control points)
2. Set all internal parameters, e.g. damping parameters, monotonic spline, etc.
3. Call fnSetRelative.
    NB: for speed, the spline interpolators are stored in the object. Only fnSetRelative will update them!
"""

def catmull_rom(xctrl, yctrl, **kwargs):
    """returns a catmull_rom spline using ctrl points xctrl, yctrl"""
    assert(len(xctrl)==len(yctrl)) # same shape

    x_xtnd = numpy.insert(xctrl, 0, xctrl[0] - (xctrl[1]-xctrl[0]))
    x_xtnd = numpy.append(x_xtnd, xctrl[-1]-xctrl[-2] + xctrl[-1])

    y_xtnd = numpy.insert(yctrl, 0, yctrl[0])
    y_xtnd = numpy.append(y_xtnd, yctrl[-1])
    dydx = 0.5 * (y_xtnd[2:] - y_xtnd[:-2]) / (x_xtnd[2:] - x_xtnd[:-2])

    return CubicHermiteSpline(xctrl, yctrl, dydx, **kwargs)


class Hermite(nSLDObj):
    def __init__(self, dnormarea, **kwargs):
        super().__init__(**kwargs)
        self.numberofcontrolpoints = 10
        self.nSLD=None
        self.normarea=dnormarea
        self.monotonic=True
        self.damping=True
        self.dampthreshold=0.001
        self.dampFWHM=0.0002
        self.damptrigger=0.04
        self.dstartposition = 0.0
        
        self.dp     = numpy.arange(self.numberofcontrolpoints)
        self.vf     = numpy.zeros(self.numberofcontrolpoints)
        self.damp   = numpy.zeros(self.numberofcontrolpoints)

        # TODO: get the interface with a previous layer correct by defining an erf function at the first control point or between the first two control points.

    def _apply_damping(self):
        """ only called from _set_area_spline, which is called from fnSetRelative """
        dampfactor = numpy.ones_like(self.vf)
        if self.damping:
            above_damptrigger = numpy.where(self.vf > self.damptrigger)[0] # find index of first point beyond damping trigger
            if len(above_damptrigger) > 0 : # if no points above trigger, damping doesn't apply
                dampfactor[above_damptrigger[0] + 1:] = 1./(1+numpy.exp(-2.1*(self.vf[above_damptrigger[0] + 1:]-self.dampthreshold)/self.dampFWHM))    # does nothing if peaked point is at the end
                dampfactor = numpy.cumprod(dampfactor)
        
        self.damp = self.vf * dampfactor

    def _set_area_spline(self):

        self._apply_damping()

        if self.monotonic:  # monotone interpolation
            self.area_spline = PchipInterpolator(self.dp, self.damp, extrapolate=False)

        else:   # catmull-rom

            self.area_spline = catmull_rom(self.dp, self.damp, extrapolate=False)
        
        self.area_spline_integral = self.area_spline.antiderivative()

    def fnGetProfiles(self, z):

        vf = self.area_spline(z)
        vf[numpy.isnan(vf)] = 0.0
        vf[vf < 0] = 0.0
        area = vf * self.normarea * self.nf

        sld = self.fnGetnSLDProfile(z)

        return area, area * numpy.gradient(z) * sld, sld

    def fnGetnSLDProfile(self, z):

        return self.nSLD * numpy.ones_like(z)

    def fnGetnSLD(self, dz): 
        return self.nSLD

    def fnGetLowerLimit(self): 
        return self.dp[0]

    def fnGetUpperLimit(self): 
        return self.dp[-1]

    def fnGetVolume(self, z1, z2):
        """ use stored antiderivatives to calculate volume """
        # make sure z1 and z2 are in the defined interval. If both are above or below, result will be zero
        z1 = max(self.fnGetLowerLimit(), z1)
        z1 = min(self.fnGetUpperLimit(), z1)
        z2 = max(self.fnGetLowerLimit(), z2)
        z2 = min(self.fnGetUpperLimit(), z2)

        return (self.area_spline_integral(z2) - self.area_spline_integral(z1)) * self.nf * self.normarea

    def fnSetNormarea(self, dnormarea): 
        self.normarea = dnormarea

    def fnSetnSLD(self, dnSLD): 
        self.nSLD = dnSLD

    def fnSetRelative(self, dSpacing, dStart, dDp, dVf, dnSLD, dnf):
        self.fnSetnSLD(dnSLD)
        self.vf = numpy.array(dVf)
        self.numberofcontrolpoints = len(self.vf)
        self.dp = dStart + dSpacing*numpy.arange(self.numberofcontrolpoints) + numpy.array(dDp)
        # make sure the control points have compatible shapes (previous command should fail if not)
        assert(self.vf.shape == self.dp.shape)
        self.dstartposition = dStart
        self.nf = dnf
        self._set_area_spline()

    def fnWritePar2File(self, fp, cName, z):
        fp.write("Hermite "+cName+" numberofcontrolpoints "+str(self.numberofcontrolpoints)+" normarea "
                 +str(self.normarea)+" nf "+str(self.nf)+" \n")
        self.fnWriteData2File(fp, cName, z)

class SLDHermite(Hermite):
    def __init__(self, dnormarea, **kwargs):
        super().__init__(dnormarea, **kwargs)
        self.sld   = numpy.zeros(self.numberofcontrolpoints)

        # TODO: get the interface with a previous layer correct by defining an erf function at the first control point or between the first two control points.

    def _set_sld_spline(self):

        if self.monotonic:  # monotone interpolation
            self.sld_spline = PchipInterpolator(self.dp, self.sld, extrapolate=False)

        else:   # catmull-rom
            self.sld_spline = catmull_rom(self.dp, self.sld, extrapolate=False)

    def fnGetnSLDProfile(self, z):

        sld = self.sld_spline(z)
        # deal with out of range values
        sld[z < self.dp[0]] = self.sld[0]
        sld[z > self.dp[-1]] = self.sld[-1]
        assert(not numpy.any(numpy.isnan(sld))) # shouldn't be any NaNs left

        return sld

    def fnSetnSLD(self, dnSLD): 
        self.nSLD = None

    def fnSetRelative(self, dSpacing, dStart, dDp, dVf, dnSLD, dnf):
        self.vf = numpy.array(dVf)
        self.numberofcontrolpoints = len(self.vf)
        self.dp = dStart + dSpacing*numpy.arange(self.numberofcontrolpoints) + numpy.array(dDp)
        # make sure the control points have compatible shapes (previous command should fail if not)
        assert(self.vf.shape == self.dp.shape)
        self.sld = numpy.array(dnSLD)
        assert(self.sld.shape == self.dp.shape)
        self.dstartposition = dStart
        self.nf = dnf
        self._set_area_spline()
        self._set_sld_spline()

    def fnWritePar2File(self, fp, cName, z):
        fp.write("Hermite "+cName+" numberofcontrolpoints "+str(self.numberofcontrolpoints)+" normarea "
                 +str(self.normarea)+" nf "+str(self.nf)+" \n")
        self.fnWriteData2File(fp, cName, z)

class ContinuousEuler(nSLDObj):
    """ Uses scipy.spatial library to do Euler rotations in real time"""
    def __init__(self, fn8col, rotcenter=None, **kwargs):
        """ requires file name fn8col containing 8-column data, any header information commented with #:
            1. residue number
            2. x coordinate
            3. y coordinate
            4. z coordinate
            5. residue volume
            6. electon scattering length
            7. neutron scattering length (H)
            8. neutron scattering length (D)

            Use helper function pdbto8col to create this data file

            The simplest approach is to provide coordinates relative to the center of mass; then,
            "z" corresponds to the absolute z position of the center of mass of the object. Otherwise,
            if the rotation center "rotcenter" [x0, y0, z0] is specified, all coordinates are translated
            by "rotcenter" and "z" is the absolute z position of the rotation center.
        """
        super().__init__(**kwargs)
        self.fn = fn8col
        resdata = numpy.loadtxt(fn8col)
        self.resnums = resdata[:,0]
        self.rescoords = resdata[:,1:4]
        if rotcenter is not None:
            rotcenter = numpy.array(rotcenter)
            assert(rotcenter.shape == self.rescoords[0,:].shape)
            self.rescoords -= rotcenter
        self.rotcoords = numpy.zeros_like(self.rescoords)
        self.resscatter = resdata[:, 4:]
        self.alpha = 0.
        self.beta = 0.
        self.sigma = 2.
        self.z = 0.
        self.nf = 1.
        self.protexchratio = 1.
        self.R = Rotation.from_euler('zy', [self.alpha, self.beta], degrees=True)

        # TODO: Would it make sense to have a concept of "normarea"? Then there could be a "volume fraction" concept so that
        # max(area) = volume_fraction * normarea

    def _apply_transform(self):
        
        self.R = Rotation.from_euler('zy', [self.alpha, self.beta], degrees=True)
        self.rotcoords = self.R.apply(self.rescoords)
        self.rotcoords[:,2] += self.z

    def fnGetProfiles(self, z, bulknsld=None, xray=False):
        # xray=True for xray profile
        # xray=False (default) for a neutron probe; bulknsld determines fraction of exchangeable hydrogens used

        # get area profile
        dz = z[1] - z[0]                            # calculate step size (MUST be uniform)
        zbins = numpy.append(z, z[-1]+dz) - 0.5 * dz    # bin edges; z is treated as bin centers here
        # TODO: Smoothing with gaussian_filter requires constant step size. What happens if z is a single point?
        h = numpy.histogram(self.rotcoords[:,2], bins=zbins, weights=self.resscatter[:,0])
        area = gaussian_filter(h[0], self.sigma / dz, order=0, mode='constant', cval=0)
        area /= dz

        if xray:
            h = numpy.histogram(self.rotcoords[:,2], bins=zbins, weights=self.resscatter[:,1])
            nsl = gaussian_filter(h[0], self.sigma / dz, order=0, mode='constant', cval=0)
        else:
            # get nslH profile
            h = numpy.histogram(self.rotcoords[:,2], bins=zbins, weights=self.resscatter[:,2])
            nslH = gaussian_filter(h[0], self.sigma / dz, order=0, mode='constant', cval=0)

            if bulknsld is not None:
                # get nslD profile
                h = numpy.histogram(self.rotcoords[:,2], bins=zbins, weights=self.resscatter[:,3])
                nslD = gaussian_filter(h[0], self.sigma / dz, order=0, mode='constant', cval=0)
                fracD = self.protexchratio * (bulknsld + 0.56e-6) / (6.36e-6 + 0.56e-6)
                nsl = fracD * nslD + (1 - fracD) * nslH

            else:
                nsl = nslH

        nsld = numpy.zeros_like(z)
        pos = (area > 0)
        nsld[pos] = nsl[pos] / (area[pos] * numpy.gradient(z)[pos])

        return area * self.nf, nsl * self.nf, nsld

    def fnGetVolume(self, z1, z2):
        """ Calculates volume based on the number of residues of the rotated molecule located between
            z positions z1 and z2 (inclusive).
            
            Note: the result is (slightly) different from integrating the area, because the roughness has already
            been applied to the area. However, this remains more accurate than the roughened value because
            the integration limits  will typically correspond to the limits of a lipid Box2Err function
            which are themselves defined before the roughness is applied."""
        # use a single bin defined by bin edges z1 and z2.
        # First [0] selects the histogram array, second [0] gets the single value from the array

        volume = numpy.histogram(self.rotcoords[:,2], bins=[z1,z2], weights=self.resscatter[:,0])[0][0]

        return volume * self.nf

    def fnSet(self, alpha, beta, zpos, sigma, nf):

        self.alpha = alpha
        self.beta = beta
        self.z = zpos
        self.nf = nf
        self.sigma = sigma
        self._apply_transform()

    def fnWritePar2File(self, fp, cName, z):
        fp.write("ContinuousEuler "+cName+" StartPosition "+str(self.z)+" Alpha "
                 +str(self.alpha)+" Beta "+str(self.beta) +" nf "+str(self.nf)+" \n")
        self.fnWriteData2File(fp, cName, z)

def pdbto8col(pdbfilename, datfilename, selection='all', center_of_mass=numpy.array([0,0,0]), deuterated_residues=None):
    """ Creates an 8-column data file for use with ContinuousEuler from a pdb file\
        with optional selection. "center_of_mass" is the position in space at which to position the
        molecule's center of mass. "deuterated_residues" is a list of residue IDs for which to use deuterated values"""
    
    # not currently used but useful
    rshort = dict({'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU', 'Q':'GLN', 'G':'GLY', 'H':'HIS',
               'I': 'ILE', 'L':'LEU', 'K':'LYS', 'M':'MET', 'F': 'PHE', 'P': 'PRO', 'S':'SER', 'T':'THR', 'Y':'TYR',
                'V': 'VAL', 'W': 'TRP'})

    # residue name : [vol Ang^3, eSL, SLprot Ang, SLdeut Ang, #exchngH]     
    # volumes from Chothia, C. (1975) Nature 254, 304-308.
    residprop = {
            'ALA' : [91.5, 38, 20.1466,  61.7942, 1],
            'ARG' : [202.1, 85, 56.9491, 129.8324, 6],
            'ASP' : [124.5, 59, 42.149,  73.3816, 1],
            'ASN' : [135.2, 60, 45.7009, 76.9366, 3],
            'CYS' : [105.6, 54, 26.7345, 57.9702, 2],
            'GLU' : [155.1, 67, 41.371,  93.372, 1],
            'GLN' : [161.1, 68, 44.8675, 96.927, 3],
            'GLY' : [66.4,  30, 20.98,   41.8038, 1],
            'HSD' : [167.3, 72, 55.0709, 107.1304, 2],
            'HIS' : [167.3, 72, 55.0709, 107.1304, 2],
            'HSE' : [167.3, 72, 55.0709, 107.1304, 2],
            'HSP' : [167.3, 72, 55.0709, 107.1304, 3],
            'ILE' : [168.8, 62, 17.6464, 121.7654, 1],
            'LEU' : [167.9, 62, 17.6464, 121.7654, 1],
            'LYS' : [171.3, 71, 30.7473, 124.4544, 4],
            'MET' : [170.8, 70, 21.3268, 104.622, 1],
            'PHE' : [203.4, 78, 45.0734, 128.3686, 1],
            'PRO' : [129.3, 52, 22.2207, 95.104, 0],
            'MLY' : [224.58, 96, 16.2, 16.2, 0],
            'MSE' : [253.56, 98, 21.0, 21.0, 0]
            } #deuteration not correct for MLY and MSE
    residprop['SER'] = [99.1,  46, 29.6925, 60.9282, 2]
    residprop['DAV'] = [99.1,  46, 290000.6925, 600000.9282, 2] # huge scattering lengths for testing!
    residprop['THR'] = [122.1, 54, 25.1182, 87.5896, 2]
    residprop['TRP'] = [237.6, 98, 67.7302, 151.0254, 2]
    residprop['TYR'] = [203.6, 86, 54.6193, 127.5026, 2]
    residprop['VAL'] = [141.7, 54, 18.4798, 101.775, 1]
    #volume of Zn2+ from Obst et al. J.Mol.Model 1997,3,224-232 (Zn-O g(r))
    residprop['ZN2'] = [100.5, 28, 5.68, 5.68, 0] 

    import MDAnalysis
    molec = MDAnalysis.Universe(pdbfilename)
    sel = molec.select_atoms(selection)
    Nres = sel.n_residues
    
    if not Nres:
        print('Warning: no atoms selected')

    sel.translate(-sel.center_of_mass() + center_of_mass)
    
    resnums = []
    rescoords = []
    resscatter = []
    deut_header = ''
    HSL = -3.7409
    DSL = 6.671
    for i in range(Nres):
        resnums.append(molec.residues[i].resid)
        rescoords.append(molec.residues[i].atoms.center_of_mass())
        resscatter.append(residprop[molec.residues[i].resname])

    resnums = numpy.array(resnums)
    rescoords = numpy.array(rescoords)
    resscatter = numpy.array(resscatter)

    # replace base value in nsl calculation with proper deuterated scattering length\
    resnsl = resscatter[:,2]
    if deuterated_residues is not None:
        deuterated_indices = numpy.searchsorted(resnums, deuterated_residues)
        resnsl = resscatter[:,2]
        resnsl[deuterated_indices] = resscatter[deuterated_indices, 3]
        deut_header = 'deuterated residues: ' + ', '.join(map(str, deuterated_residues)) + '\n'


    resesl = resscatter[:,1] * 2.8e-5
    resnslH = (resnsl + HSL * resscatter[:,4]) * 1e-5
    resnslD = (resnsl + DSL * resscatter[:,4]) * 1e-5

    numpy.savetxt(datfilename, numpy.hstack((resnums[:,None], rescoords, resscatter[:,0][:,None], resesl[:,None], resnslH[:, None], resnslD[:,None])), delimiter='\t',
            header=pdbfilename + '\n' + deut_header + 'resid\tx\ty\tz\tvol\tesl\tnslH\tnslD')

    return datfilename      # allows this to be fed into ContinuousEuler directly

class DiscreteEuler(nSLDObj):
    """ Uses precalculated Euler rotation densities"""
    def __init__(self, generic_filename, betas, gammas, betafmt='%i', gammafmt='%i', **kwargs):
        """ requires precalculated 4-column files for each angle beta and gamma. These should not have any smoothing applied.
            Each file must have same z vector and one header row. z can be negative.
            generic_filename contains the path and a generic file name with positions for angles beta and gamma are marked with <beta> and <gamma>.
            Format strings other than integer angles can be specified with betafmt and gammafmt keyword arguments.
            "betas" and "gammas" are lists or numpy arrays containing the beta and gamma points to be loaded.

        Example: DiscreteEuler('./dat/exou_beta<beta>_gamma<gamma>.txt', range(0, 190, 10), range(0, 370, 10))
        """
        super().__init__(**kwargs)
        betas = numpy.array(betas, dtype=float)
        gammas = numpy.array(gammas, dtype=float)
        areadata = None
        for i, beta in enumerate(betas):
            for j, gamma in enumerate(gammas):
                fn = generic_filename.replace('<beta>', betafmt % beta).replace('<gamma>', gammafmt % gamma)
                d = numpy.loadtxt(fn, skiprows=1)
                if areadata is None:
                    zdata = d[:,0]
                    areadata = numpy.zeros((len(betas), len(gammas), len(zdata)))
                    nslHdata = numpy.zeros_like(areadata)
                    nslDdata = numpy.zeros_like(areadata)
                areadata[i, j, :] = d[:,1]
                nslHdata[i, j, :] = d[:,2]
                nslDdata[i, j, :] = d[:,3]

        self.betas = betas
        self.gammas = gammas
        self.areadata = areadata
        self.zdata = zdata
        self.nslHdata = nslHdata
        self.nslDdata = nslDdata
        self.beta = 0.         # current beta rotation
        self.gamma = 0.          # current gamma rotation
        self.sigma = 2.         # smoothing function
        self.z = 0.             # z offset value
        self.nf = 1.            # number fraction
        self.protexchratio = 1. # proton exchange ratio

        # TODO: Would it make sense to have a concept of "normarea"? Then there could be a "volume fraction" concept so that
        # max(area) = volume_fraction * normarea

    def fnGetProfiles(self, z, bulknsld=None):

        # perform interpolation
        self._set_interpolation_points(z)
        area = interpn((self.betas, self.gammas, self.zdata + self.z), self.areadata, self._interppoint, bounds_error=False, fill_value=0.0)
        nslH = interpn((self.betas, self.gammas, self.zdata + self.z), self.nslHdata, self._interppoint, bounds_error=False, fill_value=0.0)
        nslD = interpn((self.betas, self.gammas, self.zdata + self.z), self.nslDdata, self._interppoint, bounds_error=False, fill_value=0.0)

        # apply roughness (sigma)
        dz = z[1] - z[0]
        area = gaussian_filter(area, self.sigma / dz, order=0, mode='constant', cval=0)
        #area /= dz

        nslH = gaussian_filter(nslH, self.sigma / dz, order=0, mode='constant', cval=0)

        if bulknsld is not None:
            # get nslD profile
            nslD = gaussian_filter(nslD, self.sigma / dz, order=0, mode='constant', cval=0)
            fracD = self.protexchratio * (bulknsld + 0.56e-6) / (6.36e-6 + 0.56e-6)
            nsl = fracD * nslD + (1 - fracD) * nslH

        else:
            nsl = nslH

        nsld = numpy.zeros_like(z)
        pos = (area > 0)
        nsld[pos] = nsl[pos] / (area[pos] * numpy.gradient(z)[pos])

        return area * self.nf, nsl * self.nf, nsld

    def _set_interpolation_points(self, z):
        self._interppoint = numpy.zeros((len(z), 3))
        self._interppoint[:,:-1] = [self.beta, self.gamma]
        self._interppoint[:, -1] = z

    def fnGetVolume(self, z1, z2):
        """ Calculates volume based on the number of residues of the rotated molecule located between
            z positions z1 and z2 (inclusive).
            
            Note: the result is (slightly) different from integrating the area, because the roughness has already
            been applied to the area. However, this remains more accurate than the roughened value because
            the integration limits  will typically correspond to the limits of a lipid Box2Err function
            which are themselves defined before the roughness is applied."""

        # For volume calculation, use intrinsic z vector, and no smoothing.
        zvol = self.zdata + self.z
        self._set_interpolation_points(zvol)

        area = interpn((self.betas, self.gammas, zvol), self.areadata, self._interppoint, bounds_error=False, fill_value=0.0)

        crit = (zvol > min(z1, z2)) & (zvol < max(z1, z2))
        if numpy.any(crit):
            volume= numpy.trapz(area[crit], zvol[crit])
        else:
            volume = 0.0


        return volume * self.nf

    def fnSet(self, beta, gamma, zpos, sigma, nf):

        self.beta = beta
        self.gamma = gamma
        self.z = zpos
        self.nf = nf
        self.sigma = sigma

    def fnWritePar2File(self, fp, cName, z):
        fp.write("DiscreteEuler "+cName+" StartPosition "+str(self.z)+" beta "
                 +str(self.beta)+" gamma "+str(self.gamma) +" nf "+str(self.nf)+" \n")
        self.fnWriteData2File(fp, cName, z)