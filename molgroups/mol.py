import numpy
from scipy.special import erf
from scipy.integrate import trapezoid
from scipy.interpolate import PchipInterpolator, CubicHermiteSpline, interpn
from scipy.spatial.transform import Rotation
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.interpolation import shift
from scipy.signal import peak_widths
import sys
import warnings

from refl1d.polymer import mushroom_math, smear

from periodictable.fasta import Molecule, AMINO_ACID_CODES as aa
from periodictable.core import default_table
from periodictable.fasta import xray_sld
from periodictable.fasta import D2O_SLD, H2O_SLD

from molgroups import components as cmp

D2O_SLD *= 1e-6
H2O_SLD *= 1e-6


def pdbto8col(pdbfilename, datfilename, selection='all', center_of_mass=numpy.array([0, 0, 0]),
              deuterated_residues=None, xray_wavelength=1.5418):
    """
    Creates an 8-column data file for use with ContinuousEuler from a pdb file\
    with optional selection. "center_of_mass" is the position in space at which to position the
    molecule's center of mass. "deuterated_residues" is a list of residue IDs for which to use deuterated values
    """
    import MDAnalysis
    from MDAnalysis.lib.util import convert_aa_code

    if deuterated_residues is None:
        deuterated_residues = []

    elements = default_table()

    molec = MDAnalysis.Universe(pdbfilename)
    sel = molec.select_atoms(selection)
    Nres = sel.n_residues

    if not Nres:
        print('Warning: no atoms selected')

    sel.translate(-sel.center_of_mass() + center_of_mass)

    resnums = []
    rescoords = []
    resscatter = []
    resvol = numpy.zeros(Nres)
    resesl = numpy.zeros(Nres)
    resnslH = numpy.zeros(Nres)
    resnslD = numpy.zeros(Nres)
    deut_header = ''

    for i in range(Nres):
        resnum = molec.residues[i].resid
        resnums.append(resnum)
        rescoords.append(molec.residues[i].atoms.center_of_mass())
        key = convert_aa_code(molec.residues[i].resname)
        if resnum in deuterated_residues:
            resmol = Molecule(name='Dres', formula=aa[key].formula.replace(elements.H, elements.D),
                              cell_volume=aa[key].cell_volume)
        else:
            resmol = aa[key]
        resvol[i] = resmol.cell_volume
        # TODO: Make new column for xray imaginary part (this is real part only)
        resesl[i] = resmol.cell_volume * xray_sld(resmol.formula, wavelength=xray_wavelength)[0] * 1e-6
        resnslH[i] = resmol.cell_volume * resmol.sld * 1e-6
        resnslD[i] = resmol.cell_volume * resmol.Dsld * 1e-6

    resnums = numpy.array(resnums)
    rescoords = numpy.array(rescoords)

    average_sldH = numpy.sum(resnslH[:, ]) / numpy.sum(resvol[:, ])
    average_sldD = numpy.sum(resnslD[:, ]) / numpy.sum(resvol[:, ])
    average_header = f'Average nSLD in H2O: {average_sldH}\nAverage nSLD in D2O: {average_sldD}\n'

    # replace base value in nsl calculation with proper deuterated scattering length
    # resnsl = resscatter[:, 2]
    if deuterated_residues is not None:
        deut_header = 'deuterated residues: ' + ', '.join(map(str, deuterated_residues)) + '\n'

    numpy.savetxt(datfilename, numpy.hstack((resnums[:, None], rescoords, resvol[:, None], resesl[:, None],
                                             resnslH[:, None], resnslD[:, None])), delimiter='\t',
                  header=pdbfilename + '\n' + deut_header + average_header + 'resid\tx\ty\tz\tvol\tesl\tnslH\tnslD')

    return datfilename  # allows this to be fed into ContinuousEuler directly


class nSLDObj:
    def __init__(self, name=None):
        self.bWrapping = True
        self.bConvolution = False
        self.bProtonExchange = False
        self.dSigmaConvolution = 1
        self.iNumberOfConvPoints = 7
        self.absorb = 0
        self.nf = 0
        self.sigma = 0
        self.bulknsld = None
        self.z = 0.

        # allows to flip groups, eg. inner vs. outer leaflet, outer leaflet is default
        self.flip = False

        # stored profile for derived parameter calculation without recalculation of the profile
        self.area = None
        self.sl = None
        self.sld = None
        self.zaxis = None

        if name is not None:
            self.name = name

    def fnGetAbsorb(self, z):
        return self.absorb

    def fnGetArea(self, z1=None, z2=None, recalculate=True):
        if z1 is None:
            # no zaxis range given, take entire stored self.zaxis
            if recalculate or self.area is None:
                print('Scalar-defined interval for area calculation requires using a stored profile')
                sys.exit(1)
            else:
                return self.area
        elif z2 is None:
            # z1 is a numpy array over the zaxis range to be evaluated
            if recalculate or self.area is None:
                # new evaluation instead of using any stored area profiles
                return self.fnGetProfiles(z1)[0]
            else:
                # get area from previously computed profile
                return self.area[numpy.where((z1[0] <= self.zaxis) & (self.zaxis <= z1[-1]))]
        else:
            # z1 and z2 are scalars and representing the interval over which the area is evaluated
            if recalculate or self.area is None:
                print('Scalar-defined interval for area calculation requires using a stored profile')
                sys.exit(1)
            else:
                return self.area[numpy.where((z1 <= self.zaxis) & (self.zaxis <= z2))]

    def fnGetConvolutedArea(self, dz):
        # returns an n-point gaussian interpolation of the area within 4 sigma
        # all area calculations are routed through this function, whether they use convolution or not
        # convolution works only for objects with fixed nSLD. Broadening an nSLD profile is not as direct as
        # broadening a nSL profile. For reason, however, objects report a nSLD(z) and not a nSL(z)
        # if it becomes necessary to broaden profiles with variable nSLD, structural changes to the code
        # have to be implemented.
        # TODO: use scipy image filters or numpy convolution to do this.
        if self.bConvolution:
            dnormsum = 0
            dsum = 0
            for i in range(self.iNumberOfConvPoints):
                dd = 8 / float(self.iNumberOfConvPoints * i) - 4
                dgauss = numpy.exp((-0.5) * dd * dd)  # (sigma_convolution)^2/(sigma_convolution)^2 cancels
                dnormsum += dgauss
                dsum += self.fnGetArea(dz + dd * self.dSigmaConvolution) * dgauss
            if dnormsum != 0:
                return dsum / dnormsum
            else:
                return 0
        else:
            return self.fnGetArea(dz)

    def fnGetnSLD(self, z):
        # generic area function using fnGetProfile, can be replaced with a more efficient routine in derived objects
        return self.fnGetProfiles(z)[2][z]

    def fnGetProfiles(self, z):
        # returns (area, nsl, nsld)
        raise NotImplementedError()

    def fnGetVolume(self, z1=None, z2=None, recalculate=True):
        area = self.fnGetArea(z1=z1, z2=z2, recalculate=recalculate)
        if z1 is None:
            # no zaxis range given, take entire stored self.zaxis
            if recalculate or self.area is None:
                print('Scalar-defined interval for area calculation requires using a stored profile')
                sys.exit(1)
            else:
                zaxis_gradient = numpy.gradient(self.zaxis)
        elif z2 is None:
            # z1 is a numpy array over the zaxis range to be evaluated
            zaxis_gradient = numpy.gradient(z1)
        else:
            # z1 and z2 are scalars and representing the interval over which the area is evaluated
            if recalculate or self.area is None:
                print('Scalar-defined interval for area calculation requires using a stored profile')
                sys.exit(1)
            else:
                z_interval = self.zaxis[numpy.where((z1 <= self.zaxis) & (self.zaxis <= z2))]
                if len(z_interval) > 1:
                    zaxis_gradient = numpy.gradient(z_interval)
                else:
                    return 0.0
        return numpy.sum(area * zaxis_gradient)

    def fnGetZ(self):
        return self.z

    def fnSetBulknSLD(self, bulknsld):
        self.bulknsld = bulknsld

    def fnSetConvolution(self, sigma_convolution, iNumberOfConvPoints):
        self.bConvolution = True
        self.dSigmaConvolution = sigma_convolution
        self.iNumberOfConvPoints = iNumberOfConvPoints

    @staticmethod
    def fnWriteConstant(fp, name, darea, dSLD, z):
        header = f"Constant {name} area {darea} \nz_{name} a_{name} nsl_{name}"
        area = darea * numpy.ones_like(z)
        nsl = dSLD * darea * numpy.gradient(z)
        A = numpy.vstack((z, area, nsl)).T
        numpy.savetxt(fp, A, fmt='%0.6e', delimiter=' ', comments='', header=header)
        fp.write('\n')

    @staticmethod
    def fnWriteConstant2Dict(rdict, name, darea, dSLD, z):
        rdict[name] = {}
        rdict[name]['header'] = f"Constant {name} area {darea}"
        rdict[name]['area value'] = darea
        rdict[name]['sld value'] = dSLD
        constarea = numpy.ones_like(z) * darea
        constsld = numpy.ones_like(z) * dSLD
        rdict[name]['zaxis'] = z
        rdict[name]['area'] = constarea
        rdict[name]['sl'] = constsld * numpy.gradient(z) * constarea
        rdict[name]['sld'] = constsld
        return rdict

    def fnWriteGroup2File(self, fp, cName, z):
        raise NotImplementedError()

    def fnWriteGroup2Dict(self, rdict, cName, z):
        # Same as fnWriteGroup2File, but output is saved in a dictionary
        raise NotImplementedError()

    def fnWriteResults2Dict(self, rdict, cName):
        # empty return function, children might implement their own results (combined parameters to return)
        return rdict

    def fnWriteProfile(self, z, aArea=None, anSL=None):
        # Philosophy for this first method: You simply add more and more volume and nSLD to the
        # volume and nSLD array. After all objects have filled up those arrays the maximal area is
        # determined which is the area per molecule and unfilled volume is filled with bulk solvent.
        # Hopefully the fit algorithm finds a physically meaningful solution. There has to be a global
        # hydration parameter for the bilayer.
        # Returns maximum area

        # Do we want a 0.5 * stepsize shift? I believe refl1d FunctionalLayer uses
        # z = numpy.linspace(0, dimension * stepsize, dimension, endpoint=False)
        # z, aArea, anSL must be numpy arrays with the same shape

        # TODO implement wrapping
        # TODO implement absorption

        aArea = numpy.zeros_like(z) if aArea is None else aArea
        anSL = numpy.zeros_like(z) if anSL is None else anSL

        assert (aArea.shape == z.shape)
        assert (anSL.shape == z.shape)

        area, nsl, _ = self.fnGetProfiles(z)
        dMaxArea = area.max()
        aArea += area
        anSL += nsl

        return dMaxArea, aArea, anSL

    def fnWriteProfile2File(self, f, cName, z):
        header = f"z{cName} a{cName} nsl{cName}"
        area, nsl, _ = self.fnGetProfiles(z)
        A = numpy.vstack((z, area, nsl)).T
        numpy.savetxt(f, A, fmt='%0.6e', delimiter=' ', comments='', header=header)
        f.write('\n')
        # TODO: implement wrapping

    def fnWriteProfile2Dict(self, rdict, z):
        area, sl, sld = self.fnGetProfiles(z)
        rdict['zaxis'] = z
        rdict['area'] = area
        rdict['sl'] = sl
        rdict['sld'] = sld
        return rdict

    def fnOverlayProfile(self, z, aArea, anSL, dMaxArea):
        def overlay_profiles(dMaxArea, area1, area2, nsl1, nsl2):
            """
            Overlays area2 (and nsl1) onto area1 (nsl2), replacing area1 if the total area is larger than dMaxArea.
            """
            temparea = area1 + area2

            # find any area for which the final area will be greater than dMaxArea
            overmax = temparea > dMaxArea
            # note: unphysical overfill will trigger the following assertion error
            # TODO: implement a correction instead of throwing an error
            #assert (not numpy.any((temparea - dMaxArea) > area1))
            nsl1[overmax] = nsl1[overmax] * (1 - ((temparea[overmax] - dMaxArea) / area1[overmax])) + nsl2[overmax]
            area1[overmax] = dMaxArea

            # deal with area for which the final area is not greater than dMaxArea
            area1[~overmax] += area2[~overmax]
            nsl1[~overmax] += nsl2[~overmax]

            return area1, nsl1

        # z, aArea, anSL must be numpy arrays with the same shape
        assert (aArea.shape == z.shape)
        assert (anSL.shape == z.shape)

        # TODO implement wrapping
        # TODO implement absorption

        area, nsl, _ = self.fnGetProfiles(z)

        return overlay_profiles(dMaxArea, aArea, area, anSL, nsl)


class CompositenSLDObj(nSLDObj):
    def __init__(self, items=None,  **kwargs):
        super().__init__(**kwargs)
        if items is not None:
            self.subgroups = items
        else:
            # abstract object without subgroups
            self.subgroups = []
        self.nf = 1.0

    def __getitem__(self, index):
        return self.subgroups[index]

    def __setitem__(self, index, item):
        self.subgroups[index] = item

    def __contains__(self, item):
        return item in self.subgroups

    def fnFindSubgroups(self):
        # this should be run at the end of init. Could also store keys if memory is an issue
        # and then use self.__dict__[key]
        # Could use a "subgroup adder" function that adds to subgroups
        self.subgroups = [getattr(self, attr) for attr in dir(self) if isinstance(getattr(self, attr), nSLDObj)]

    def fnSetBulknSLD(self, bulknsld):
        super().fnSetBulknSLD(bulknsld)
        for g in self.subgroups:
            g.fnSetBulknSLD(bulknsld)

    def fnGetVolume(self, z1=None, z2=None, recalculate=True):
        volume = 0.0
        for g in self.subgroups:
            volume += g.fnGetVolume(z1, z2, recalculate)

        return volume * self.nf

    def fnGetProfiles(self, z):
        area = numpy.zeros_like(z)
        nsl = numpy.zeros_like(z)
        nsld = numpy.zeros_like(z)

        for g in self.subgroups:
            newarea, newnsl, _ = g.fnGetProfiles(z)
            area += newarea
            nsl += newnsl

        pos = (area > 0)
        nsld[pos] = nsl[pos] / (area[pos] * numpy.gradient(z)[pos])

        # store profile for post-processing such as derived parameters
        self.zaxis = z
        self.area = area * self.nf
        self.sl = nsl * self.nf
        self.sld = nsld

        return self.area, self.sl, self.sld

    def fnGetnSL(self):
        nSL = 0.
        for group in self.subgroups:
            nSL += group.fnGetnSL()
        return nSL

    def fnWriteGroup2File(self, f, cName, z):
        for g in self.subgroups:
            # this allows objects with the same names from multiple bilayers
            g.fnWriteGroup2File(f, f"{cName}.{g.name}", z)

    def fnWriteGroup2Dict(self, rdict, cName, z):
        for g in self.subgroups:
            # this allows objects with the same names from multiple bilayers
            rdict = g.fnWriteGroup2Dict(rdict, f"{cName}.{g.name}", z)
        return rdict

    def fnWriteResults2Dict(self, rdict, cName):
        for g in self.subgroups:
            rdict = g.fnWriteResults2Dict(rdict, f"{cName}.{g.name}")
        return rdict


class Box2Err(nSLDObj):
    def __init__(self, dz=20, dsigma1=2, dsigma2=2, dlength=10, dvolume=10, dnSL=0, dnumberfraction=1, name=None):
        super().__init__(name=name)
        self.z = dz
        self.sigma1 = dsigma1
        self.sigma2 = dsigma2
        self.length = dlength
        self.init_l = dlength
        self.vol = dvolume
        self.nSL = dnSL
        self.nf = dnumberfraction
        self.nSL2 = None
        self.flip = False
        self.flipcenter = 0.

    @staticmethod
    def _flip_shift(arr, z, flipcenter):
        """
        Flips array and shifts it back that self.z is at same position. Used to invert groups from outer to inner
        leaflet.
        """
        ret = numpy.flip(arr)
        shiftvalue = -1 * (z[-1] - (flipcenter - z[0]) - flipcenter) / (z[1] - z[0])
        ret = shift(ret, shiftvalue, mode='constant')
        return ret

    def fnGetVolume(self, z1=None, z2=None, recalculate=True):
        if (z1 is not None) & (z2 is not None) & recalculate:

            def antiderivative(z):
                # antiderivative of error function as currently implemented
                # https://nvlpubs.nist.gov/nistpubs/jres/73b/jresv73bn1p1_a1b.pdf
                # checked against scipy.integrate.cumtrapz
                ad = (z - self.z + 0.5 * self.length) * erf((z - self.z + 0.5 * self.length) / (numpy.sqrt(2) * self.sigma1)) + 1.0 / (numpy.sqrt(numpy.pi) / (numpy.sqrt(2) * self.sigma1)) * numpy.exp(-((z - self.z + 0.5 * self.length) / (numpy.sqrt(2) * self.sigma1)) ** 2)
                ad -= (z - self.z - 0.5 * self.length) * erf((z - self.z - 0.5 * self.length) / (numpy.sqrt(2) * self.sigma2)) + 1.0 / (numpy.sqrt(numpy.pi) / (numpy.sqrt(2) * self.sigma2)) * numpy.exp(-((z - self.z - 0.5 * self.length) / (numpy.sqrt(2) * self.sigma2)) ** 2)
                ad *= (self.vol / self.length) * 0.5
                ad += self.vol * 0.5

                return ad

            return abs(antiderivative(z2) - antiderivative(z1)) * self.nf
        
        else:
            
            return super().fnGetVolume(z1, z2, recalculate)

    def fnGetProfiles(self, z):
        # calculate area
        # Gaussian function definition, integral is volume, return value is area at positions z
        if (self.length != 0) and (self.sigma1 != 0) and (self.sigma2 != 0):
            area = erf((z - self.z + 0.5 * self.length) / (numpy.sqrt(2) * self.sigma1))
            area -= erf((z - self.z - 0.5 * self.length) / (numpy.sqrt(2) * self.sigma2))
            area *= (self.vol / self.length) * 0.5
            area *= self.nf
        else:
            area = numpy.zeros_like(z)

        # calculate nSLD
        nsld = self.fnGetnSL() / self.vol * numpy.ones_like(z) if self.vol != 0 else numpy.zeros_like(z)
        # calculate nSL.
        nsl = area * nsld * numpy.gradient(z)

        if self.flip:
            area = self._flip_shift(area, z, self.flipcenter)
            nsl = self._flip_shift(nsl, z, self.flipcenter)
            nsld = self._flip_shift(nsld, z, self.flipcenter)

        self.zaxis = z
        self.area = area
        self.sl = nsl
        self.sld = nsld

        return area, nsl, nsld

    def fnGetnSL(self):
        if self.bProtonExchange & (self.bulknsld is not None):
            if self.vol != 0:
                return ((self.bulknsld - H2O_SLD) * self.nSL2 + (D2O_SLD - self.bulknsld) * self.nSL) / (D2O_SLD -
                                                                                                         H2O_SLD)
            else:
                return 0.
        else:
            return self.nSL

    def fnGetnSLD(self, z):
        _, _, nsld = self.fnGetProfiles(z)
        return nsld

    # Gaussians are cut off below and above 3 sigma double
    def fnGetLowerLimit(self):
        return self.z - 0.5 * self.length - 3 * self.sigma1

    def fnGetUpperLimit(self):
        return self.z + 0.5 * self.length + 3 * self.sigma2

    # 7/6/2021 new feature: only use proton exchange if nSL2 is explicitly set
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

    def fnSet(self, volume=None, length=None, position=None, nSL=None, sigma=None, nf=None):
        if volume is not None:
            self.vol = volume
        if length is not None:
            self.length = length
        if position is not None:
            self.fnSetZ(position)
        if nSL is not None:
            nSL2 = None
            if isinstance(nSL, (list, tuple, numpy.ndarray)):
                nSL1 = nSL[0]
                if len(nSL) == 2:
                    nSL2 = nSL[1]
            else:
                nSL1 = nSL
            self.fnSetnSL(nSL1, nSL2)
        if sigma is not None:
            sigma2 = None
            if isinstance(sigma, (list, tuple, numpy.ndarray)):
                sigma1 = sigma[0]
                if len(sigma) == 2:
                    sigma2 = sigma[1]
            else:
                sigma1 = sigma
            self.fnSetSigma(sigma1, sigma2)
        if nf is not None:
            self.nf = nf

    def fnWriteGroup2File(self, fp, cName, z):
        if self.flip:
            position = 2 * self.flipcenter + - self.z
        else:
            position = self.z
        fp.write(f"Box2Err {cName} z {position} sigma1 {self.sigma1} sigma2 {self.sigma2} l {self.length} vol "
                 f"{self.vol} nSL {self.nSL} nSL2 {self.nSL2} nf {self.nf} flip {self.flip}\n")
        self.fnWriteProfile2File(fp, cName, z)

    def fnWriteGroup2Dict(self, rdict, cName, z):
        if self.flip:
            position = 2 * self.flipcenter + - self.z
        else:
            position = self.z
        rdict[cName] = {}
        rdict[cName]['header'] = f"{self.__class__.__name__} {cName} z {position} sigma1 {self.sigma1} " \
                                 f"sigma2 {self.sigma2} l {self.length} vol {self.vol} " \
                                 f"nSL {self.nSL} nSL2 {self.nSL2} nf {self.nf} flip {self.flip}"
        rdict[cName]['z'] = self.z
        rdict[cName]['sigma1'] = self.sigma1
        rdict[cName]['sigma2'] = self.sigma2
        rdict[cName]['l'] = self.length
        rdict[cName]['vol'] = self.vol
        rdict[cName]['nSL'] = self.nSL
        rdict[cName]['nSL2'] = self.nSL2
        rdict[cName]['nf'] = self.nf
        rdict[cName] = self.fnWriteProfile2Dict(rdict[cName], z)
        return rdict

    def fnWriteResults2Dict(self, rdict, cName):
        rdict = super().fnWriteResults2Dict(rdict, cName)
        if cName not in rdict:
            rdict[cName] = {}
        rdict[cName]['COM'] = self.z
        rdict[cName]['INT'] = self.vol
        return rdict


class ComponentBox(Box2Err):
    """
    Box2Err from a components.Component object
    -
    If present, the diffmolecule will be subtracted from molecule. This is needed when, for example, subtracting
    methyl groups from lipid tails to obtain the methylene group
    -
    The arguments molecule and diffmolecule can either be a single instance of Molecule or a list, but the same number
    of elements must be given each. If given as a list, an average length is calculated.
    """
    def __init__(self, components=None, diffcomponents=None, xray_wavelength=None, **kwargs):
        if not isinstance(components, (list, tuple)):
            components = [components]
        if diffcomponents is not None and not isinstance(diffcomponents, (list, tuple)):
            diffcomponents = [diffcomponents]

        dvolume = sum(m.cell_volume for m in components)
        dlength = sum(m.length for m in components) / float(len(components))
        nSL = sum(m.fnGetnSL(xray_wavelength) for m in components)
        if diffcomponents is not None:
            dvolume -= sum(m.cell_volume for m in diffcomponents)
            dlength -= sum(m.length for m in diffcomponents) / float(len(diffcomponents))
            nSL -= sum(m.fnGetnSL(xray_wavelength) for m in diffcomponents)

        super().__init__(dvolume=dvolume, dlength=dlength, **kwargs)
        self.fnSetnSL(*nSL)


class CompositeHeadgroup(CompositenSLDObj):
    def __init__(self, name='headgroup', components=None, innerleaflet=False, xray_wavelength=None,
                 sigma1=None, sigma2=None, rel_pos=None, length=None, position=None, num_frac=None, **kwargs):
        super().__init__(name=name, **kwargs)
        self.flip = innerleaflet
        self.xray_wavelength = xray_wavelength

        if not isinstance(components, (list, tuple)):
            components = [components]

        self.components = []
        for i, component in enumerate(components):
            comp_name = component.name
            # check if componnent is being used more than once, and update name
            if components.count(component) > 1:
                comp_name += f"_{i}"
            comp_obj = ComponentBox(name=comp_name, components=[component], xray_wavelength=self.xray_wavelength)
            self.__setattr__(comp_name, comp_obj)
            self.components.append(comp_obj)

        if length is None:
            self.length = 10.0
        else:
            self.length = length
        self.init_l = self.length

        if position is None:
            self.z = 0.
        else:
            self.z = position

        if num_frac is None:
            self.nf = 1.0
        else:
            self.nf = num_frac

        # will be set in fnAdjustParameters()
        self.vol = 0.

        if sigma1 is None:
            self.sigma1 = numpy.full(len(components), 2.0)
        else:
            self.sigma1 = numpy.array(sigma1)

        if sigma2 is None:
            self.sigma2 = numpy.full(len(components), 2.0)
        else:
            self.sigma2 = numpy.array(sigma2)

        if rel_pos is None:
            self.rel_pos = numpy.linspace(0, self.length, len(components))
        else:
            self.rel_pos = rel_pos

        assert (len(self.sigma1) == len(components)), 'Number of sigma1 values must equal number of components'
        assert (len(self.sigma2) == len(components)), 'Number of sigma2 values must equal number of components'
        assert (len(self.rel_pos) == len(components)), 'Number of rel_pos values must equal number of components'

        self.fnFindSubgroups()
        self.fnAdjustParameters()

    def fnAdjustParameters(self):
        # make sure no group is larger than the entire length of the headgroup
        for g in self.subgroups:
            g.length = min(self.init_l, g.length)

        vol = 0.
        for i, component in enumerate(self.components):
            # Set rel_pos for first and last group to 0.0 and 1.0 respectively if they are supposed to sit
            # flush with the end of the headgroup
            # if i == 0:
            #     pos = 0.5 * component.l
            # elif i == len(self.components)-1:
            #     pos = self.l - 0.5 * component.l
            if self.rel_pos[i] * self.length < component.length * 0.5:
                pos = 0.5 * component.length
            elif self.rel_pos[i] * self.length > self.length - component.length * 0.5:
                pos = self.length - 0.5 * component.length
            else:
                pos = self.rel_pos[i] * self.length
            component.z = self.z - 0.5 * self.length + pos
            component.sigma1 = self.sigma1[i]
            component.sigma2 = self.sigma2[i]
            vol += component.vol
            if self.flip:
                component.flip = True
                component.flipcenter = self.z

        self.vol = vol

    def fnGetLowerLimit(self):
        return self.z - 0.5 * self.length

    def fnGetUpperLimit(self):
        return self.z + 0.5 * self.length

    def fnSet(self, length=None, rel_pos=None, position=None, num_frac=None, bulknsld=None):
        if length is not None:
            self.length = length
        if rel_pos is not None:
            self.rel_pos = rel_pos
        if position is not None:
            self.z = position
        if num_frac is not None:
            self.nf = num_frac
        if bulknsld is not None:
            self.fnSetBulknSLD(bulknsld)
        self.fnAdjustParameters()

    def fnSetSigma(self, sigma=2.0, sigma1=None, sigma2=None):
        # either fill all sigmas with sigma or provide arrays of sigma1 and sigma2

        if (sigma1 is not None) and (sigma2 is not None):
            if len(sigma1) == len(self.sigma1):
                self.sigma1 = numpy.array(sigma1)
            else:
                print('Length of fnSetSigma attribute sigma1 does not match number of subgroups.')
            if len(sigma2) == len(self.sigma2):
                self.sigma2 = numpy.array(sigma2)
            else:
                print('Length of fnSetSigma attribute sigma2 does not match number of subgroups.')
        else:
            self.sigma1.fill(sigma)
            self.sigma2.fill(sigma)

        self.fnAdjustParameters()

    def fnSetZ(self, dz):
        self.z = dz
        self.fnAdjustParameters()

    def fnWriteGroup2File(self, fp, cName, z):
        prefix = "m" if self.flip else ""
        fp.write(f"{prefix}CompositeHeadgroup {cName} z {self.z} l {self.length} vol {self.vol} nf {self.nf} \n")
        self.fnWriteProfile2File(fp, cName, z)
        super().fnWriteGroup2File(fp, cName, z)

    def fnWriteGroup2Dict(self, rdict, cName, z):
        prefix = "m" if self.flip else ""
        rdict[cName] = {}
        rdict[cName]['header'] = f"{prefix}CompositeHeadgroup {cName} z {self.z} l {self.length} vol {self.vol} nf " \
                                 f"{self.nf}"
        rdict[cName]['z'] = self.z
        rdict[cName]['l'] = self.length
        rdict[cName]['vol'] = self.vol
        rdict[cName]['nf'] = self.nf
        rdict[cName] = self.fnWriteProfile2Dict(rdict[cName], z)
        rdict = super().fnWriteGroup2Dict(rdict, cName, z)
        return rdict

# ------------------------------------------------------------------------------------------------------
# Bilayer models
# ------------------------------------------------------------------------------------------------------


class BLM(CompositenSLDObj):
    def __init__(self, inner_lipids=None, inner_lipid_nf=None, outer_lipids=None, outer_lipid_nf=None, lipids=None,
                 lipid_nf=None, xray_wavelength=None, name='blm', **kwargs):
        """
        Free bilayer object. Requires:
            o lipids definition:
                - lipids: list of components.Lipid objects. If set, creates a symmetric bilayer and overrides
                          inner_lipids and outer_lipids. If not set, both inner_lipids and outer_lipids are required.
                - inner_lipids: list of components.Lipid objects. Ignored if lipids is set, required if outer_lipids is
                                set.
                - outer_lipids: list of components.Lipid objects. Ignored if lipids is set, required if inner_lipids is
                                set.
            o number fraction vector: lipid_nf, inner_lipid_nf, outer_lipid_nf: a list of number fractions (not
                                      necessarily normalized) of equal length to 'lipids', 'inner_lipids', or
                                      'outer_lipids', respectively.
        To use an xray probe, set xray_wavelength to the appropriate value in Angstroms.
        """

        super().__init__(name=name, **kwargs)

        # symmetric bilayers can provide only one list of lipids and number fractions
        if lipids is not None:
            inner_lipids = lipids
            inner_lipid_nf = lipid_nf
            outer_lipids = lipids
            outer_lipid_nf = lipid_nf

        assert ((inner_lipids is not None) & (outer_lipids is not None)), 'Inner and outer leaflet lipids must be set.'
        assert len(inner_lipids) == len(inner_lipid_nf), \
            'List of inner lipids and number fractions must be of equal length, not %i and %i' % (len(inner_lipids),
                                                                                                  len(inner_lipid_nf))
        assert len(inner_lipids) > 0, 'Must specify at least one inner lipid'
        assert len(outer_lipids) == len(outer_lipid_nf), \
            'List of inner lipids and number fractions must be of equal length, not %i and %i' % (len(outer_lipids),
                                                                                                  len(outer_lipid_nf))
        assert len(outer_lipids) > 0, 'Must specify at least one inner lipid'

        # normalize number fractions. This allows ratios of lipids to be given instead of number fractions
        self.inner_lipid_nf = numpy.array(inner_lipid_nf) / numpy.sum(inner_lipid_nf)
        self.outer_lipid_nf = numpy.array(outer_lipid_nf) / numpy.sum(outer_lipid_nf)

        self.inner_lipids = inner_lipids
        self.outer_lipids = outer_lipids
        self.xray_wavelength = xray_wavelength

        # unpack lipids
        h, nh, m, mm = self._unpack_lipids(inner_lipids, 'headgroup1', 'methylene1', 'methyl1', innerleaflet=True)
        self.headgroups1 = h
        self.null_hg1 = nh
        self.methylenes1 = m
        self.methyls1 = mm
        h, nh, m, mm = self._unpack_lipids(outer_lipids, 'headgroup2', 'methylene2', 'methyl2', innerleaflet=False)
        self.headgroups2 = h
        self.null_hg2 = nh
        self.methylenes2 = m
        self.methyls2 = mm

        self.initial_hg1_lengths = numpy.array([hg1.length for hg1 in self.headgroups1])
        self.defect_hydrocarbon = Box2Err(name='defect_hc')
        self.defect_headgroup = Box2Err(name='defect_hg')
        self.nf = 1.
        self.vf_bilayer = 1.0
        self.absorb = 0.
        self.l_lipid1 = 11.
        self.l_lipid2 = 11.
        self.normarea = 60.
        self.startz = 50.
        self.sigma = 2.
        # allow for different methyl sigmas if accessed directly (not providing interface)
        self.methyl_sigma = numpy.full_like(self.inner_lipid_nf, 2.)
        self.radius_defect = 100.
        self.hc_substitution_1 = 0
        self.hc_substitution_2 = 0

        self._calc_av_hg()
        self.initial_hg1_l = self.av_hg1_l

        self.fnFindSubgroups()
        self.fnSetBulknSLD(-0.56e-6)
        self.fnAdjustParameters()

    def _adjust_outer_lipids(self):
        self.vol_methyl_outer, self.nsl_methyl_outer = self._unpack_component_pars(self.methyls2)
        self.vol_methylene_outer, self.nsl_methylene_outer = self._unpack_component_pars(self.methylenes2)
        self.l_lipid2 = max(self.l_lipid2, 0.01)

        # make sure number fractions are zero or greater and normalize to one
        self.outer_lipid_nf = self.outer_lipid_nf.clip(min=0.)
        self.outer_lipid_nf /= numpy.sum(self.outer_lipid_nf)

        self.vf_bilayer = max(self.vf_bilayer, 1E-5)
        self.vf_bilayer = min(self.vf_bilayer, 1)

        # TODO: Check if still appropriate here after splitting _adjust_lipids
        self._calc_av_hg()

        # outer hydrocarbons
        self.l_ohc = self.l_lipid2
        self.nf_ohc_lipid = self.outer_lipid_nf
        self.V_ohc = numpy.sum(self.nf_ohc_lipid * self.vol_methylene_outer)
        self.nsl_ohc = numpy.sum(self.nf_ohc_lipid * self.nsl_methylene_outer)
        self.normarea = self.V_ohc / self.l_ohc
        c_s_ohc = self.vf_bilayer
        c_A_ohc = 1
        c_V_ohc = 1
        for i, methylene in enumerate(self.methylenes2):
            methylene.length = self.l_ohc
            methylene.nf = self.nf_ohc_lipid[i] * c_s_ohc * c_A_ohc * c_V_ohc

        # outer methyls
        self.nf_om_lipid = self.nf_ohc_lipid
        self.V_om = numpy.sum(self.nf_om_lipid * self.vol_methyl_outer)
        self.l_om = self.l_ohc * self.V_om / self.V_ohc
        self.nsl_om = numpy.sum(self.nf_om_lipid * self.nsl_methyl_outer)
        c_s_om = c_s_ohc
        c_A_om = 1
        c_V_om = 1
        for i, methyl in enumerate(self.methyls2):
            methyl.length = self.l_om
            methyl.nf = self.nf_om_lipid[i] * c_s_om * c_A_om * c_V_om

        # headgroup size and position
        for hg2, nf_ohc in zip(self.headgroups2, self.nf_ohc_lipid):
            hg2.nf = c_s_ohc * c_A_ohc * nf_ohc * (1 - self.hc_substitution_2)
            if hasattr(hg2, 'fnAdjustParameters'):
                hg2.fnAdjustParameters()

    def _adjust_inner_lipids(self):
        self.vol_methylene_inner, self.nsl_methylene_inner = self._unpack_component_pars(self.methylenes1)
        self.vol_methyl_inner, self.nsl_methyl_inner = self._unpack_component_pars(self.methyls1)
        self.l_lipid1 = max(self.l_lipid1, 0.01)
        # make sure number fractions are zero or greater and normalize to one
        self.inner_lipid_nf = self.inner_lipid_nf.clip(min=0.)
        self.inner_lipid_nf /= numpy.sum(self.inner_lipid_nf)

        # inner hydrocarbons
        self.l_ihc = self.l_lipid1
        self.nf_ihc_lipid = self.inner_lipid_nf
        self.V_ihc = numpy.sum(self.nf_ihc_lipid * self.vol_methylene_inner)
        self.nsl_ihc = numpy.sum(self.nf_ihc_lipid * self.nsl_methylene_inner)
        c_s_ihc = self.vf_bilayer
        c_A_ihc = self.normarea * self.l_ihc / self.V_ihc
        c_V_ihc = 1
        for i, methylene in enumerate(self.methylenes1):
            methylene.length = self.l_ihc
            methylene.nf = self.nf_ihc_lipid[i] * c_s_ihc * c_A_ihc * c_V_ihc

        # inner methyl
        self.nf_im_lipid = self.nf_ihc_lipid
        self.V_im = numpy.sum(self.nf_im_lipid * self.vol_methyl_inner)
        self.l_im = self.l_ihc * self.V_im / self.V_ihc
        self.nsl_im = numpy.sum(self.nf_im_lipid * self.nsl_methyl_inner)
        c_s_im = c_s_ihc
        c_A_im = c_A_ihc
        c_V_im = 1
        for i, methyl in enumerate(self.methyls1):
            methyl.length = self.l_im
            methyl.nf = self.nf_im_lipid[i] * c_s_im * c_A_im * c_V_im

        # headgroup size and position
        for hg1, nf_ihc, in zip(self.headgroups1, self.nf_ihc_lipid):
            hg1.nf = c_s_ihc * c_A_ihc * nf_ihc * (1 - self.hc_substitution_1)
            if hasattr(hg1, 'fnAdjustParameters'):
                hg1.fnAdjustParameters()

        # TODO: Check on usage
        self._calc_av_hg()

    def _adjust_defects(self):
        hclength = self.l_ihc + self.l_im + self.l_om + self.l_ohc
        hglength = self.av_hg1_l + self.av_hg2_l

        if self.radius_defect < (0.5 * (hclength + hglength)):
            self.radius_defect = 0.5 * (hclength + hglength)

        volhalftorus = numpy.pi ** 2 * (self.radius_defect - (2. * hclength / 3. / numpy.pi)) * hclength * hclength / 4.
        volcylinder = numpy.pi * self.radius_defect * self.radius_defect * hclength
        defectarea = volhalftorus / volcylinder * (1 - self.vf_bilayer) * self.normarea

        self.defect_hydrocarbon.vol = defectarea * hclength
        self.defect_hydrocarbon.length = hclength
        self.defect_hydrocarbon.z = self.z_ihc - 0.5 * self.l_ihc + 0.5 * hclength
        self.defect_hydrocarbon.nSL = (self.nsl_ohc + self.nsl_om) / (self.V_ohc + self.V_om) *\
            self.defect_hydrocarbon.vol
        self.defect_hydrocarbon.fnSetSigma(self.sigma)
        self.defect_hydrocarbon.nf = 1

        defectratio = self.defect_hydrocarbon.vol / self.V_ohc
        self.defect_headgroup.vol = defectratio * numpy.sum([hg.nf * hg.vol for hg in self.headgroups2])
        self.defect_headgroup.length = hclength + hglength
        self.defect_headgroup.z = self.z_ihc - 0.5 * self.l_ihc - 0.5 * self.av_hg1_l + 0.5 * (hclength + hglength)
        self.defect_headgroup.nSL = defectratio * numpy.sum([hg.nf * hg.fnGetnSL() for hg in self.headgroups2])
        self.defect_headgroup.fnSetSigma(self.sigma)
        self.defect_headgroup.nf = 1

    def _adjust_z(self, startz):
        # startz is the position of the hg1/lipid1 interface.
        # change here: use average headgroup length instead of a specific headgroup
        self.z_ihc = startz + 0.5 * self.l_ihc
        self.z_im = self.z_ihc + 0.5 * (self.l_ihc + self.l_im)
        self.z_om = self.z_im + 0.5 * (self.l_im + self.l_om)
        self.z_ohc = self.z_om + 0.5 * (self.l_om + self.l_ohc)

        for m1, m2 in zip(self.methylenes1, self.methylenes2):
            m1.fnSetZ(self.z_ihc)
            m2.fnSetZ(self.z_ohc)

        for m1, m2 in zip(self.methyls1, self.methyls2):
            m1.fnSetZ(self.z_im)
            m2.fnSetZ(self.z_om)

        for hg1, hg2 in zip(self.headgroups1, self.headgroups2):
            hg1.fnSetZ(self.z_ihc - 0.5 * self.l_ihc - 0.5 * hg1.length)
            hg2.fnSetZ(self.z_ohc + 0.5 * self.l_ohc + 0.5 * hg2.length)

    def _calc_av_hg(self):
        # calculate average headgroup lengths, ignore zero volume (e.g. cholesterol)
        self.av_hg1_l = numpy.sum(numpy.array([hg.length for hg, use in zip(self.headgroups1, ~self.null_hg1) if use]) *
                                  self.inner_lipid_nf[~self.null_hg1]) / numpy.sum(self.inner_lipid_nf[~self.null_hg1])
        self.av_hg2_l = numpy.sum(numpy.array([hg.length for hg, use in zip(self.headgroups2, ~self.null_hg2) if use]) *
                                  self.outer_lipid_nf[~self.null_hg2]) / numpy.sum(self.outer_lipid_nf[~self.null_hg2])

    @staticmethod
    def _unpack_component_pars(components):
        n_components = len(components)
        vol_components = numpy.zeros(n_components)
        nsl_components = numpy.zeros(n_components)

        for i, component in enumerate(components):
            vol_components[i] = component.vol
            nsl_components[i] = component.fnGetnSL()

        return vol_components, nsl_components

    def _unpack_lipids(self, _lipids, hgprefix, methyleneprefix, methylprefix, innerleaflet=True):
        """ Helper function for BLM classes that unpacks a lipid list into headgroup objects
            and lists of acyl chain and methyl volumes and nSLs. Creates the following attributes:
            o headgroups1: a list of inner leaflet headgroup objects
            o headgroups2: a list of outer leaflet headgroup objects
            NB: each object in headgroups1 and headgroups2 is created as a standalone attribute in "self"
                so that searching for all nSLDObj in self will return each headgroup individually as
                headgroup1_1, headgroup1_2, ... headgroup1_n for n lipids. Similarly, for headgroup2.
            o vol_acyl_lipids: a numpy array of acyl chain volumes
            o vol_methyl_lipids: a numpy array of methyl volumes
            o nsl_acyl_lipids: a numpy array of acyl chain nSL
            o nsl_methyl_lipids: a numpy array of methyl nSL
            """
        # create objects for each lipid headgroup and a composite tail object
        headgroups = []
        methylenes = []
        methyls = []

        for i, lipid in enumerate(_lipids):
            hg_name = f"{hgprefix}_{i+1}"
            methylene_name = f"{methyleneprefix}_{i+1}"
            methyl_name = f"{methylprefix}_{i+1}"

            if isinstance(lipid.headgroup, cmp.Component):
                # populates nSL, nSL2, vol, and l
                hg_obj = ComponentBox(name=hg_name, components=[lipid.headgroup], xray_wavelength=self.xray_wavelength)
            elif isinstance(lipid.headgroup, list):
                hg_obj = lipid.headgroup[0](name=hg_name, innerleaflet=innerleaflet,
                                            xray_wavelength=self.xray_wavelength, **(lipid.headgroup[1]))
            else:
                raise TypeError('Lipid.hg must be a Headgroup object or a subclass of CompositenSLDObj')

            methylene_obj = ComponentBox(name=methylene_name, components=lipid.tails, diffcomponents=lipid.methyls,
                                         xray_wavelength=self.xray_wavelength)
            methyl_obj = ComponentBox(name=methyl_name, components=lipid.methyls, xray_wavelength=self.xray_wavelength)

            self.__setattr__(hg_name, hg_obj)
            headgroups.append(self.__getattribute__(hg_name))
            self.__setattr__(methylene_name, methylene_obj)
            methylenes.append(self.__getattribute__(methylene_name))
            self.__setattr__(methyl_name, methyl_obj)
            methyls.append(self.__getattribute__(methyl_name))

        # find null headgroups to exclude from averaging over headgroup properties
        null_hg = numpy.array([hg.vol <= 0.0 for hg in headgroups], dtype=bool)

        return headgroups, null_hg, methylenes, methyls

    def fnAdjustParameters(self):
        self._adjust_outer_lipids()
        self._adjust_inner_lipids()
        self._adjust_z(self.startz + self.av_hg1_l)
        self._adjust_defects()
        self.fnSetSigma(self.sigma)

    # return value of center of the membrane
    def fnGetCenter(self):
        return self.methyls1[0].z + 0.5 * self.methyls1[0].length

    # Use limits of molecular subgroups
    def fnGetLowerLimit(self):
        return min([hg.fnGetLowerLimit() for hg in self.headgroups1])

    def fnGetUpperLimit(self):
        return max([hg.fnGetUpperLimit() for hg in self.headgroups2])

    def fnSet(self, sigma=2.0, bulknsld=-0.56e-6, startz=20., l_lipid1=11., l_lipid2=11.0, vf_bilayer=1.0,
              nf_inner_lipids=None, nf_outer_lipids=None, nf_lipids=None, hc_substitution_1=0.,
              hc_substitution_2=0., radius_defect=100., **kwargs):
        self.sigma = sigma
        self.fnSetBulknSLD(bulknsld)
        self.startz = startz
        self.l_lipid1 = l_lipid1
        self.l_lipid2 = l_lipid2
        self.vf_bilayer = vf_bilayer

        if nf_lipids is not None:
            nf_inner_lipids = nf_outer_lipids = nf_lipids
        if nf_inner_lipids is not None:  # pass None to keep lipid_nf unchanged
            assert len(nf_inner_lipids) == len(self.inner_lipids), \
                'nf_lipids must be same length as number of lipids in bilayer, not %i and %i' % (len(nf_inner_lipids),
                                                                                                 len(self.inner_lipids))
            self.inner_lipid_nf = numpy.array(nf_inner_lipids)
        if nf_outer_lipids is not None:  # pass None to keep lipid_nf unchanged
            assert len(nf_outer_lipids) == len(self.outer_lipids), \
                'nf_lipids must be same length as number of lipids in bilayer, not %i and %i' % (len(nf_outer_lipids),
                                                                                                 len(self.outer_lipids))
            self.outer_lipid_nf = numpy.array(nf_outer_lipids)

        self.hc_substitution_1 = hc_substitution_1
        self.hc_substitution_2 = hc_substitution_2
        self.radius_defect = radius_defect

        self.fnAdjustParameters()

    def fnSetSigma(self, sigma):
        self.sigma = sigma
        for hg1, hg2 in zip(self.headgroups1, self.headgroups2):
            hg1.fnSetSigma(sigma)
            hg2.fnSetSigma(sigma)

        for i, (methylene1, methyl1, methyl2, methylene2) in enumerate(zip(self.methylenes1, self.methyls1,
                                                                           self.methyls2, self.methylenes2)):
            methylene1.fnSetSigma(sigma, numpy.sqrt(sigma ** 2 + self.methyl_sigma[i] ** 2))
            methyl1.fnSetSigma(numpy.sqrt(sigma ** 2 + self.methyl_sigma[i] ** 2), numpy.sqrt(sigma ** 2 +
                                                                                           self.methyl_sigma[i] ** 2))
            methyl2.fnSetSigma(numpy.sqrt(sigma ** 2 + self.methyl_sigma[i] ** 2), numpy.sqrt(sigma ** 2 +
                                                                                           self.methyl_sigma[i] ** 2))
            methylene2.fnSetSigma(numpy.sqrt(sigma ** 2 + self.methyl_sigma[i] ** 2), sigma)

        self.defect_hydrocarbon.fnSetSigma(sigma)
        self.defect_headgroup.fnSetSigma(sigma)

    def fnWriteGroup2File(self, fp, cName, z):
        super().fnWriteGroup2File(fp, cName, z)
        self.fnWriteConstant(fp, f"{cName}_normarea", self.normarea, 0, z)

    def fnWriteGroup2Dict(self, rdict, cName, z):
        rdict = super().fnWriteGroup2Dict(rdict, cName, z)
        rdict = self.fnWriteConstant2Dict(rdict, f"{cName}.normarea", self.normarea, 0, z)
        return rdict

    def fnWriteResults2Dict(self, rdict, cName):
        rdict = super().fnWriteResults2Dict(rdict, cName)
        if cName not in rdict:
            rdict[cName] = {}
        rdict[cName]['area_per_lipid'] = self.normarea
        rdict[cName]['volume_fraction'] = self.vf_bilayer
        rdict[cName]['thickness_inner_leaflet'] = self.l_ihc + self.l_im
        rdict[cName]['thickness_outer_leaflet'] = self.l_ohc + self.l_om
        rdict[cName]['thickness_total'] = self.l_ohc + self.l_om + self.l_ihc + self.l_im

        if self.normarea != 0:
            p2 = self.headgroups1[0].z - 0.5 * self.headgroups1[0].length
            p3 = self.methylenes1[0].z - 0.5 * self.methylenes1[0].length
            p5 = self.methylenes2[0].z + 0.5 * self.methylenes2[0].length
            p6 = self.headgroups2[0].z + 0.5 * self.headgroups2[0].length

            rdict[cName]['water in inner headgroups'] = self.fnGetVolume(p2, p3, recalculate=False)
            rdict[cName]['water in inner headgroups'] /= (self.normarea * (p3 - p2))
            rdict[cName]['water in inner headgroups'] = 1 - rdict[cName]['water in inner headgroups']
            rdict[cName]['water in hydrocarbons'] = self.fnGetVolume(p3, p5, recalculate=False)
            rdict[cName]['water in hydrocarbons'] /= (self.normarea * (p5 - p3))
            rdict[cName]['water in hydrocarbons'] = 1 - rdict[cName]['water in hydrocarbons']
            rdict[cName]['water in outer headgroups'] = self.fnGetVolume(p5, p6, recalculate=False)
            rdict[cName]['water in outer headgroups'] /= (self.normarea * (p6 - p5))
            rdict[cName]['water in outer headgroups'] = 1 - rdict[cName]['water in outer headgroups']

        return rdict

'''
class Monolayer(CompositenSLDObj, BLM):
    def __init__(self, lipids=None, lipid_nf=None, xray_wavelength=None, name='monolayer', **kwargs):
        """
        Free monolayer object. Requires:
            o lipids definition: lipids: list of components.Lipid objects. If set, creates a symmetric bilayer and
            overrides inner_lipids and outer_lipids.
            o number fraction vector: lipid_nf, : a list of number fractions (not necessarily normalized) of equal
            length to 'lipids'.
        To use an xray probe, set xray_wavelength to the appropriate value in Angstroms.
        """

        super().__init__(name=name, **kwargs)

        # symmetric bilayers can provide only one list of lipids and number fractions
        outer_lipids = lipids
        outer_lipid_nf = lipid_nf
        assert len(outer_lipids) > 0, 'Must specify at least one inner lipid'

        # normalize number fractions. This allows ratios of lipids to be given instead of number fractions
        self.outer_lipid_nf = numpy.array(outer_lipid_nf) / numpy.sum(outer_lipid_nf)
        self.outer_lipids = outer_lipids
        self.xray_wavelength = xray_wavelength

        # unpack lipids
        h, nh, m, mm = self._unpack_lipids(outer_lipids, 'headgroup2', 'methylene2', 'methyl2', innerleaflet=False)
        self.headgroups2 = h
        self.null_hg2 = nh
        self.methylenes2 = m
        self.methyls2 = mm

        self.nf = 1.
        self.vf_bilayer = 1.0
        self.absorb = 0.
        self.l_lipid2 = 11.
        self.normarea = 60.
        self.startz = 50.
        self.sigma = 2.
        # allow for different methyl sigmas if accessed directly (not providing interface)
        self.methyl_sigma = numpy.full_like(self.inner_lipid_nf, 2.)
        self.hc_substitution_2 = 0

        self._calc_av_hg()

        self.fnFindSubgroups()
        self.fnSetBulknSLD(-0.56e-6)
        self.fnAdjustParameters()

    def _adjust_z(self, startz):
        # startz is the position of the air/methyl2 interface.
        # change here: use average headgroup length instead of a specific headgroup
        self.z_om = startz - 0.5 * self.l_om
        self.z_ohc = self.z_om + 0.5 * (self.l_om + self.l_ohc)
        for m2 in self.methylenes2:
            m2.fnSetZ(self.z_ohc)
        for m2 in self.methyls2:
            m2.fnSetZ(self.z_om)
        for hg2 in self.headgroups2:
            hg2.fnSetZ(self.z_ohc + 0.5 * self.l_ohc + 0.5 * hg2.length)

    def _calc_av_hg(self):
        # calculate average headgroup lengths, ignore zero volume (e.g. cholesterol)
        self.av_hg2_l = numpy.sum(numpy.array([hg.length for hg, use in zip(self.headgroups2, ~self.null_hg2) if use]) *
                                  self.outer_lipid_nf[~self.null_hg2]) / numpy.sum(self.outer_lipid_nf[~self.null_hg2])

    def fnAdjustParameters(self):
        self._adjust_outer_lipids()
        self._adjust_z(self.startz)
        self.fnSetSigma(self.sigma)

    # return value of center of the membrane
    def fnGetCenter(self):
        return self.methyls2[0].z - 0.5 * self.methyls2[0].length + 0.5 * \
            (self.methyls2[0].length + self.methylenes2[0].length + self.avg_hg2_l)

    # Use limits of molecular subgroups
    def fnGetLowerLimit(self):
        return self.methyls2[0].z - 0.5 * self.methyls2[0].length

    def fnGetUpperLimit(self):
        return max([hg.fnGetUpperLimit() for hg in self.headgroups2])

    def fnSet(self, sigma=2.0, bulknsld=-0.56e-6, startz=20., l_lipid2=11.0, vf_bilayer=1.0,
              nf_lipids=None, hc_substitution_2=0., **kwargs):
        self.sigma = sigma
        self.fnSetBulknSLD(bulknsld)
        self.startz = startz
        self.l_lipid2 = l_lipid2
        self.vf_bilayer = vf_bilayer

        if nf_lipids is not None:
            nf_outer_lipids = nf_lipids
            self.outer_lipid_nf = numpy.array(nf_outer_lipids)

        self.hc_substitution_2 = hc_substitution_2

        self.fnAdjustParameters()

    def fnSetSigma(self, sigma):
        self.sigma = sigma
        for hg2 in self.headgroups2:
            hg2.fnSetSigma(sigma)

        for i, (methyl2, methylene2) in enumerate(zip(self.methyls2, self.methylenes2)):
            methyl2.fnSetSigma(numpy.sqrt(sigma ** 2 + self.methyl_sigma[i] ** 2),
                               numpy.sqrt(sigma ** 2 + self.methyl_sigma[i] ** 2))
            methylene2.fnSetSigma(numpy.sqrt(sigma ** 2 + self.methyl_sigma[i] ** 2), sigma)

    def fnWriteGroup2File(self, fp, cName, z):
        BLM.fnWriteGroup2File(self, fp, cName, z)

    def fnWriteGroup2Dict(self, rdict, cName, z):
        return BLM.fnWriteGroup2Dict(self, rdict, cName, z)

    def fnWriteResults2Dict(self, rdict, cName):
        rdict = super().fnWriteResults2Dict(rdict, cName)
        if cName not in rdict:
            rdict[cName] = {}
        rdict[cName]['area_per_lipid'] = self.normarea
        rdict[cName]['volume_fraction'] = self.vf_bilayer
        rdict[cName]['thickness_lipid_leaflet'] = self.l_ohc + self.l_om

        if self.normarea != 0:
            p3 = self.methyls2[0].z - 0.5 * self.methyls2[0].length
            p5 = self.methylenes2[0].z + 0.5 * self.methylenes2[0].length
            p6 = self.headgroups2[0].z + 0.5 * self.headgroups2[0].length

            rdict[cName]['water in hydrocarbons'] = self.fnGetVolume(p3, p5, recalculate=False)
            rdict[cName]['water in hydrocarbons'] /= (self.normarea * (p5 - p3))
            rdict[cName]['water in hydrocarbons'] = 1 - rdict[cName]['water in hydrocarbons']
            rdict[cName]['water in outer headgroups'] = self.fnGetVolume(p5, p6, recalculate=False)
            rdict[cName]['water in outer headgroups'] /= (self.normarea * (p6 - p5))
            rdict[cName]['water in outer headgroups'] = 1 - rdict[cName]['water in outer headgroups']

        return rdict
'''

class Monolayer(BLM):
    def __init__(self, lipids=None, lipid_nf=None, xray_wavelength=None, name='monolayer', **kwargs):
        """
        Monolayer object at an interface. Requires:
            o lipids definition:
                - lipids: list of components.Lipid objects. If set, creates a symmetric bilayer and overrides
                    inner_lipids and outer_lipids. If not set, both inner_lipids and outer_lipids are required.
            o number fraction vector: lipid_nf: a list of number fractions (not necessarily normalized) of equal length
                to 'lipids'.
        To use an xray probe, set xray_wavelength to the appropriate value in Angstroms.
        """

        super().__init__(inner_lipids=lipids, inner_lipid_nf=lipid_nf, outer_lipids=lipids, outer_lipid_nf=lipid_nf,
                         lipids=lipids, lipid_nf=lipid_nf, xray_wavelength=xray_wavelength, name=name, **kwargs)
        
        self.methyl_sigma = numpy.zeros_like(self.outer_lipid_nf)

    def fnGetCenter(self):
        return self.methyls2[0].z - 0.5 * self.methyls2[0].length

    def _adjust_defects(self):
        hclength = self.l_om + self.l_ohc
        hglength = self.av_hg2_l

        if self.radius_defect < (0.5 * (hclength + hglength)):
            self.radius_defect = 0.5 * (hclength + hglength)

        volhalftorus = numpy.pi ** 2 * (self.radius_defect - (2. * hclength / 3. / numpy.pi)) * hclength * hclength / 4.
        volcylinder = numpy.pi * self.radius_defect * self.radius_defect * hclength
        defectarea = volhalftorus / volcylinder * (1 - self.vf_bilayer) * self.normarea

        self.defect_hydrocarbon.vol = defectarea * hclength
        self.defect_hydrocarbon.length = hclength
        self.defect_hydrocarbon.z = self.z_ohc + 0.5 * self.l_ohc - 0.5 * hclength
        self.defect_hydrocarbon.nSL = (self.nsl_ohc + self.nsl_om) / (self.V_ohc + self.V_om) *\
            self.defect_hydrocarbon.vol
        self.defect_hydrocarbon.fnSetSigma(self.sigma)
        self.defect_hydrocarbon.nf = 1

        defectratio = self.defect_hydrocarbon.vol / self.V_ohc
        self.defect_headgroup.vol = defectratio * numpy.sum([hg.nf * hg.vol for hg in self.headgroups2])
        self.defect_headgroup.length = hclength + hglength
        self.defect_headgroup.z = self.z_ohc + 0.5 * self.l_ihc + 0.5 * self.av_hg2_l - 0.5 * (hclength + hglength)
        self.defect_headgroup.nSL = defectratio * numpy.sum([hg.nf * hg.fnGetnSL() for hg in self.headgroups2])
        self.defect_headgroup.fnSetSigma(self.sigma)
        self.defect_headgroup.nf = 1

    def fnAdjustParameters(self):
        self._adjust_outer_lipids()
        self._adjust_inner_lipids()

        # turn off inner lipids
        for gp in self.headgroups1 + self.methyls1 + self.methylenes1:
            gp.nf = 0.0

        # reference this to outer leaflet hydrophobic interface
        self._adjust_z(self.startz - self.l_ihc - self.l_im - self.l_ohc - self.l_om)
        self._adjust_defects()
        self.fnSetSigma(self.sigma)

    def fnSet(self, startz = 20.0, sigma=2.0, bulknsld=-0.56e-6, l_lipid2=11.0, vf_bilayer=1.0,
              nf_outer_lipids=None, nf_lipids=None, hc_substitution_2=0., radius_defect=100.):

        if startz is not None:
            self.startz = startz
        if sigma is not None:
            self.sigma = sigma
        if bulknsld is not None:
            self.fnSetBulknSLD(bulknsld)
        if l_lipid2 is not None:
            self.l_lipid2 = l_lipid2
        if vf_bilayer is not None:
            self.vf_bilayer = vf_bilayer

        if nf_lipids is not None:
            nf_outer_lipids = nf_lipids
        if nf_outer_lipids is not None:  # pass None to keep lipid_nf unchanged
            assert len(nf_outer_lipids) == len(self.outer_lipids), \
                'nf_lipids must be same length as number of lipids in bilayer, not %i and %i' % (len(nf_outer_lipids),
                                                                                                 len(self.outer_lipids))
            self.outer_lipid_nf = numpy.array(nf_outer_lipids)

        self.hc_substitution_2 = hc_substitution_2
        self.radius_defect = radius_defect

        self.fnAdjustParameters()

    def fnWriteResults2Dict(self, rdict, cName):
        rdict = super().fnWriteResults2Dict(rdict, cName)
        if cName not in rdict:
            rdict[cName] = {}

        rdict[cName]['area_per_lipid'] = self.normarea
        rdict[cName]['volume_fraction'] = self.vf_bilayer
        rdict[cName]['hydrocarbon_thickness'] = self.l_ohc + self.l_om

        if self.normarea != 0:
            p3 = self.fnGetCenter()
            p5 = self.methylenes2[0].z + 0.5 * self.methylenes2[0].length
            p6 = self.headgroups2[0].z + 0.5 * self.headgroups2[0].length

            #rdict[cName]['water in hydrocarbons'] = self.fnGetVolume(p3, p5, recalculate=False)
            #rdict[cName]['water in hydrocarbons'] /= (self.normarea * (p5 - p3))
            #rdict[cName]['water in hydrocarbons'] = 1 - rdict[cName]['water in hydrocarbons']
            rdict[cName]['water in outer headgroups'] = self.fnGetVolume(p5, p6, recalculate=False)
            rdict[cName]['water in outer headgroups'] /= (self.normarea * (p6 - p5))
            rdict[cName]['water in outer headgroups'] = 1 - rdict[cName]['water in outer headgroups']

        return rdict

class ssBLM(BLM):
    def __init__(self, inner_lipids=None, inner_lipid_nf=None, outer_lipids=None, outer_lipid_nf=None, lipids=None,
                 lipid_nf=None, xray_wavelength=None, name='ssblm', **kwargs):
        """
        Solid supported bilayer object. Requires:
            o lipids definition:
                - lipids: list of components.Lipid objects. If set, creates a symmetric bilayer and overrides
                    inner_lipids and outer_lipids. If not set, both inner_lipids and outer_lipids are required.
                - inner_lipids: list of components.Lipid objects. Ignored if lipids is set, required if outer_lipids
                    is set.
                - outer_lipids: list of components.Lipid objects. Ignored if lipids is set, required if inner_lipids
                    is set.
            o number fraction vector: lipid_nf, inner_lipid_nf, outer_lipid_nf: a list of number fractions
                (not necessarily normalized) of equal length to 'lipids', 'inner_lipids', or 'outer_lipids',
                respectively.
        To use an xray probe, set xray_wavelength to the appropriate value in Angstroms.
        """

        # add ssBLM-specific subgroups
        self.substrate = Box2Err(name='substrate')
        self.substrate.length = 40
        self.substrate.z = 0
        self.substrate.nf = 1
        self.rho_substrate = 2.07e-6
        self.l_siox = 1
        self.rho_siox = 3.55e-6

        self.siox = Box2Err(name='siox')
        self.siox.length = 20
        self.siox.z = self.substrate.z + self.substrate.length / 2.0 + self.siox.length / 2.0
        self.siox.nf = 1

        self.l_submembrane = 10.
        self.global_rough = 2.0

        super().__init__(inner_lipids, inner_lipid_nf, outer_lipids=outer_lipids, outer_lipid_nf=outer_lipid_nf,
                         lipids=lipids, lipid_nf=lipid_nf, xray_wavelength=xray_wavelength, name=name, **kwargs)

    def _adjust_substrate(self):
        self.substrate.vol = self.normarea * self.substrate.length
        self.substrate.nSL = self.rho_substrate * self.substrate.vol
        self.siox.length = self.l_siox
        self.siox.vol = self.normarea * self.siox.length
        self.siox.nSL = self.rho_siox * self.siox.vol
        self.siox.z = self.substrate.z + self.substrate.length / 2 + 0.5 * self.siox.length

    def fnAdjustParameters(self):
        self._adjust_outer_lipids()
        self._adjust_inner_lipids()
        self._adjust_substrate()
        self._adjust_z(self.substrate.z + self.substrate.length / 2 + self.siox.length + self.l_submembrane +
                       self.av_hg1_l)
        self._adjust_defects()
        self.fnSetSigma(self.sigma)

    def fnSetSigma(self, sigma):
        super().fnSetSigma(sigma)
        self.substrate.fnSetSigma(self.global_rough)
        self.siox.fnSetSigma(self.global_rough)

    def fnGetLowerLimit(self):
        return self.substrate.fnGetLowerLimit()
        # does this make sense since this goes to negative z and isn't intended to be used?

    def fnSet(self, global_rough=2.0, rho_substrate=2.07e-6, rho_siox=3.55e-6, l_siox=20, l_submembrane=10, **kwargs):
        self.global_rough = global_rough
        self.rho_substrate = rho_substrate
        self.rho_siox = rho_siox
        self.l_siox = l_siox
        self.l_submembrane = l_submembrane
        super().fnSet(**kwargs)

    def fnWriteResults2Dict(self, rdict, cName):
        rdict = super().fnWriteResults2Dict(rdict, cName)
        if cName not in rdict:
            rdict[cName] = {}

        if self.normarea != 0:
            p1 = self.siox.z + 0.5 * self.siox.length
            p2 = self.headgroups1[0].z - 0.5 * self.headgroups1[0].length
            rdict[cName]['water in submembrane'] = self.fnGetVolume(p1, p2, recalculate=False)
            rdict[cName]['water in submembrane'] /= (self.normarea * (p2 - p1))
            rdict[cName]['water in submembrane'] = 1 - rdict[cName]['water in submembrane']

        return rdict


class tBLM(BLM):
    def __init__(self, tether, filler, inner_lipids=None, inner_lipid_nf=None, outer_lipids=None, outer_lipid_nf=None,
                 lipids=None, lipid_nf=None, xray_wavelength=None, name='tblm', **kwargs):
        """
        Tethered lipid bilayer. Requires:

        o inner_lipids, outer_lipids: a list of components.Lipid objects. See BLM documentation.
        o inner_lipid_nf, outer_lipid_nf: a list of number fractions (not necessarily normalized) of 
                        equal length to 'lipids'
        o tether: a components.Tether object
        o filler: a components.Component object describing the filler molecule, including its length

        If outer_lipids is not specified, inner_lipids is used for both leaflets.

        To use an xray probe, set xray_wavelength to the appropriate value in Angstroms.
        """
        self.substrate = Box2Err(name='substrate')
        self.substrate.length = 40
        self.substrate.z = 0
        self.substrate.nf = 1
        self.rho_substrate = 2.07e-6
        self.global_rough = 2.0

        self.nf_tether = 0.3
        self.l_tether = 8.0
        self.mult_tether = 7.0 / 3.0
        self.tether_methyl_sigma = 2.0
        self.bME = ComponentBox(name='bME', components=filler, xray_wavelength=xray_wavelength)
        self.initial_bME_l = self.bME.length
        self.tether_bme = Box2Err(name='tether_bme')
        self.tether_free = Box2Err(name='tether_free')
        self.tether_hg = Box2Err(name='tether_hg')
        self.tether = cmp.Component(name='tether', formula=tether.tether.formula,
                                    cell_volume=tether.tether.cell_volume, xray_wavelength=xray_wavelength,
                                    length=self.l_tether)
        self.tetherg = cmp.Component(name='tetherg', formula=tether.tetherg.formula,
                                     cell_volume=tether.tetherg.cell_volume, xray_wavelength=xray_wavelength,
                                     length=10.0)
        self.tether_methylene = ComponentBox(name='tether_methylene', components=tether.tails,
                                             diffcomponents=tether.methyls, xray_wavelength=xray_wavelength)
        self.tether_methyl = ComponentBox(name='tether_methyl', components=tether.methyls,
                                          xray_wavelength=xray_wavelength)

        super().__init__(inner_lipids=inner_lipids, inner_lipid_nf=inner_lipid_nf, outer_lipids=outer_lipids,
                         outer_lipid_nf=outer_lipid_nf, lipids=lipids, lipid_nf=lipid_nf,
                         xray_wavelength=xray_wavelength, name=name, **kwargs)

    def _adjust_inner_lipids(self):
        self.vol_methylene_inner, self.nsl_methylene_inner = self._unpack_component_pars(self.methylenes1)
        self.vol_methyl_inner, self.nsl_methyl_inner = self._unpack_component_pars(self.methyls1)
        self.vol_methylene_tether = self.tether_methylene.vol
        self.vol_methyl_tether = self.tether_methyl.vol
        self.nsl_methylene_tether = self.tether_methylene.fnGetnSL()
        self.nsl_methyls_tether = self.tether_methyl.fnGetnSL()

        self.l_lipid1 = max(self.l_lipid1, 0.01)
        # make sure number fractions are zero or greater and normalize to one
        self.inner_lipid_nf = self.inner_lipid_nf.clip(min=0.)
        self.inner_lipid_nf /= numpy.sum(self.inner_lipid_nf)

        # inner hydrocarbons
        self.l_ihc = self.l_lipid1
        self.nf_ihc_lipid = self.inner_lipid_nf * (1 - self.nf_tether)
        self.nf_ihc_tether = self.nf_tether
        self.V_ihc = numpy.sum(self.nf_ihc_lipid * self.vol_methylene_inner) + self.nf_ihc_tether * \
                     self.vol_methylene_tether
        self.nsl_ihc = numpy.sum(self.nf_ihc_lipid * self.nsl_methylene_inner) + self.nf_ihc_tether * \
                  self.nsl_methylene_tether
        c_s_ihc = self.vf_bilayer
        c_A_ihc = self.normarea * self.l_ihc / self.V_ihc
        c_V_ihc = 1
        self.tether_methylene.length = self.l_ihc
        self.tether_methylene.nf = self.nf_ihc_tether * c_s_ihc * c_A_ihc * c_V_ihc
        for i, methylene in enumerate(self.methylenes1):
            methylene.length = self.l_ihc
            methylene.nf = self.nf_ihc_lipid[i] * c_s_ihc * c_A_ihc * c_V_ihc

        # inner methyl
        self.nf_im_lipid = self.nf_ihc_lipid
        self.nf_im_tether = self.nf_ihc_tether
        self.V_im = numpy.sum(self.nf_im_lipid * self.vol_methyl_inner) + self.nf_im_tether * self.vol_methyl_tether
        self.l_im = self.l_ihc * self.V_im / self.V_ihc
        self.nsl_im = numpy.sum(self.nf_im_lipid * self.nsl_methyl_inner) + self.nf_im_tether * self.nsl_methyls_tether
        c_s_im = c_s_ihc
        c_A_im = c_A_ihc
        c_V_im = 1
        self.tether_methyl.length = self.l_im
        self.tether_methyl.nf = self.nf_im_tether * c_s_im * c_A_im * c_V_im
        for i, methyl in enumerate(self.methyls1):
            methyl.length = self.l_im
            methyl.nf = self.nf_im_lipid[i] * c_s_im * c_A_im * c_V_im

        # headgroup size and position
        for hg1, nf_ihc, in zip(self.headgroups1, self.nf_ihc_lipid):
            hg1.nf = c_s_ihc * c_A_ihc * nf_ihc * (1 - self.hc_substitution_1)
            if hasattr(hg1, 'fnAdjustParameters'):
                hg1.fnAdjustParameters()

        self._calc_av_hg()

        # tether glycerol part
        c_s_tg = c_s_ihc
        c_A_tg = c_A_ihc
        c_V_tg = self.nf_ihc_tether
        # self.tetherg.l = self.tetherg.cell_volume / ((self.vol_acyl_tether -
        # self.vol_methyl_tether) / self.lipid1.l) / 0.9
        self.tetherg_nf = c_s_tg * c_A_tg * c_V_tg
        self.tetherg.nf = self.tetherg_nf

        # tether EO part
        c_s_EO = c_s_ihc
        c_A_EO = c_A_ihc
        c_V_EO = self.nf_ihc_tether
        self.tether_nf = c_s_EO * c_A_EO * c_V_EO
        self.tether.nf = self.tether_nf

    def _adjust_submembrane(self):
        """
        Approach:
        1. Divide submembrane space into three regions: tether_bme, tether_free, and tether_hg
        2. Attempt to evenly distribute density of everything in these three regions.
        3. Algorithm:
            Case 1: sparse (l_tether > BME.l + av_hg1_l): -> l_tether_free>0
            Case 2: compact (l_tether < BME.l + av_hg1_l):
                -> Proportionately compress BME.l and av_hg1_l, so typically 1/3 BME compression, 2/3 hg compression,
                such that l_tether=BME.l + av_hg1_l -> l_tether_free=0
            Then:
                a. Calculate A_bme (tether in bME plus bME itself) assuming min tether area equals that of bME
                b. Calculate A_hg (tether_g in hg plus hgs themselves) assuming min tether area is that from tether_g
                c. Fill up  A_tether_bme, A_tether_free, A_tether_hg with remaining tether volume using a bucket filling
                   algorithm such that A_tether_bme + A_bme = A_tether_free = A_tether_hg + A_hg + A_tetherg, if enough
                   tether volume remains
        """

        def _adjust_mult_tether():
            """
            bME + tether in the bME region can't be bigger than normarea. If it is, reduce mult_tether until it isn't.
            """
            # minimum amount of area in tether region is assumed to be equal to that of bMe
            min_A_tether_bme = self.tether_nf * self.bME.vol / self.bME.length
            # area in bme region
            A_bme = self.mult_tether * self.tether_nf * self.bME.vol / self.bME.length + min_A_tether_bme
            if A_bme > self.normarea * self.vf_bilayer:
                self.mult_tether = max(0, (self.normarea * self.vf_bilayer - min_A_tether_bme) /
                                       (self.tether_nf * self.bME.vol / self.bME.length))
                # area in bme region
                A_bme = self.mult_tether * self.tether_nf * self.bME.vol / self.bME.length + min_A_tether_bme
            return A_bme, min_A_tether_bme

        # Total volume of tether molecules, including the glycerol
        V_tether = self.tether.cell_volume * self.tether_nf + self.tetherg.cell_volume * self.tetherg_nf
        total_submembrane_V = V_tether + self.mult_tether * self.tether_nf * self.bME.vol
        total_submembrane_V += numpy.sum([hg.nf * hg.vol for hg in self.headgroups1])
        # If there is too much volume in the submembrane space, increase l_tether to accommodate it
        if total_submembrane_V > self.l_tether * self.normarea * self.vf_bilayer:
            self.l_tether = total_submembrane_V / (self.normarea * self.vf_bilayer)

        # Reset bME length and headgroup length
        # TODO: find a robust way to keep track of these initial values. If a user changes them between init and this
        #  point, the user value will be overwritten unless they specifically change initial_bME_l as well. This is
        #  already done for headgroups (fnSetHeadgroupLength)
        self.bME.length = self.initial_bME_l
        for hg1, initial_length in zip(self.headgroups1, self.initial_hg1_lengths):
            hg1.length = initial_length
        self._calc_av_hg()

        # If too much bME is present, adjust the mult_tether parameter
        # TODO: I do not think that this is good. Mult_tether is a parameter that is conserved across multiple
        #  data sets that might vary in sub-membrane thickness. This approach breaks this link. We should think of
        #  a different way how the structure reacts to an overfilling of the bMe region. (F.H.)
        A_bme, min_A_tether_bme = _adjust_mult_tether()
        l_tether_free = self.l_tether - (self.initial_bME_l + self.av_hg1_l)

        # TODO: Currently bMe and headgroups are squished to a maximum such that no hydration water remains in either
        #   group. This is unrealistic. Also, composite PC headgroups do not have a homogeneous space filling, and the
        #   current algorithm leads to an overfilling of the available volume in regions where the composite PC area
        #   is above average.
        if l_tether_free < 0:
            # squish headgroups and bME proportionately to their existing size.
            d1s = self.l_tether - (self.initial_bME_l + self.initial_hg1_lengths)
            self.bME.length += l_tether_free * self.initial_bME_l / (self.initial_bME_l + self.initial_hg1_l)
            for hg1, d1, initial_length in zip(self.headgroups1, d1s, self.initial_hg1_lengths):
                d1 = self.l_tether - (self.bME.length + initial_length)
                # only squish headgroups if individual length is too big.
                hg1.length = initial_length + min(d1, 0.0)
            self._calc_av_hg()
            A_bme, min_A_tether_bme = _adjust_mult_tether()
            l_tether_free = 0

        # Calculate minimum area that has to reside in the headgroup region
        # NOTE: right now this is just the tetherg volume. It should probably be larger (at least 2 EO groups). This
        # can be adjusted in the Tether molecule so tetherg has more volume.
        # TODO: Not good. The glycerol is a small group with an area close to that of the tether double
        #   chain. This is why its area traditionally had been modeled as a fixed fraction of the double chain, allowing
        #   for some water. The remedy to include EO volume in the tether_g group is not preferable as it does not
        #   divide the molecule based on chemistry/scattering properties but on later usage. (F.H.)
        min_A_tether_hg = self.tetherg.cell_volume * self.tetherg_nf / self.av_hg1_l
        A_hg = numpy.sum([hg.nf * hg.vol / hg.length for hg in self.headgroups1]) + min_A_tether_hg

        V_tether_bme = min_A_tether_bme * self.bME.length
        V_tether_hg = min_A_tether_hg * self.av_hg1_l
        V_tether_remainder = V_tether - V_tether_bme - V_tether_hg

        bucket_vol = numpy.array([(self.normarea-A_bme+min_A_tether_bme) * self.bME.length, self.normarea * l_tether_free,
                                  (self.normarea-A_hg+min_A_tether_hg) * self.av_hg1_l])
        bucket_fill = numpy.array([min_A_tether_bme * self.bME.length, 0, min_A_tether_hg * self.av_hg1_l])
        V_tether_remainder, bucket_fill = self._fill_bucket(bucket_vol, V_tether_remainder, bucket_fill)
        V_tether_bme, V_tether_free, V_tether_hg = bucket_fill

        self.tether_bme.fnSet(volume=V_tether_bme, length=self.bME.length,
                              position=0.5 * self.bME.length + self.substrate.z + self.substrate.length * 0.5,
                              nSL=self.tether.nSLs / self.tether.cell_volume * V_tether_bme)
        self.tether_free.fnSet(volume=V_tether_free, length=l_tether_free,
                               position=self.tether_bme.z + 0.5 * self.tether_bme.length + 0.5 * l_tether_free,
                               nSL=self.tether.nSLs / self.tether.cell_volume * V_tether_free)
        frac_tether = 1 - self.tetherg.cell_volume / V_tether_hg
        self.tether_hg.fnSet(volume=V_tether_hg, length=self.av_hg1_l, position=self.tether_free.z + 0.5 *
                             l_tether_free + 0.5 * self.av_hg1_l, nSL=self.tether.nSLs /
                             self.tether.cell_volume * frac_tether * V_tether_hg + self.tetherg.nSLs)
        self.bME.fnSet(position=self.tether_bme.z, nf=self.mult_tether * self.tether_nf)

    def _adjust_substrate(self):
        self.substrate.vol = self.normarea * self.substrate.length
        self.substrate.nSL = self.rho_substrate * self.substrate.vol

    def _adjust_z(self, startz):
        # startz is the position of the hg1/lipid1 interface.
        self.z_ihc = self.substrate.length * 0.5 + self.l_tether + 0.5 * self.l_ihc
        self.z_im = self.z_ihc + 0.5 * (self.l_ihc + self.l_im)
        self.z_om = self.z_im + 0.5 * (self.l_im + self.l_om)
        self.z_ohc = self.z_om + 0.5 * (self.l_om + self.l_ohc)

        for m1, m2 in zip(self.methylenes1, self.methylenes2):
            m1.fnSetZ(self.z_ihc)
            m2.fnSetZ(self.z_ohc)

        for m1, m2 in zip(self.methyls1, self.methyls2):
            m1.fnSetZ(self.z_im)
            m2.fnSetZ(self.z_om)

        for hg1, hg2 in zip(self.headgroups1, self.headgroups2):
            hg1.fnSetZ(self.z_ihc - 0.5 * self.l_ihc - 0.5 * hg1.length)
            hg2.fnSetZ(self.z_ohc + 0.5 * self.l_ohc + 0.5 * hg2.length)

        self.tether_methylene.fnSetZ(self.z_ihc)
        self.tether_methyl.fnSetZ(self.z_im)

    @staticmethod
    def _fill_bucket(bucket_volume, volume, fill_level=None):
        """
        Bucket filling algorithm for n buckets with bucket_volume and total liquid volume to fill.
        Buckets can be prefilled. First underfilled buckets are filled up to the level of the
        prefilled ones. Supports non-overlapping buckets concerning bucket volume and fill level.
        """
        n_bckts = bucket_volume.shape[0]
        if fill_level is None:
            fill_level = numpy.zeros_like(bucket_volume)
        current_fill_level = fill_level
        minlevel = numpy.amin(current_fill_level)
        filltolevel = minlevel

        while volume > 0:
            fillnow = numpy.zeros_like(bucket_volume)
            for i in range(n_bckts):
                if current_fill_level[i] == minlevel:
                    if bucket_volume[i] > minlevel:
                        fillnow[i] = 1.0
                        if filltolevel == minlevel or bucket_volume[i] < filltolevel:
                            filltolevel = bucket_volume[i]
                if current_fill_level[i] > minlevel:
                    if filltolevel == minlevel or current_fill_level[i] < filltolevel:
                        filltolevel = current_fill_level[i]

            # buckets are full, volume remaining
            if minlevel == filltolevel:
                break
            if numpy.sum(fillnow) != 0:
                filling_volume = min(filltolevel - minlevel, volume / numpy.sum(fillnow))
                current_fill_level += filling_volume * fillnow
                volume -= filling_volume
            minlevel = filltolevel

        return volume, current_fill_level

    def fnAdjustParameters(self):
        self._adjust_outer_lipids()
        self._adjust_inner_lipids()
        self._adjust_submembrane()
        self._adjust_substrate()
        self._adjust_z(self.tether_hg.z + 0.5 * self.tether_hg.length)
        self._adjust_defects()
        self.fnSetSigma(self.sigma)

    def fnSetHeadgroupLength(self, hg, value):
        """
            Sets specific headgroup to length 'value'.
            Does not recalculate values using fnAdjustParameters.
        """
        hg.length = value
        # only for inner leaflets
        if hg in self.headgroups1:
            self.initial_hg1_lengths[self.headgroups1.index(hg)] = value

    def fnSetSigma(self, sigma):
        super().fnSetSigma(sigma)

        tether_methyl_sigma = numpy.sqrt(self.sigma ** 2 + self.tether_methyl_sigma ** 2)
        self.tether_methylene.fnSetSigma(self.sigma, tether_methyl_sigma)
        self.tether_methyl.fnSetSigma(tether_methyl_sigma, tether_methyl_sigma)

        self.substrate.fnSetSigma(self.global_rough)
        if self.tether_free.vol > 0:
            self.bME.fnSetSigma(self.global_rough)
            self.tether_bme.fnSetSigma(self.global_rough)
            self.tether_free.fnSetSigma(self.global_rough, sigma)
            self.tether_hg.fnSetSigma(sigma)
        else:
            self.bME.fnSetSigma(self.global_rough, sigma)
            self.tether_bme.fnSetSigma(self.global_rough, sigma)
            self.tether_hg.fnSetSigma(sigma)

    def fnGetLowerLimit(self):
        return self.substrate.fnGetLowerLimit()
        # does this make sense since this goes to negative z and isn't intended to be used?

    def fnSet(self, global_rough=2.0, rho_substrate=4.55e-6, nf_tether=0.5, mult_tether=3., l_tether=20.,  **kwargs):
        self.global_rough = global_rough
        self.rho_substrate = rho_substrate
        self.nf_tether = nf_tether
        self.mult_tether = mult_tether
        self.l_tether = l_tether
        super().fnSet(**kwargs)

    def fnWriteResults2Dict(self, rdict, cName):
        rdict = super().fnWriteResults2Dict(rdict, cName)
        if cName not in rdict:
            rdict[cName] = {}
        rdict[cName]['thickness_tether'] = self.l_tether

        if self.normarea != 0:
            p1 = self.substrate.z + 0.5 * self.substrate.length
            p2 = self.headgroups1[0].z - 0.5 * self.headgroups1[0].length
            rdict[cName]['water in submembrane'] = self.fnGetVolume(p1, p2, recalculate=False)
            rdict[cName]['water in submembrane'] /= (self.normarea * (p2 - p1))
            rdict[cName]['water in submembrane'] = 1 - rdict[cName]['water in submembrane']

        return rdict


# ------------------------------------------------------------------------------------------------------
# Protein models
# ------------------------------------------------------------------------------------------------------

class ProteinBox(Box2Err):
    """Special Box2Err designed for use with protein densities that are expressed as volume fractions.
        Unlike Box2Err, has a "normarea" attribute; setting this attribute changes only the total volume
        but not the volume fraction nor the length. (For fixed volume objects, use Box2Err.)
    """

    def __init__(self, dz=20, dsigma1=2, dsigma2=2, dlength=10, dvolume_fraction=0, dnSLD_H=0, dnSLD_D=0, protexchratio=1.0, dnumberfraction=1, normarea=1.0, name=None):
        super().__init__(dz, dsigma1, dsigma2, dlength, 0, 0, dnumberfraction, name)
        self.nSLD_H = dnSLD_H
        self.nSLD_D = dnSLD_D
        self.protexchratio = protexchratio
        self.bProtonExchange = True
        self.normarea = normarea
        self._update_volume(dvolume_fraction)

    def _update_volume(self, vf: float):
        """check volume calculation if length, normarea, or volume fraction are changed"""
        self.vol = self.length * self.normarea * vf
    
    def fnGetnSL(self):
        """Overrides base class for simplicity."""

        if self.bProtonExchange & (self.bulknsld is not None):
            fracD = self.protexchratio * (self.bulknsld - H2O_SLD) / (D2O_SLD - H2O_SLD)
            sld = fracD * self.nSLD_D + (1 - fracD) * self.nSLD_H
        else:
            sld = self.nSLD_H

        return self.vol * sld

    def fnSetnSL(self, _nSL, _nSL2=None):
        raise NotImplementedError

    def fnSetNormarea(self, dnormarea):
        old_vf = self.vol / (self.length * self.normarea)
        self.normarea = dnormarea
        self._update_volume(old_vf)

    def fnSet(self, volume_fraction=None, length=None, position=None, nSLD=None, sigma=None, nf=None):
        vf = self.vol / (self.length * self.normarea)
        if volume_fraction is not None:
            vf = volume_fraction
        if length is not None:
            self.length = length
        if position is not None:
            self.fnSetZ(position)
        if nSLD is not None:
            if isinstance(nSLD, (list, tuple, numpy.ndarray)):
                self.nSLD_H = nSLD[0]
                if len(nSLD) == 2:
                    self.nSLD_D = nSLD[1]
            else:
                self.nSLD_H = nSLD
                self.nSLD_D = nSLD
        if sigma is not None:
            sigma2 = None
            if isinstance(sigma, (list, tuple, numpy.ndarray)):
                sigma1 = sigma[0]
                if len(sigma) == 2:
                    sigma2 = sigma[1]
            else:
                sigma1 = sigma
            self.fnSetSigma(sigma1, sigma2)
        if nf is not None:
            self.nf = nf

        self._update_volume(vf)

    def fnWriteGroup2Dict(self, rdict, cName, z):
        rdict = super().fnWriteGroup2Dict(rdict, cName, z)
        rdict[cName]['nSLD_H'] = self.nSLD_H
        rdict[cName]['nSLD_D'] = self.nSLD_D
        return rdict

class TetheredBoxDouble(Box2Err):
    """Special Box2Err that is tethered to, but can move freely around, two z positions.
        Allows specification of the two tether points, the length, and the fractional position.
        The fractional position ranges from 0 to 1; 0 places the rightmost edge of the box at the
        right tether point; 1 places the leftmost edge of the box at the left tether point.

        Intended for use with, e.g. protein loops that are not found in X-ray structures
    """

    def __init__(self, z_tether1=20, z_tether2=20, frac_position=0.5, dsigma1=2, dsigma2=2, dlength=10, dvolume=10, dnSL=0, dnumberfraction=1, name=None):
        super().__init__(0, dsigma1, dsigma2, dlength, dvolume, dnSL, dnumberfraction, name)
        self.z_tether1 = z_tether1
        self.z_tether2 = z_tether2
        self.frac_position = frac_position

        self._update_z()

    def _update_z(self):
        # check order of tether points
        z_tether1 = min(self.z_tether1, self.z_tether2)
        z_tether2 = max(self.z_tether1, self.z_tether2)

        # calculate tether position
        z_tether = z_tether2 - (z_tether2 - z_tether1) * self.frac_position

        # set center z of box
        self.fnSetZ(z_tether + self.length * (self.frac_position - 0.5))

    def fnSet(self, volume=None, length=None, tether_position1=None, tether_position2=None, frac_position=None, nSL=None, sigma=None, nf=None):
        super().fnSet(volume, length, None, nSL, sigma, nf)

        if tether_position1 is not None:
            self.z_tether1 = tether_position1
        if tether_position1 is not None:
            self.z_tether2 = tether_position2
        if frac_position is not None:
            self.frac_position = frac_position

        self._update_z()

    def fnWriteGroup2Dict(self, rdict, cName, z):
        rdict = super().fnWriteGroup2Dict(rdict, cName, z)
        rdict[cName]['z_tether1'] = self.z_tether1
        rdict[cName]['z_tether2'] = self.z_tether2
        rdict[cName]['frac_position'] = self.frac_position

        return rdict

class TetheredBox(TetheredBoxDouble):
    """Special case of TetheredBoxDouble where z_tether1 = z_tether2.

        Intended use is for missing N- or C-terminal domains of proteins
    """

    def __init__(self, z_tether=20, frac_position=0.5, dsigma1=2, dsigma2=2, dlength=10, dvolume=10, dnSL=0, dnumberfraction=1, name=None):
        super().__init__(z_tether, z_tether, frac_position, dsigma1, dsigma2, dlength, dvolume, dnSL, dnumberfraction, name)

    def fnSet(self, volume=None, length=None, tether_position=None, frac_position=None, nSL=None, sigma=None, nf=None):
        super().fnSet(volume, length, tether_position, tether_position, frac_position, nSL, sigma, nf)

    def fnWriteGroup2Dict(self, rdict, cName, z):
        rdict = super().fnWriteGroup2Dict(rdict, cName, z)
        rdict[cName]['z_tether'] = self.z_tether

        return rdict


class Hermite(nSLDObj):
    """
    Hermite splines

    Notes on Hermite usage:
    1. Instantiate the spline object, e.g. h = SLDHermite()
        NB: the only initialization variable that isn't overwritten by fnSetRelative is dnormarea, so the others have
        been removed and dnSLD has been moved to fnSetRelative instead of the initialization. This is to keep the
        fnSetRelative function calls similar between Hermite and SLDHermite (except in Hermite nSLD is a constant and in
        SLDHermite it's a list of control points)
    2. Set all internal parameters, e.g. damping parameters, monotonic spline, etc.
    3. Call fnSetRelative.
        NB: for speed, the spline interpolators are stored in the object. Only fnSetRelative will update them!
    """
    def __init__(self, dnormarea=60, name='hermite'):
        super().__init__(name=name)
        self.numberofcontrolpoints = 10
        self.nSLD = None
        self.normarea = dnormarea
        self.monotonic = True
        self.damping = True
        self.dampthreshold = 0.001
        self.dampFWHM = 0.0002
        self.damptrigger = 0.04
        self.dstartposition = 0.0

        self.dp = numpy.arange(self.numberofcontrolpoints)
        self.vf = numpy.zeros(self.numberofcontrolpoints)
        self.damp = numpy.zeros(self.numberofcontrolpoints)

    def _apply_damping(self):
        dampfactor = numpy.ones_like(self.vf)
        if self.damping:
            # find index of first point beyond damping trigger
            above_damptrigger = numpy.where(self.vf > self.damptrigger)[0]
            # if no points above trigger, damping doesn't apply
            if len(above_damptrigger) > 0:
                # does nothing if peaked point is at the end
                dampfactor[above_damptrigger[0] + 1:] = 1. / (1 + numpy.exp(-2.1 * (self.vf[above_damptrigger[0] + 1:]
                                                                                    - self.dampthreshold) /
                                                                            self.dampFWHM))
                dampfactor = numpy.cumprod(dampfactor)

        self.damp = self.vf * dampfactor

    def _set_area_spline(self):
        self._apply_damping()

        if self.monotonic:  # monotone interpolation
            self.area_spline = PchipInterpolator(self.dp, self.damp, extrapolate=False)
        else:  # catmull-rom
            self.area_spline = self._catmull_rom(self.dp, self.damp, extrapolate=False)

        self.area_spline_integral = self.area_spline.antiderivative()

    @staticmethod
    def _catmull_rom(xctrl, yctrl, **kwargs):
        """returns a catmull_rom spline using ctrl points xctrl, yctrl"""
        assert (len(xctrl) == len(yctrl))  # same shape

        x_xtnd = numpy.insert(xctrl, 0, xctrl[0] - (xctrl[1] - xctrl[0]))
        x_xtnd = numpy.append(x_xtnd, xctrl[-1] - xctrl[-2] + xctrl[-1])

        y_xtnd = numpy.insert(yctrl, 0, yctrl[0])
        y_xtnd = numpy.append(y_xtnd, yctrl[-1])
        dydx = 0.5 * (y_xtnd[2:] - y_xtnd[:-2]) / (x_xtnd[2:] - x_xtnd[:-2])

        return CubicHermiteSpline(xctrl, yctrl, dydx, **kwargs)

    def fnGetProfiles(self, z):
        vf = self.area_spline(z)
        vf[numpy.isnan(vf)] = 0.0
        vf[vf < 0] = 0.0

        self.zaxis = z
        self.area = vf * self.normarea * self.nf
        self.sld = self.fnGetnSLDProfile(z)
        self.sl = self.area * numpy.gradient(z) * self.sld

        return self.area, self.sl, self.sld

    def fnGetnSLDProfile(self, z):
        return self.nSLD * numpy.ones_like(z)

    def fnGetnSLD(self, dz):
        return self.nSLD

    def fnGetLowerLimit(self):
        return self.dp[0]

    def fnGetUpperLimit(self):
        return self.dp[-1]

    def fnGetVolume(self, z1=None, z2=None, recalculate=True):
        if recalculate:
            # use stored antiderivatives to calculate volume
            # make sure z1 and z2 are in the defined interval. If both are above or below, result will be zero
            z1 = max(self.fnGetLowerLimit(), z1)
            z1 = min(self.fnGetUpperLimit(), z1)
            z2 = max(self.fnGetLowerLimit(), z2)
            z2 = min(self.fnGetUpperLimit(), z2)
            return (self.area_spline_integral(z2) - self.area_spline_integral(z1)) * self.nf * self.normarea
        else:
            return super().fnGetVolume(z1=z1, z2=z2, recalculate=recalculate)

    def fnSetNormarea(self, dnormarea):
        self.normarea = dnormarea

    def fnSetnSLD(self, dnSLD):
        self.nSLD = dnSLD

    def fnSetRelative(self, dSpacing, dStart, dDp, dVf, dnSLD, dnf):
        self.fnSetnSLD(dnSLD)
        self.vf = numpy.array(dVf)
        self.numberofcontrolpoints = len(self.vf)
        self.dp = dStart + dSpacing * numpy.arange(self.numberofcontrolpoints) + numpy.array(dDp)
        # make sure the control points have compatible shapes (previous command should fail if not)
        assert (self.vf.shape == self.dp.shape)
        self.dstartposition = dStart
        self.nf = dnf
        self._set_area_spline()

    def fnWriteGroup2File(self, fp, cName, z):
        fp.write(f"Hermite {cName} numberofcontrolpoints {self.numberofcontrolpoints} normarea {self.normarea} nf "
                 f"{self.nf} \n")
        self.fnWriteProfile2File(fp, cName, z)

    def fnWriteGroup2Dict(self, rdict, cName, z):
        rdict[cName] = {}
        rdict[cName]['header'] = f"Hermite {cName} numberofcontrolpoints {self.numberofcontrolpoints} normarea " \
                                 f"{self.normarea} nf {self.nf}"
        rdict[cName]['numberofcontrolpoints'] = self.numberofcontrolpoints
        rdict[cName]['normarea'] = self.normarea
        rdict[cName]['nf'] = self.nf
        rdict[cName] = self.fnWriteProfile2Dict(rdict[cName], z)
        return rdict

    def fnWriteResults2Dict(self, rdict, cName):
        rdict = super().fnWriteResults2Dict(rdict, cName)
        if cName not in rdict:
            rdict[cName] = {}
        # to avoid costly recalculations, the results are derived from self.area stored after evaluating the profiles
        # this requires evaluating the profiles before calculation of the results, which is usually done
        pos = numpy.argmax(self.area)
        # peak_widths often produces a PeakPropertyWarning when the peak prominence is zero
        # we do not want to handle this situation or alert the user to it
        # a meaningless result will be obvious (hopefully)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            results = peak_widths(self.area, [pos], rel_height=0.5)
        rdict[cName]['peak position'] = self.zaxis[pos]
        rdict[cName]['FWHM'] = results[0] * (self.zaxis[1] - self.zaxis[0])
        rdict[cName]['COM'] = numpy.sum(self.area * self.zaxis) / numpy.sum(self.area) if numpy.sum(self.area) != 0 else 0
        rdict[cName]['INT'] = numpy.sum(self.area * numpy.gradient(self.zaxis))
        return rdict

class BoxHermite(Hermite):
    """Special Hermite class where the first control points define the first half of an error function.
        Note that this may behave badly if sigma is more than about 1/2 the control point spacing.
    """

    def __init__(self, dnormarea=60, n_box=21, name='boxhermite'):
        super().__init__(dnormarea, name)

        # width of spline at previous interface
        self.sigma = 2.0

        # set number of control points for box. n = 9 keeps errors below 0.1 % (empirically)
        self.n_box = n_box

    def _get_box_spline(self):
        # finds extra control points corresponding to a generic error function with width sigma
        new_dp = numpy.linspace(-3, 0, self.n_box, endpoint=True) * self.sigma
        new_vf = 0.5 * (erf(new_dp / (self.sigma * numpy.sqrt(2))) + 1)

        return new_dp, new_vf

    def _set_area_spline(self):
        new_dp, new_vf = self._get_box_spline()

        # shift to correct start position
        new_dp += self.dstartposition

        # scale to volume fraction of first control point
        new_vf *= self.vf[0]

        # add new control points to user-defined list
        if len(self.dp) > 1:
            #self.dp = numpy.hstack((new_dp, self.dp[1:]))
            #self.vf = numpy.hstack((new_vf, self.vf[1:]))
            self.vf = numpy.hstack((new_vf[new_dp < (self.dp[1] - self.sigma)], self.vf[1:]))
            self.dp = numpy.hstack((new_dp[new_dp < (self.dp[1] - self.sigma)], self.dp[1:]))

        else:
            self.dp = new_dp
            self.vf = new_vf * self.vf[0]

        return super()._set_area_spline()

    def fnSetRelative(self, dSpacing, dStart, dDp, dVf, dnSLD, dnf, sigma):
        self.sigma = sigma
        return super().fnSetRelative(dSpacing, dStart, dDp, dVf, dnSLD, dnf)
    

class SLDHermite(Hermite):
    def __init__(self, dnormarea, **kwargs):
        super().__init__(dnormarea, **kwargs)
        self.sld = numpy.zeros(self.numberofcontrolpoints)

    def _set_sld_spline(self):
        if self.monotonic:  # monotone interpolation
            self.sld_spline = PchipInterpolator(self.dp, self.sld, extrapolate=False)
        else:  # catmull-rom
            self.sld_spline = self._catmull_rom(self.dp, self.sld, extrapolate=False)

    def fnGetnSLDProfile(self, z):
        sld = self.sld_spline(z)
        # deal with out-of-range values
        sld[z < self.dp[0]] = self.sld[0]
        sld[z > self.dp[-1]] = self.sld[-1]
        assert (not numpy.any(numpy.isnan(sld)))  # shouldn't be any NaNs left

        return sld

    def fnSetnSLD(self, dnSLD):
        self.nSLD = None

    def fnSetRelative(self, dSpacing, dStart, dDp, dVf, dnSLD, dnf):
        self.vf = numpy.array(dVf)
        self.numberofcontrolpoints = len(self.vf)
        self.dp = dStart + dSpacing * numpy.arange(self.numberofcontrolpoints) + numpy.array(dDp)
        # make sure the control points have compatible shapes (previous command should fail if not)
        assert (self.vf.shape == self.dp.shape)
        self.sld = numpy.array(dnSLD)
        assert (self.sld.shape == self.dp.shape)
        self.dstartposition = dStart
        self.nf = dnf
        self._set_area_spline()
        self._set_sld_spline()


class ContinuousEuler(nSLDObj):
    """
    Uses scipy.spatial library to do Euler rotations in real time
    """
    def __init__(self, fn8col, rotcenter=None, xray=False, **kwargs):
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
        self.resnums = resdata[:, 0]
        self.rescoords = resdata[:, 1:4]
        if rotcenter is not None:
            rotcenter = numpy.array(rotcenter)
            assert rotcenter.shape == self.rescoords[0, :].shape
            self.rescoords -= rotcenter
        self.rotcoords = numpy.zeros_like(self.rescoords)
        self.resscatter = resdata[:, 4:]
        self.gamma = 0.
        self.beta = 0.
        self.sigma = 2.
        self.z = 0.
        self.nf = 1.
        self.protexchratio = 1.
        self.xray = xray
        self.fnSetBulknSLD(None)
        self.R = Rotation.from_euler('zy', [self.gamma, self.beta], degrees=True)

        # TODO: Would it make sense to have a concept of "normarea"? Then there could be a "volume fraction" concept
        #  so that max(area) = volume_fraction * normarea

    aa3to1 = dict({'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                   'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                   'THR': 'T', 'TYR': 'Y', 'VAL': 'V', 'TRP': 'W'})

    def _apply_transform(self):
        self.R = Rotation.from_euler('zy', [self.gamma, self.beta], degrees=True)
        self.rotcoords = self.R.apply(self.rescoords)
        self.rotcoords[:, 2] += self.z

    def fnGetProfiles(self, z):
        # self.xray=True for xray profile
        # self.xray=False (default) for a neutron probe; bulknsld determines fraction of exchangeable hydrogens used

        # get area profile
        dz = z[1] - z[0]  # calculate step size (MUST be uniform)
        zbins = numpy.append(z, z[-1] + dz) - 0.5 * dz  # bin edges; z is treated as bin centers here
        # TODO: Smoothing with gaussian_filter requires constant step size. What happens if z is a single point?
        h = numpy.histogram(self.rotcoords[:, 2], bins=zbins, weights=self.resscatter[:, 0])
        area = gaussian_filter(h[0], self.sigma / dz, order=0, mode='constant', cval=0)
        area /= dz

        if self.xray:
            h = numpy.histogram(self.rotcoords[:, 2], bins=zbins, weights=self.resscatter[:, 1])
            nsl = gaussian_filter(h[0], self.sigma / dz, order=0, mode='constant', cval=0)
        else:
            # get nslH profile
            h = numpy.histogram(self.rotcoords[:, 2], bins=zbins, weights=self.resscatter[:, 2])
            nslH = gaussian_filter(h[0], self.sigma / dz, order=0, mode='constant', cval=0)

            if self.bulknsld is not None:
                # get nslD profile
                h = numpy.histogram(self.rotcoords[:, 2], bins=zbins, weights=self.resscatter[:, 3])
                nslD = gaussian_filter(h[0], self.sigma / dz, order=0, mode='constant', cval=0)
                fracD = self.protexchratio * (self.bulknsld - H2O_SLD) / (D2O_SLD - H2O_SLD)
                nsl = fracD * nslD + (1 - fracD) * nslH

            else:
                nsl = nslH

        nsld = numpy.zeros_like(z)
        pos = (area > 0)
        nsld[pos] = nsl[pos] / (area[pos] * numpy.gradient(z)[pos])

        self.area = area * self.nf
        self.zaxis = z
        self.sl = nsl * self.nf
        self.sld = nsld

        return self.area, self.sl, self.sld


    def fnGetVolume(self, z1=None, z2=None, recalculate=True):
        """
        Calculates volume based on the number of residues of the rotated molecule located between
        z positions z1 and z2 (inclusive).
            
        Note: the result is (slightly) different from integrating the area, because the roughness has already
        been applied to the area. However, this remains more accurate than the roughened value because
        the integration limits  will typically correspond to the limits of a lipid Box2Err function
        which are themselves defined before the roughness is applied.
        """
        if z1 is None or z2 is None or recalculate is False:
            return super().fnGetVolume(z1, z2, recalculate)
        else:
            # use a single bin defined by bin edges z1 and z2.
            # First [0] selects the histogram array, second [0] gets the single value from the array
            volume = numpy.histogram(self.rotcoords[:, 2], bins=[z1, z2], weights=self.resscatter[:, 0])[0][0]
            return volume * self.nf

    def fnSet(self, gamma=0., beta=0., zpos=0., sigma=2.0, nf=1.0, bulknsld=None):
        self.gamma = gamma
        self.beta = beta
        self.z = zpos
        self.nf = nf
        self.sigma = sigma
        self.fnSetBulknSLD(bulknsld)
        self._apply_transform()

    def fnWriteGroup2File(self, fp, cName, z):
        fp.write(f"ContinuousEuler {cName} StartPosition {self.z} Gamma {self.gamma} Beta {self.beta} nf {self.nf} \n")
        self.fnWriteProfile2File(fp, cName, z)

    def fnWriteGroup2Dict(self, rdict, cName, z):
        rdict[cName] = {}
        rdict[cName]['header'] = f"ContinuousEuler {cName} StartPosition {self.z} Gamma {self.gamma} Beta {self.beta} nf {self.nf}"
        rdict[cName]['startposition'] = self.z
        rdict[cName]['gamma'] = self.gamma
        rdict[cName]['beta'] = self.beta
        rdict[cName]['nf'] = self.nf
        rdict[cName] = self.fnWriteProfile2Dict(rdict[cName], z)
        return rdict


class DiscreteEuler(nSLDObj):
    """
        Uses precalculated Euler rotation densities
    """
    def __init__(self, generic_filename, betas, gammas, betafmt='%i', gammafmt='%i', **kwargs):
        """
            requires precalculated 4-column files for each angle beta and gamma. These should not have any smoothing
            applied. Each file must have same z vector and one header row. z can be negative.
            generic_filename contains the path and a generic file name with positions for angles beta and gamma are
            marked with <beta> and <gamma>. Format strings other than integer angles can be specified with betafmt and
            gammafmt keyword arguments. "betas" and "gammas" are lists or numpy arrays containing the beta and gamma
            points to be loaded.
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
                    zdata = d[:, 0]
                    areadata = numpy.zeros((len(betas), len(gammas), len(zdata)))
                    nslHdata = numpy.zeros_like(areadata)
                    nslDdata = numpy.zeros_like(areadata)
                areadata[i, j, :] = d[:, 1]
                nslHdata[i, j, :] = d[:, 2]
                nslDdata[i, j, :] = d[:, 3]

        self.betas = betas
        self.gammas = gammas
        self.areadata = areadata
        # TODO: self.zdata and self.zaxis from the parent might be redundant. Needs consolidation.
        self.zdata = zdata
        self.nslHdata = nslHdata
        self.nslDdata = nslDdata
        self.beta = 0.  # current beta rotation
        self.gamma = 0.  # current gamma rotation
        self.sigma = 2.  # smoothing function
        self.z = 0.  # z offset value
        self.nf = 1.  # number fraction
        self.fnSetBulknSLD(None)
        self.protexchratio = 1.  # proton exchange ratio

        # TODO: Would it make sense to have a concept of "normarea"? Then there could be a "volume fraction"
        # concept so that max(area) = volume_fraction * normarea

    def fnGetProfiles(self, z):
        # perform interpolation
        self._set_interpolation_points(z)
        area = interpn((self.betas, self.gammas, self.zdata + self.z), self.areadata, self._interppoint,
                       bounds_error=False, fill_value=0.0)
        nslH = interpn((self.betas, self.gammas, self.zdata + self.z), self.nslHdata, self._interppoint,
                       bounds_error=False, fill_value=0.0)
        nslD = interpn((self.betas, self.gammas, self.zdata + self.z), self.nslDdata, self._interppoint,
                       bounds_error=False, fill_value=0.0)

        # apply roughness (sigma)
        dz = z[1] - z[0]
        area = gaussian_filter(area, self.sigma / dz, order=0, mode='constant', cval=0)
        # area /= dz

        nslH = gaussian_filter(nslH, self.sigma / dz, order=0, mode='constant', cval=0)

        if self.bulknsld is not None:
            # get nslD profile
            nslD = gaussian_filter(nslD, self.sigma / dz, order=0, mode='constant', cval=0)
            fracD = self.protexchratio * (self.bulknsld - H2O_SLD) / (D2O_SLD - H2O_SLD)
            nsl = fracD * nslD + (1 - fracD) * nslH

        else:
            nsl = nslH

        self.zaxis = z
        self.area = area * self.nf
        self.sl = nsl * self.nf
        self.sld = numpy.divide(nsl, area * numpy.gradient(z), out=numpy.zeros_like(area), where=area!=0)

        return self.area, self.sl, self.sld

    def _set_interpolation_points(self, z):
        self._interppoint = numpy.zeros((len(z), 3))
        self._interppoint[:, :-1] = [self.beta, self.gamma]
        self._interppoint[:, -1] = z

    def fnGetVolume(self, z1=None, z2=None, recalculate=True):
        """
            Calculates volume based on the number of residues of the rotated molecule located between
            z positions z1 and z2 (inclusive).
            
            Note: the result is (slightly) different from integrating the area, because the roughness has already
            been applied to the area. However, this remains more accurate than the roughened value because
            the integration limits  will typically correspond to the limits of a lipid Box2Err function
            which are themselves defined before the roughness is applied.
        """
        if z1 is None or z2 is None or recalculate is False:
            return super().fnGetVolume(z1, z2, recalculate)
        else:
            # For volume calculation, use intrinsic z vector, and no smoothing.
            zvol = self.zdata + self.z
            self._set_interpolation_points(zvol)

            area = interpn((self.betas, self.gammas, zvol), self.areadata, self._interppoint, bounds_error=False,
                           fill_value=0.0)

            crit = (zvol > min(z1, z2)) & (zvol < max(z1, z2))
            if numpy.any(crit):
                volume = numpy.trapz(area[crit], zvol[crit])
            else:
                volume = 0.0

            return volume * self.nf

    def fnSet(self, beta, gamma, zpos, sigma, nf, bulknsld=None):
        self.beta = beta
        self.gamma = gamma
        self.z = zpos
        self.nf = nf
        self.sigma = sigma
        self.fnSetBulknSLD(bulknsld)

    def fnWriteGroup2File(self, fp, cName, z):
        fp.write(f"DiscreteEuler {cName} StartPosition {self.z} Gamma {self.gamma} Beta {self.beta} nf {self.nf} \n")
        self.fnWriteProfile2File(fp, cName, z)

    def fnWriteGroup2Dict(self, rdict, cName, z):
        rdict[cName] = {}
        rdict[cName]['header'] = f"DiscreteEuler {cName} StartPosition {self.z} Gamma {self.gamma} Beta {self.beta} nf {self.nf}"
        rdict[cName]['startposition'] = self.z
        rdict[cName]['beta'] = self.beta
        rdict[cName]['gamma'] = self.gamma
        rdict[cName]['nf'] = self.nf
        rdict[cName] = self.fnWriteProfile2Dict(rdict[cName], z)
        return rdict


class ContinuousEulerMissingResidues(CompositenSLDObj):
    """Composite protein model with arbitrary number of missing residue loops attached
        to specific residues in a ContinuousEuler model
        Requires a Continuous Euler model.

        Missing residues are specified using the add_missing_residues method; the
        resulting object is returned as a TetheredBoxDouble object
        
        """

    def __init__(self, euler: ContinuousEuler, **kwargs):
        super().__init__([euler], **kwargs)
        self.euler = euler
        self.missing_residues = []

    def fnFindSubgroups(self):
        # override default behavior so don't have to store objects separately
        pass

    def add_missing_residues(self, sequence, attachment_residues, deut_res=[], name=None):
        """
        Add missing residues.

        Usage: missing_residues_box = add_missing_residues(...)

        Requires:
        sequence -- string of single-residue codes
        attachment residues -- int or (int, int): residue numbers where loops or disordered regions
            are attached
        deut_res -- None if no deuterated residues; else indices of residues in "sequence"
            that are deuterated

        Returns:
        new_obj -- new TetheredBoxDouble object. Use fnSet to set length and frac_position
        """

        default_name = f'missing{len(self.missing_residues)}'

        name = name if name is not None else default_name

        elements = default_table()

        # Analyze sequence
        volume = 0.0
        nslH = 0.0
        nslD = 0.0
        for i, iaa in enumerate(sequence):
            if i in deut_res:
                resmol = Molecule(name='Dres', formula=aa[iaa].formula.replace(elements.H, elements.D),
                                  cell_volume=aa[iaa].cell_volume)
            else:
                resmol = aa[iaa]

            volume += resmol.cell_volume
            nslH += resmol.cell_volume * resmol.sld * 1e-6
            nslD += resmol.cell_volume * resmol.Dsld * 1e-6

        if isinstance(attachment_residues, float):
            attachment_residues = int(attachment_residues)

        if isinstance(attachment_residues, int):
            attachment_residues = (attachment_residues, attachment_residues)

        new_obj = TetheredBoxDouble()

        self.missing_residues.append({'sequence': sequence,
                                      'volume': volume,
                                      'nslH': nslH,
                                      'nslD': nslD,
                                      'attach': attachment_residues,
                                      'object': new_obj})
        
        self.subgroups.append(new_obj)

        return new_obj

    def _update_missing_residues(self):
        """Updates all missing residues with current Euler coordinates and
            number fraction"""
        resnums = list(self.euler.resnums)

        for mr in self.missing_residues:
            box: TetheredBoxDouble = mr['object']
            tp1 = self.euler.rotcoords[resnums.index(mr['attach'][0]),2]
            tp2 = self.euler.rotcoords[resnums.index(mr['attach'][1]),2]
            
            # set Euler-specific positions. Other quantities (length, frac_position)
            # must be set separately
            box.fnSet(volume=mr['volume'],
                      tether_position1=tp1,
                      tether_position2=tp2,
                      nf=self.euler.nf)
            box.fnSetnSL(mr['nslH'], mr['nslD'])
            box.fnSetSigma(self.euler.sigma)

    def fnSet(self, gamma, beta, zpos, sigma, nf, bulknsld=None):
        self.euler.fnSet(gamma, beta, zpos, sigma, nf, bulknsld)
        self.fnSetBulknSLD(bulknsld)
        self._update_missing_residues()


# ------------------------------------------------------------------------------------------------------
# Bilayer + protein models
# ------------------------------------------------------------------------------------------------------

class BLMProteinComplex(CompositenSLDObj):
    """
    Composite bilayer/protein model

    Requires:
    blm: list of at least one bilayer object (BLM or its subclasses)
    proteins: list of any number of protein or protein-like objects (i.e. has fnGetVolume),
                e.g. Hermite, SLDHermite, DiscreteEuler, ContinuousEuler

    Note: Bilayers should be non-overlapping for hc substitution to work. This is currently not enforced or tested for.
    
    Example usage:
    blm = BLM(...)
    prot1 = Hermite(...)
    prot2 = Hermite(...)
    blmprot = BLMProteinComplex(blms=[blm], proteins=[prot1, prot2])
    
    Setting parameters should occur as usual on subgroups:
    blmprot.blms[0].fnSet(...) or blmprot.blms[0].fnSet(...)
    blmprot.proteins[0].fnSet(...)
    blmprot.proteins[1].fnSet(...)

    Adjusts bilayer when protein is present:
    blmprot.fnAdjustBLMs() # adjusts bilayers for presence of protein
    """
    def __init__(self, blms=None, proteins=None, name=None):
        super().__init__(name=name)

        assert len(blms) > 0, 'At least one bilayer is required in BLMProteinComplex'
        if blms is not None:
            self.blms = CompositenSLDObj(blms, name='blms')
        if proteins is not None:
            self.proteins = CompositenSLDObj(proteins, name='proteins')

        self.nf = 1.0
        self.normarea = 0.0
        self.fnFindSubgroups()

    def fnAdjustBLMs(self):
        """ Adjust bilayers for the presence of proteins.
        """

        dMaxArea = 0
        for blm in self.blms:
            # intial bilayer adjustment with respect to the current parameters and withouth hc substitution
            blm.hc_substitution_1 = 0
            blm.hc_substitution_2 = 0
            blm.fnAdjustParameters()
            normarea = blm.normarea
            dMaxArea = max(normarea, dMaxArea)

        self.normarea = dMaxArea

        for prot in self.proteins:
            if hasattr(prot, 'normarea'):
                prot.fnSetNormarea(dMaxArea)

        for blm in self.blms:
            # inner leaflet
            z1 = blm.methylenes1[0].z - 0.5 * blm.methylenes1[0].length
            z2 = blm.methyls1[0].z + 0.5 * blm.methyls1[0].length
            lipidvol = sum(methylene.vol for methylene in blm.methylenes1) + sum(methyl.vol for methyl in blm.methyls1)
            lipidvol *= blm.vf_bilayer
            # TODO: Create numerical volume function on the level of nSLDObj, which might be faster than the spline
            #  integration and more flexible
            blm.hc_substitution_1 = self.proteins.fnGetVolume(z1, z2, True) / lipidvol

            # outer leaflet
            lipidvol = sum(methylene.vol for methylene in blm.methylenes2) + sum(methyl.vol for methyl in blm.methyls2)
            lipidvol *= blm.vf_bilayer
            z1 = blm.methyls2[0].z - 0.5 * blm.methyls2[0].length
            z2 = blm.methylenes2[0].z + 0.5 * blm.methylenes2[0].length
            blm.hc_substitution_2 = self.proteins.fnGetVolume(z1, z2, True) / lipidvol

            # final adjustement of the bilayer after determining the hc substitution values.
            blm.fnAdjustParameters()

    def fnGetProfiles(self, z):
        # calculate total area from bilayers
        blm_area, blm_nsl, _ = self.blms.fnGetProfiles(z)
        dMaxArea = blm_area.max()
        # calculate total area from proteins and overlay proteins on bilayers
        area, nsl = self.proteins.fnOverlayProfile(z, blm_area, blm_nsl, dMaxArea)
        nsld = numpy.divide(nsl, area * numpy.gradient(z), out=numpy.zeros_like(area), where=area > 0)

        self.zaxis = z
        self.area = area * self.nf
        self.sl = nsl * self.nf
        self.sld = nsld

        return self.area, self.sl, self.sld

    def fnWriteResults2Dict(self, rdict, cName):
        rdict = super().fnWriteResults2Dict(rdict, cName)
        if cName not in rdict:
            rdict[cName] = {}

        vprot = self.proteins.fnGetVolume(recalculate=False)

        if self.normarea != 0 and vprot != 0:
            p1 = self.blms[0].substrate.z + 0.5 * self.blms[0].substrate.length
            p2 = self.blms[0].headgroups1[0].z - 0.5 * self.blms[0].headgroups1[0].length
            p3 = self.blms[0].methylenes1[0].z - 0.5 * self.blms[0].methylenes1[0].length
            p4 = self.blms[0].methyls1[0].z + 0.5 * self.blms[0].methyls1[0].length
            p5 = self.blms[0].methylenes2[0].z + 0.5 * self.blms[0].methylenes2[0].length
            p6 = self.blms[0].headgroups2[0].z + 0.5 * self.blms[0].headgroups2[0].length

            rdict[cName]['protein in submembrane'] = self.proteins.fnGetVolume(p1, p2, recalculate=False) / vprot
            rdict[cName]['protein in inner headgroup'] = self.proteins.fnGetVolume(p2, p3, recalculate=False) / vprot
            rdict[cName]['protein in inner hydrocarbon'] = self.proteins.fnGetVolume(p3, p4, recalculate=False) / vprot
            rdict[cName]['protein in outer hydrocarbon'] = self.proteins.fnGetVolume(p4, p5, recalculate=False) / vprot
            rdict[cName]['protein in outer headgroup'] = self.proteins.fnGetVolume(p5, p6, recalculate=False) / vprot
            rdict[cName]['protein in headgroups'] = rdict[cName]['protein in inner headgroup']
            rdict[cName]['protein in headgroups'] += rdict[cName]['protein in outer headgroup']
            rdict[cName]['protein in hydrocarbons'] = rdict[cName]['protein in inner hydrocarbon']
            rdict[cName]['protein in hydrocarbons'] += rdict[cName]['protein in outer hydrocarbon']
            rdict[cName]['protein in bulk'] = 1 - rdict[cName]['protein in submembrane']
            rdict[cName]['protein in bulk'] -= rdict[cName]['protein in headgroups']
            rdict[cName]['protein in bulk'] -= rdict[cName]['protein in hydrocarbons']

            rdict[cName]['water in submembrane'] = self.fnGetVolume(p1, p2, recalculate=False)
            rdict[cName]['water in submembrane'] /= (self.normarea * (p2 - p1))
            rdict[cName]['water in submembrane'] = 1 - rdict[cName]['water in submembrane']
            rdict[cName]['water in inner headgroups'] = self.fnGetVolume(p2, p3, recalculate=False)
            rdict[cName]['water in inner headgroups'] /= (self.normarea * (p3 - p2))
            rdict[cName]['water in inner headgroups'] = 1 - rdict[cName]['water in inner headgroups']
            rdict[cName]['water in hydrocarbons'] = self.fnGetVolume(p3, p5, recalculate=False)
            rdict[cName]['water in hydrocarbons'] /= (self.normarea * (p5 - p3))
            rdict[cName]['water in hydrocarbons'] = 1 - rdict[cName]['water in hydrocarbons']
            rdict[cName]['water in outer headgroups'] = self.fnGetVolume(p5, p6, recalculate=False)
            rdict[cName]['water in outer headgroups'] /= (self.normarea * (p6 - p5))
            rdict[cName]['water in outer headgroups'] = 1 - rdict[cName]['water in outer headgroups']

            rdict[cName]['total protein [cvo]'] = vprot / self.normarea

            pos = numpy.argmax(self.proteins.area)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                results = peak_widths(self.area, [pos], rel_height=0.5)
            rdict[cName]['total protein peak from hg'] = self.zaxis[pos] - p6
            rdict[cName]['total protein FWHM'] = results[0] * (self.zaxis[1] - self.zaxis[0])

        return rdict

# ------------------------------------------------------------------------------------------------------
# Polymer models
# ------------------------------------------------------------------------------------------------------

class PolymerMushroom(nSLDObj):
    """Polymer mushroom model (relatively low grafting density)

        Uses refl1d.polymer.mushroom_math
    """
    def __init__(self, startz=0.0, rho=0.7e-6, vf=0.1, Rg=7.0, delta=1.0, normarea=1.0, sigma=2.0, name=None):
        super().__init__(name=name)
        self.startz = startz
        self.rho = rho
        self.vf = vf
        self.Rg = Rg
        self.delta = delta
        self.normarea = normarea
        self.sigma = sigma
        self.nf = 1.0

    def fnGetProfiles(self, z):

        v = numpy.zeros_like(z)
        x = (z - self.startz) / self.Rg

        # protect against divide by zero error (from refl1d.polymer)
        delta_thresh = 1e-10

        if abs(self.delta) > delta_thresh:
            v[x > 0] = mushroom_math(x[x>0], self.delta, self.vf)
        else: # we should RARELY get here
            scale = (self.delta+delta_thresh)/2.0/delta_thresh
            v[x > 0] = (scale*mushroom_math(x[x>0], delta_thresh, self.vf)
                                + (1.0-scale)*mushroom_math(x[x>0], -delta_thresh, self.vf))

        # convolve with roughness        
        self.area = self.normarea * smear(z, v, self.sigma) * self.nf
        self.zaxis = z
        self.sl = self.area * self.rho
        self.sld = numpy.ones_like(self.area) * self.rho

        return self.area, self.sl, self.sld
    
    def fnWriteResults2Dict(self, rdict, cName):
        rdict = super().fnWriteResults2Dict(rdict, cName)
        if cName not in rdict:
            rdict[cName] = {}
        vol = trapezoid(self.area, self.zaxis)
        rdict[cName]['INT'] = vol
        rdict[cName]['COM'] = trapezoid(self.area * self.zaxis, self.zaxis) / vol if vol != 0 else 0
        rdict[cName]['max area'] = max(self.area)
        maxpos = numpy.argmax(self.area)
        rdict[cName]['starting position'] = self.startz
        rdict[cName]['position of maximum area'] = self.zaxis[maxpos]
        rdict[cName]['position of half-height'] = numpy.interp(0.5 * max(self.area), self.area[maxpos:][::-1], self.zaxis[maxpos:][::-1])
        rdict[cName]['distance to maximum area'] = rdict[cName]['position of maximum area'] - self.startz
        rdict[cName]['distance to half-height'] = rdict[cName]['position of half-height'] - self.startz
        return rdict

class PolymerBrush(nSLDObj):
    """Polymer brush model (parabolic profile)
    
        Adapted from refl1d.polymer.PolymerBrush.profile
    """
    def __init__(self, startz=0.0, base_length=10, interface_length=10, thinning_power=1, rho=0.7e-6, vf=0.1, normarea=1.0, sigma=2.0, name=None):
        super().__init__(name=name)
        self.startz = startz
        self.rho = rho
        self.vf = vf
        self.base_length = base_length
        self.interface_length = interface_length
        self.thinning_power = thinning_power
        self.normarea = normarea
        self.sigma = sigma
        self.nf = 1.0

    def fnGetProfiles(self, z):

        L0 = self.startz + self.base_length
        L1 = L0 + self.interface_length
        if self.interface_length == 0:
            v = numpy.ones_like(z)
        else:
            v = (1 - ((z-L0)/(L1-L0))**2)
        v[z < L0] = 1
        v[z > L1] = 0
        brush_profile = self.vf * v**self.thinning_power
        brush_profile[z < self.startz] = 0

        # convolve with roughness
        vf = smear(z, brush_profile, self.sigma)
        self.area = self.normarea * vf * self.nf
        self.zaxis = z
        self.sl = self.area * self.rho
        self.sld = numpy.ones_like(self.area) * self.rho

        return self.area, self.sl, self.sld
    
    def fnWriteResults2Dict(self, rdict, cName):
        rdict = super().fnWriteResults2Dict(rdict, cName)
        if cName not in rdict:
            rdict[cName] = {}
        vol = trapezoid(self.area, self.zaxis)
        rdict[cName]['INT'] = vol
        rdict[cName]['COM'] = trapezoid(self.area * self.zaxis, self.zaxis) / vol if vol != 0 else 0
        rdict[cName]['max area'] = max(self.area)
        maxpos = numpy.argmax(self.area)
        rdict[cName]['starting position'] = self.startz
        rdict[cName]['position of half-height'] = numpy.interp(0.5 * max(self.area), self.area[maxpos:][::-1], self.zaxis[maxpos:][::-1])
        rdict[cName]['distance to half-height'] = rdict[cName]['position of half-height'] - self.startz
        return rdict
