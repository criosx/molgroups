# distutils: language = c++

from molgroups cimport Box2Err, Hermite, ssBLM_POPC_POPS_PIP_chol, BLM_quaternary
from cpython cimport array
import array

cdef class PyBox2Err:
    cdef Box2Err c_box2err

    def __cinit__(self, double d0, double d1, double d2, double d3, double d4, double d5):
        self.c_box2err = Box2Err(d0, d1, d2, d3, d4, d5)
    def get_lowerlimit(self):
        return self.c_box2err.fnGetLowerLimit()

cdef class PyHermite:
    # Hermite does not have a nullary constructor and therefore needs to be 
    # initialized dynamically using new and del (see below)
    cdef Hermite *c_hermite

    def __cinit__(self):
        self.c_hermite = new Hermite()
    def __cinit__(self, int n, double dstartposition, double dnSLD, double dnormarea):
        self.c_hermite = new Hermite(n, dstartposition, dnSLD, dnormarea)
    def __dealloc__(self):
        del self.c_hermite
    def GetArea(self, double dz):
        return self.c_hermite.fnGetArea(dz)
    def GetnSLD(self, double dz):
        return self.c_hermite.fnGetnSLD(dz)
    def GetLowerLimit(self):
        return self.c_hermite.fnGetLowerLimit()
    def GetUpperLimit(self):
        return self.c_hermite.fnGetUpperLimit()
    def GetVolume(self, double dz1, double dz2):
        return self.c_hermite.fnGetVolume(dz1,dz2)
    def SetNormarea(self, double dnormarea):
        return self.c_hermite.fnSetNormarea(dnormarea)
    def SetnSLD(self, double dnSLD):
        return self.c_hermite.fnSetnSLD(dnSLD)
    def SetRelative(self, double dSpacing, double dStart, list dDp, list dVf, double dnf):
        # the python side expects a list, which is converted into an array and
        # then into a memoryview object, whose first element's address is passed on
        # also see the imports at the beginning of the file
        # might be simplified, ToDo
        cdef array.array Dp = array.array('d', dDp)
        cdef array.array Vf = array.array('d', dVf)
        cdef double[::1] Dp_mv = Dp
        cdef double[::1] Vf_mv = Vf
        return self.c_hermite.fnSetRelative(dSpacing, dStart, &Dp_mv[0], &Vf_mv[0], dnf)
    def SetSigma(self, double sigma):
        return self.c_hermite.fnSetSigma(sigma)

cdef class PyBLM_quaternary:
    cdef BLM_quaternary *c_BLM_quaternary

    def __cinit__(self):
        self.c_BLM_quaternary = new BLM_quaternary()
    def __dealloc__(self):
        del self.c_BLM_quaternary
    def AdjustParameters(self):
        return self.c_BLM_quaternary.fnAdjustParameters()
    def GetLowerLimit(self):
        return self.c_BLM_quaternary.fnGetLowerLimit()
    def GetUpperLimit(self):
        return self.c_BLM_quaternary.fnGetUpperLimit()
    def GetArea(self, double z):
        return self.c_BLM_quaternary.fnGetArea(z)
    def GetnSLD(self, double z):
        return self.c_BLM_quaternary.fnGetnSLD(z)
    def Init(self, double va1, double na1, double vm1, double nm1, double vh1,
             double nh1, double lh1, double va2, double na2, double vm2,
             double nm2, double vh2, double nh2, double lh2, double va3,
             double na3, double vm3, double nm3, double vh3, double nh3,
             double lh3, double vc, double nc):
        return self.c_BLM_quaternary.fnInit(va1, na1, vm1, nm1, vh1, nh1,
                                            lh1, va2, na2, vm2, nm2, vh2,
                                            nh2, lh2, va3, na3, vm3, nm3,
                                            vh3, nh3, lh3, vc, nc)
    def fnSet(self, double sigma, double bulknsld, double startz,
                double l_lipid1, double l_lipid2, double vf_bilayer,
                double nf_lipid_2=0,double nf_lipid3=0,double nf_chol=0,
                double hc_substitution_1=0, double hc_substitution_2=0,
                double radius_defect=100):
        return self.c_BLM_quaternary.fnSet(sigma, bulknsld, startz, l_lipid1,
                    l_lipid2, vf_bilayer, nf_lipid_2,
                    nf_lipid3,nf_chol, hc_substitution_1,
                    hc_substitution_2, radius_defect)
    def SetSigma(self, double sigma):
        return self.c_BLM_quaternary.fnSetSigma(sigma)
    def WriteProfile(self, list aArea, list anSLD,
                        int dimension, double stepsize, double dMaxArea):
        cdef array.array Area = array.array('d', aArea)
        cdef array.array nSLD = array.array('d', anSLD)
        cdef double[::1] Area_mv = Area
        cdef double[::1] nSLD_mv = nSLD
        return self.c_BLM_quaternary.fnWriteProfile(&Area_mv[0], &nSLD_mv[0],
                                                    dimension, stepsize,
                                                    dMaxArea)
    def OverlayProfile(self, list aArea, list anSLD, int dimension,
                       double stepsize,double dMaxArea):
        cdef array.array Area = array.array('d', aArea)
        cdef array.array nSLD = array.array('d', anSLD)
        cdef double[::1] Area_mv = Area
        cdef double[::1] nSLD_mv = nSLD
        return self.c_BLM_quaternary.fnOverlayProfile(&Area_mv[0], &nSLD_mv[0],
                                                    dimension, stepsize,
                                                    dMaxArea)

cdef class Py_ssBLM_POPC_POPS_PIP_chol:
    cdef ssBLM_POPC_POPS_PIP_chol c_ssBLM_POPC_POPS_PIP_chol

   #  def __cinit__(self):
   #     self.c_ssBLM_POPC_POPS_PIP_chol = ssBLM_POPC_POPS_PIP_chol()
    def AdjustParameters(self):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnAdjustParameters()
    def GetLowerLimit(self):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnGetLowerLimit()
    def GetUpperLimit(self):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnGetUpperLimit()
    def GetArea(self, double z):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnGetArea(z)
    def GetnSLD(self, double z):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnGetnSLD(z)
    def Set(self, double sigma, double global_rough, double rho_substrate,
            double bulknsld, double rho_siox, double l_siox,
            double l_submembrane, double l_lipid1, double l_lipid2,
            double vf_bilayer, double nf_lipid_2=0, double nf_lipid3=0,
            double nf_chol=0, double hc_substitution_1=0,
            double hc_substitution_2=0, double radius_defect=100):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnSet(sigma, global_rough, rho_substrate,
        bulknsld, rho_siox, l_siox, l_submembrane, l_lipid1, l_lipid2, vf_bilayer,
        nf_lipid_2, nf_lipid3, nf_chol, hc_substitution_1, hc_substitution_2, radius_defect)
    def SetSigma(self, double sigma):
        return self.c_ssBLM_POPC_POPS_PIP_chol.fnSetSigma(sigma)
