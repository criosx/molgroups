# distutils: language = c++

from molgroups cimport BLM_quaternary
from cpython cimport array
import array


cdef class PyBLM_quaternary:
    cdef BLM_quaternary c_BLM_quaternary
    def AdjustParameters(self):
        return self.c_BLM_quaternary.fnAdjustParameters()
    def GetLowerLimit(self):
        return self.c_BLM_quaternary.fnGetLowerLimit()
    def GetUpperLimit(self):
        return self.c_BLM_quaternary.fnGetUpperLimit()
    def GetArea(self, double z):
        return self.c_BLM_quaternary.fnGetArea(z)
    def GetCenter(self):
        return self.c_BLM_quaternary.fnGetCenter()
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
    def Set(self, double sigma, double bulknsld, double startz,
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
        cdef double dd = self.c_BLM_quaternary.fnWriteProfile(&Area_mv[0], &nSLD_mv[0], dimension, stepsize, dMaxArea)
        aArea = list(Area_mv)
        anSLD = list(nSLD_mv)
        return dd, aArea, anSLD
    def OverlayProfile(self, list aArea, list anSLD, int dimension,
                       double stepsize,double dMaxArea):
        cdef array.array Area = array.array('d', aArea)
        cdef array.array nSLD = array.array('d', anSLD)
        cdef double[::1] Area_mv = Area
        cdef double[::1] nSLD_mv = nSLD
        return self.c_BLM_quaternary.fnOverlayProfile(&Area_mv[0], &nSLD_mv[0],
                                                    dimension, stepsize,
                                                    dMaxArea)
