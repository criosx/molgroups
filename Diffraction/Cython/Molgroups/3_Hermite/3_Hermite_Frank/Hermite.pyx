from molgroups cimport Hermite
from cython.operator cimport dereference as deref
from cython.view cimport array as cvarray


cdef class PyHermite:
    cdef Hermite c_hermite

    def Hermite(self, int n, double dstartposition, double dnSLD, double dnormarea):
        self.c_hermite = self.Hermite(n, dstartposition, dnSLD, dnormarea)
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
    def SetRelative(self, double dSpacing, double dStart, double[:] dDp, double[:] dVf, double dnf):
        return self.c_hermite.fnSetRelative(dSpacing, dStart, &dDp[0], &dVf[0], dnf)
    def SetSigma(self, double sigma):
        return self.c_hermite.fnSetSigma(sigma)

    def GetSplineAntiDerivative(self, double dz, double[:] dp, double[:] dh):
        return self.c_hermite.fnGetSplineAntiDerivative(dz, &dp[0], &dh[0])
    def GetSplineArea(self, double dz, double[:] dp, double[:] dh, int damping):
        return self.c_hermite.fnGetSplineArea(dz, &dp[0], &dh[0], damping)
    def GetSplinePars(self, double d, double[:] dp, double[:] dh, double &m0, double &m1, double &p0, double &p1):
        return self.c_hermite.fnGetSplinePars(d, &dp[0], &dh[0], m0, m1, p0, p1)
    def GetSplineIntegral(self, double dz1, double dz2, double[:] dp, double[:] dh, int damping):
        return self.c_hermite.fnGetSplineIntegral(dz1, dz2, &dp[0], &dh[0], damping)
    def GetSplineProductIntegral(self, double dz1, double dz2, double[:] dp, double[:] dh1, double[:] dh2, int damping1, int damping2):
        return self.c_hermite.fnGetSplineProductIntegral(dz1, dz2, &dp[0], &dh1[0], &dh2[0], damping1, damping2)