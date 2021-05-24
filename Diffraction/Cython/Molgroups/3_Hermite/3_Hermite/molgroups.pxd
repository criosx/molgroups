cdef extern from "molgroups.cc":
    pass

# Declare the class with cdef

# Declare the class with cdef
cdef extern from "molgroups.h":

    # cdef cppclass Box2Err:
    #     Box2Err() except +
    #     Box2Err(double, double, double, double, double, double) except +
    #     double fnGetLowerLimit()

    # cdef cppclass ssBLM_POPC_POPS_PIP_chol:
    #     ssBLM_POPC_POPS_PIP_chol() except +
    #     void   fnAdjustParameters()
    #     double fnGetLowerLimit()
    #     double fnGetUpperLimit()
    #     double fnGetArea(double)
    #     double fnGetnSLD(double)
    #     void   fnSet(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);
    #     void   fnSetSigma(double)

    cdef cppclass Hermite:
        Hermite() except +
        Hermite(int, double, double, double) except +
        double fnGetArea(double)
        double fnGetnSLD(double)
        double fnGetLowerLimit()
        double fnGetUpperLimit()
        double fnGetVolume(double, double)
        void fnSetNormarea(double)
        void fnSetnSLD(double)
        void fnSetRelative(double, double, double[], double[], double)
        void fnSetSigma(double)
        
    # protected:
        
        double fnGetSplineAntiDerivative(double, double[], double[])
        double fnGetSplineArea(double, double[], double[], int)
        int    fnGetSplinePars(double, double[], double[], double &m0, double &m1, double &p0, double &p1)
        double fnGetSplineIntegral(double, double, double[], double[], int)
        double fnGetSplineProductIntegral(double, double, double[], double[], double[], int, int)