cdef extern from "molgroups.cc":
    pass

# Declare the class with cdef

# Declare the class with cdef
cdef extern from "molgroups.h":

    cdef cppclass Box2Err:
        Box2Err() except +
        Box2Err(double, double, double, double, double, double) except +
        double fnGetLowerLimit()

    cdef cppclass ssBLM_POPC_POPS_PIP_chol:
        # ssBLM_POPC_POPS_PIP_chol() except +
        void   fnAdjustParameters();
        double fnGetLowerLimit();
        double fnGetUpperLimit();
        double fnGetArea(double);
        double fnGetnSLD(double);
        void   fnSet(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);
        void   fnSetSigma(double);
