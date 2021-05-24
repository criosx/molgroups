cdef extern from "molgroups.cc":
    pass

# Declare the class with cdef
cdef extern from "molgroups.h":

    cdef cppclass BLM_quaternary:
        # BLM_quaternary() except +
        void   fnAdjustParameters()
        double fnGetLowerLimit()
        double fnGetUpperLimit()
        double fnGetArea(double)
        double fnGetCenter()
        double fnGetnSLD(double)
        void   fnInit(double, double, double double, double, double, double,
                    double, double, double, double, double, double, double,
                    double, double, double, double, double, double, double,
                    double, double, double)
        void   fnSet(double, double, double, double,double, double, double,
                     double, double, double, double, double)
        void   fnSetSigma(double sigma)
        double fnWriteProfile(double[], double[], int, double, double)
        void   fnOverlayProfile(double[], double[], int, double, double)
