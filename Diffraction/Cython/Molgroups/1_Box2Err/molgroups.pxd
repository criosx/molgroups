cdef extern from "molgroups.cc":
    pass

# Declare the class with cdef
cdef extern from "molgroups.h":
    cdef cppclass Box2Err:
        Box2Err() except +
        Box2Err(double, double, double, double, double, double) except +
        double fnGetLowerLimit()
