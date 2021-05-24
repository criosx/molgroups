# distutils: language = c++

from molgroups cimport Box2Err

cdef class PyBox2Err:
    cdef Box2Err c_box2err

    def __cinit__(self, double d0, double d1, double d2, double d3, double d4, double d5):
        self.c_box2err = Box2Err(d0, d1, d2, d3, d4, d5)

    def get_lowerlimit(self):
        return self.c_box2err.fnGetLowerLimit()