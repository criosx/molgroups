cdef extern from "rectangle.cc":
    pass

# Declare the class with cdef
cdef extern from "rectangle.h":
    cdef cppclass Square:
        Square() except +
        Square(int, int) except +
        int x0, x1
        # int x0, y0, x1, y1
        int getArea()
        void getSize(int* width, int* height)
        void move(int, int)
