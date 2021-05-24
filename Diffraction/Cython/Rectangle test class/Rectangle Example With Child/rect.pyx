# distutils: language = c++

from rectangle cimport Square

cdef class PySquare:
    cdef Square c_square

    def __cinit__(self, int x0, int x1):
        self.c_square = Square(x0, x1)

    def get_area(self):
        return self.c_square.getArea()

    def get_size(self):
        cdef int width, height
        self.c_square.getSize(&width, &height)
        return width, height

    def move(self, dx, dy):
        self.c_square.move(dx, dy)

    # Attribute access
    @property
    def x0(self):
        return self.c_square.x0
    @x0.setter
    def x0(self, x0):
        self.c_square.x0 = x0

    # Attribute access
    @property
    def x1(self):
        return self.c_square.x1
    @x1.setter
    def x1(self, x1):
        self.c_square.x1 = x1

#    # Attribute access
#    @property
#    def y0(self):
#        return self.c_square.y0
#    @y0.setter
#    def y0(self, y0):
#        self.c_square.y0 = y0

#    # Attribute access
#    @property
#    def y1(self):
#        return self.c_square.y1
#    @y1.setter
#    def y1(self, y1):
#        self.c_square.y1 = y1
