import libadjoint
import numpy

class Vector(libadjoint.Vector):
    def __init__(self, vec):
        '''Vector(vec)
        class to wrap the numpy vector vec in a libadjoint vector.'''

        self.vec=vec

    def duplicate(self):
        
        return Vector(numpy.zeros(self.vec.size))

    
    
