import libadjoint
import numpy

class Vector(libadjoint.Vector):
  def __init__(self, vec):
    '''Vector(vec)
    class to wrap the numpy vector vec in a libadjoint vector.'''
    self.vec = vec

  def duplicate(self):
    return Vector(numpy.zeros(self.vec.size))

  def axpy(self, alpha, x):
    self.vec += alpha*x.vec

class Matrix(libadjoint.Matrix):
  def __init__(self, mat):
    '''Matrix(mat)
    class to wrap the numpy matrix mat in a libadjoint matrix.'''
    self.mat = mat

  def duplicate(self):
    return Matrix(numpy.zeros(self.mat.shape))

  def axpy(self, alpha, x):
    self.mat += alpha*x.mat

