from __future__ import absolute_import
from . import libadjoint
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

  def set_values(self, scalars):
    self.vec[:] = scalars

  def size(self):
    return self.vec.size
      
  def norm(self):
    return numpy.linalg.norm(self.vec)

  def set_random(self):
    self.vec = numpy.random.rand(self.get_size())

  def dot_product(self,b):
    return numpy.dot(self.vec, b.vec)

class Matrix(libadjoint.Matrix):
  def __init__(self, mat):
    '''Matrix(mat)
    class to wrap the numpy matrix mat in a libadjoint matrix.'''
    self.mat = mat

  def solve(self, var, bb):


    return Vector(numpy.linalg.solve(self.mat, bb.vec))
  

  def axpy(self, alpha, x):
    self.mat += alpha*x.mat

class RHS(libadjoint.RHS):
  def __init__(self, v):
    self.v = v
  def __call__(self, dependencies, values):
    assert len(dependencies) == 0
    return self.v
  def dependencies(self):
    return []

