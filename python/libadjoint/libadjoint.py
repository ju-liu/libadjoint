import clibadjoint_constants as constants
import ctypes
import exceptions
import clibadjoint as clib
import python_utils

def handle_error(ierr):
  if ierr != 0:
    exception = exceptions.get_exception(ierr)
    errstr  = clib.adj_error_msg.value
    raise exception, errstr

def list_to_carray(vars, klass):
  listtype = klass * len(vars)
  listt = listtype()
  for i in range(len(vars)):
    listt[i] = vars[i].c_object
  return listt

# Let's make handle_error the restype for all of our libadjoint functions
for member in dir(clib):
  # Looping over all of the objects this module offers us ...
  if member.startswith("adj_"):
    obj = getattr(clib, member)
    if hasattr(obj, "restype"):
      if obj.restype == ctypes.c_int:
        obj.restype = handle_error

class Variable(object):
  def __init__(self, name="", timestep=-1, iteration=0, auxiliary=False, var=None):
    if var is None:
      self.var = clib.adj_variable()
      assert name != ""
      assert timestep >= 0
      self.name = name
      clib.adj_create_variable(name, timestep, iteration, auxiliary, self.var)
      self.c_object = self.var
    else:
      self.var = var
      self.c_object = self.var
      self.name = var.name

  def __str__(self):
    buf = ctypes.create_string_buffer(255)
    clib.adj_variable_str(self.var, buf, 255)
    return buf.value

  def __getattr__(self, name):
    if name == "timestep":
      timestep = ctypes.c_int()
      clib.adj_variable_get_timestep(self.var, timestep)
      return timestep.value
    elif name == "iteration":
      iteration = ctypes.c_int()
      clib.adj_variable_get_iteration(self.var, iteration)
      return iteration.value
    else:
      raise AttributeError

class NonlinearBlock(object):
  def __init__(self, name, dependencies, context=None, coefficient=None):
    self.nblock = clib.adj_nonlinear_block()
    c_context = None
    if context is not None:
      c_context = ctypes.byref(context)
    self.c_object = self.nblock

    clib.adj_create_nonlinear_block(name, len(dependencies), list_to_carray(dependencies, clib.adj_variable), c_context, 1.0, self.nblock)

    if coefficient is not None:
      self.set_coefficient(coefficient)

  def __del__(self):
    clib.adj_destroy_nonlinear_block(self.nblock)

  def set_coefficient(self, c):
    clib.adj_nonlinear_block_set_coefficient(self.nblock, c)

class Block(object):
  def __init__(self, name, nblock=None, context=None, coefficient=None, hermitian=False, dependencies=None):
    self.block = clib.adj_block()
    self.c_object = self.block
    c_context = None
    if context is not None:
      c_context = ctypes.byref(context)

    if nblock is not None and dependencies is not None:
      raise LibadjointErrorInvalidInput, "Cannot have both nblock and dependencies"

    if nblock is None:
      if dependencies is not None and len(dependencies) > 0:
        c_nblock = NonlinearBlock(name, dependencies, context=context).nblock
      else:
        c_nblock = None
    else:
      c_nblock = nblock.nblock

    clib.adj_create_block(name, c_nblock, c_context, 1.0, self.block)

    if coefficient is not None:
      self.set_coefficient(coefficient)

    if hermitian:
      self.set_hermitian(hermitian)

  def __del__(self):
    clib.adj_destroy_block(self.block)

  def set_coefficient(self, c):
    clib.adj_block_set_coefficient(self.block, c)

  def set_hermitian(self, hermitian):
    clib.adj_block_set_hermitian(self.block, hermitian)

class Equation(object):
  def __init__(self, var, blocks, targets, rhs_deps=None, rhs_context=None, rhs_cb=None):
    self.equation = clib.adj_equation()
    self.c_object = self.equation

    assert len(blocks) == len(targets)

    clib.adj_create_equation(var.var, len(blocks), list_to_carray(blocks, clib.adj_block), list_to_carray(targets, clib.adj_variable), self.equation)

    if rhs_deps is not None and len(rhs_deps) > 0:
      clib.adj_equation_set_rhs_dependencies(self.equation, len(rhs_deps), list_to_carray(rhs_deps, clib.adj_variable), rhs_context)

    if rhs_cb is not None:
      self.rhs_fn = self.__cfunc_from_rhs__(rhs_cb)
      clib.adj_equation_set_rhs_callback(self.equation, self.rhs_fn)

  def __del__(self):
    clib.adj_destroy_equation(self.equation)

  def __cfunc_from_rhs__(self, rhs_cb):
    '''Given a rhs function defined using the Pythonic interface, we want to translate that into a function that
    can be called from C. This routine does exactly that.'''

    def cfunc(adjointer_c, variable_c, ndepends_c, dependencies_c, values_c, context_c, output_c, has_output_c):
      # build the Python objects from the C objects
      adjointer = Adjointer(adjointer_c)
      variable  = Variable(var=variable_c)
      dependencies = [Variable(var=dependencies_c[i]) for i in range(ndepends_c)]
      values = [vector(values_c[i]) for i in range(ndepends_c)]
      context = context_c

      # Now call the callback we've been given
      output = rhs_cb(adjointer, variable, dependencies, values, context)

      # Now cast the outputs back to C
      has_output_c.value = (output is None)
      output_c[0].ptr = python_utils.c_ptr(output)
      python_utils.incref(output)

    rhs_func_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_adjointer), clib.adj_variable, ctypes.c_int, ctypes.POINTER(clib.adj_variable), ctypes.POINTER(clib.adj_vector), ctypes.POINTER(None),
                                     ctypes.POINTER(clib.adj_vector), ctypes.POINTER(ctypes.c_int))

    return rhs_func_type(cfunc)


class Storage(object):
  '''Wrapper class for Vectors that contains additional information about how and where libadjoint 
  stores the vector values. This should never be used directly, instead use MemoryStorage or DiskStorage.'''

  def set_compare(self, tol):
    if tol > 0.0:
      clib.adj_storage_set_compare(self.storage_data, 1, tol)
    else:
      clib.adj_storage_set_compare(self.storage_data, 0, tol)

  def set_overwrite(self, flag=True):
    if flag:
      clib.adj_storage_set_overwrite(self.storage_data, 1)
    else:
      clib.adj_storage_set_overwrite(self.storage_data, 0)

  def set_checkpoint(self, flag=True):
    if flag:
      clib.adj_storage_set_checkpoint(self.storage_data, 1)
    else:
      clib.adj_storage_set_checkpoint(self.storage_data, 0)

class MemoryStorage(Storage):
  '''Wrapper class for Vectors that contains additional information for storing the vector values in memory.'''
  def __init__(self, vec):
    self.storage_data = clib.adj_storage_data()
    self.c_object = self.storage_data
    clib.adj_storage_memory_incref(vec.as_adj_vector(), self.storage_data)

    # Ensure that the storage object always holds a reference to the vec
    self.vec = vec


class DiskStorage(Storage):
  '''Wrapper class for Vectors that contains additional information for storing the vector values in memory.'''
  def __init__(self, vec):
    self.storage_data = clib.adj_storage_data()
    self.c_object = self.storage_data
    clib.adj_storage_disk_incref(vec.as_adj_vector(), self.storage_data)

    # Ensure that the storage object always holds a reference to the vec
    self.vec = vec

class Adjointer(object):
  def __init__(self, adjointer=None):
    if adjointer is None:
      self.adjointer = clib.adj_adjointer()
      clib.adj_create_adjointer(self.adjointer)
      self.c_object = self.adjointer
      self.__register_data_callbacks__()
    else:
      self.adjointer = adjointer
      self.c_object = self.adjointer

    self.functions_registered = []

    self.block_assembly_type = ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.POINTER(clib.adj_variable), ctypes.POINTER(clib.adj_vector), ctypes.c_int, ctypes.c_double, ctypes.POINTER(None),
                                                ctypes.POINTER(clib.adj_matrix), ctypes.POINTER(clib.adj_vector))

  def __del__(self):
    clib.adj_destroy_adjointer(self.adjointer)

  def __getattr__(self, name):
    if name == "equation_count":
      equation_count = ctypes.c_int()
      clib.adj_equation_count(self.adjointer, equation_count)
      return equation_count.value
    else:
      raise AttributeError

  def register_equation(self, equation):
    cs = ctypes.c_int()
    clib.adj_register_equation(self.adjointer, equation.equation, cs)
    assert cs.value == 0

    if hasattr(equation, 'rhs_fn'):
      self.functions_registered.append(equation.rhs_fn)

  def record_variable(self, var, storage):
    '''record_variable(self, var, storage)

    This method records the provided variable according to the settings in storage.'''

    clib.adj_record_variable(self.adjointer, var.var, storage.storage_data)

    python_utils.incref(storage.vec)

  def to_html(self, filename, viztype):
    try:
      typecode = {"forward": 1, "adjoint": 2, "tlm": 3}[viztype]
    except IndexError:
      raise exceptions.LibadjointErrorInvalidInput, "Argument viztype has to be one of the following: 'forward', 'adjoint', 'tlm'"

    if viztype == 'tlm':
      raise exceptions.LibadjointErrorNotImplemented, "HTML output for TLM is not implemented"

    clib.adj_adjointer_to_html(self.adjointer, filename, typecode)

  def get_forward_equation(self, equation):
    lhs = clib.adj_matrix()
    rhs = clib.adj_vector()
    fwd_var = clib.adj_variable()
    clib.adj_get_forward_equation(self.adjointer, equation, lhs, rhs, fwd_var)
    lhs_py = python_utils.c_deref(lhs.ptr)
    #python_utils.decref(lhs_py)
    rhs_py = python_utils.c_deref(rhs.ptr)
    #python_utils.decref(rhs_py)

    return (lhs_py, rhs_py)

  def get_adjoint_equation(self, equation, functional):
    lhs = clib.adj_matrix()
    rhs = clib.adj_vector()
    adj_var = clib.adj_variable()
    clib.adj_get_adjoint_equation(self.adjointer, equation, functional, lhs, rhs, adj_var)

  def __register_data_callbacks__(self):
    self.__register_data_callback__('ADJ_VEC_DUPLICATE_CB', self.__vec_duplicate_callback__)
    self.__register_data_callback__('ADJ_VEC_DESTROY_CB', self.__vec_destroy_callback__)
    self.__register_data_callback__('ADJ_VEC_AXPY_CB', self.__vec_destroy_callback__)
    self.__register_data_callback__('ADJ_MAT_DUPLICATE_CB', self.__mat_duplicate_callback__)
    self.__register_data_callback__('ADJ_MAT_DESTROY_CB', self.__mat_destroy_callback__)
    self.__register_data_callback__('ADJ_MAT_AXPY_CB', self.__mat_destroy_callback__)

  def __register_data_callback__(self, type_name, func):
    type_id = int(constants.adj_constants[type_name])

    try:
      cfunc = ctypes.CFUNCTYPE(None)
      clib.adj_register_data_callback(self.adjointer, ctypes.c_int(type_id), cfunc(func))
    except ctypes.ArgumentError:
      raise exceptions.LibadjointErrorInvalidInputs, 'Wrong function interface in register_data_callback for "%s".' % type_name 

  def __cfunc_from_block_assembly__(self, bassembly_cb):
    '''Given a block assembly function defined using the Pythonic interface, we want to translate that into a function
    that can be called from C. This routine does exactly that.'''

    def cfunc(ndepends_c, variables_c, dependencies_c, hermitian_c, coefficient_c, context_c, output_c, rhs_c):
      print "We've gotten into cfunc in block assembly"
      # build the Python objects from the C objects
      variables = [Variable(var=variables_c[i]) for i in range(ndepends_c)]
      dependencies = [vector(dependencies_c[i]) for i in range(ndepends_c)]
      hermitian = (hermitian_c == 1)
      coefficient = coefficient_c
      context = context_c

      # Now call the callback we've been given
      (matrix, rhs) = bassembly_cb(variables, dependencies, hermitian, coefficient, context)

      # Now cast the outputs back to C
      assert matrix is not None
      output_c[0].ptr = python_utils.c_ptr(matrix)
      python_utils.incref(matrix)

      assert rhs is not None
      rhs_c[0].ptr = python_utils.c_ptr(rhs)
      python_utils.incref(rhs)

    return self.block_assembly_type(cfunc)

  def register_block_assembly_callback(self, name, bassembly_cb):
    clib.adj_register_operator_callback.argtypes = [ctypes.POINTER(clib.adj_adjointer), ctypes.c_int, ctypes.c_char_p, self.block_assembly_type]
    fn = self.__cfunc_from_block_assembly__(bassembly_cb)
    self.functions_registered.append(fn)
    clib.adj_register_operator_callback(self.adjointer, int(constants.adj_constants['ADJ_BLOCK_ASSEMBLY_CB']), name, fn)
    clib.adj_register_operator_callback.argtypes = [ctypes.POINTER(clib.adj_adjointer), ctypes.c_int, ctypes.c_char_p, ctypes.CFUNCTYPE(None)]

  @staticmethod
  def __vec_duplicate_callback__(adj_vec, adj_vec_ptr):
    vec = vector(adj_vec)
    new_vec = vec.duplicate()

    # Increase the reference counter of the new object to protect it from deallocation at the end of the callback
    python_utils.incref(new_vec)
    adj_vec_ptr.ptr = python_utils.c_ptr(new_vec)

  @staticmethod
  def __vec_destroy_callback__(adj_vec_ptr):
    vec = vector(adj_vec_ptr.ptr)
    vec.destroy()

    # And do the corresponding decref of the object, so that the Python GC can pick it up
    python_utils.decref(vec)

  @staticmethod
  def __vec_axpy_callback__(adj_vec_ptr, alpha, adj_vec):
    y = vector(adj_vec_ptr.ptr)
    x = vector(adj_vec)
    y.axpy(alpha, x)

  @staticmethod
  def __mat_duplicate_callback__(adj_mat, adj_mat_ptr):
    mat = matrix(adj_mat)
    new_mat = mat.duplicate()

    # Increase the reference counter of the new object to protect it from deallocation at the end of the callback
    python_utils.incref(new_mat)
    adj_mat_ptr.ptr = python_utils.c_ptr(new_mat)

  @staticmethod
  def __mat_destroy_callback__(adj_mat_ptr):
    mat = matrix(adj_mat_ptr.ptr)
    mat.destroy()

    # And do the corresponding decref of the object, so that the Python GC can pick it up
    python_utils.decref(mat)

  @staticmethod
  def __mat_axpy_callback__(adj_mat_ptr, alpha, adj_mat):
    y = matrix(adj_mat_ptr.ptr)
    x = matrix(adj_mat)
    y.axpy(alpha, x)

class LinAlg(object):
  '''Base class for adjoint vector or matrix objects. In libadjoint,
  the operations performed on these are quite similar, so the common ones
  are factored out here. User applications should subclass Vector and Matrix
  instead of this directly.'''

  def __init__(self):
    raise LibadjointErrorNotImplemented("Shouldn't ever instantiate a LinAlg object directly")

  def destroy(self):
    '''destroy(self)

    Any special tasks to clean up the object. If this is a normal Python object, there is no need to override this method.'''
    pass

  def duplicate(self):
    '''duplicate(self)

    This method must return a newly allocated duplicate of its
    parent. The value of every entry of the duplicate must be zero.'''
    raise LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no copy() method')

  def axpy(self, alpha, x):
    '''axpy(self, alpha, x)

    This method must update the Vector with self=self+alpha*x where
    alpha is a scalar and x is a Vector'''

    raise LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no axpy(alpha,x) method')

class Vector(LinAlg):
  '''Base class for adjoint vector objects. User applications should
  subclass this and provide their own data and methods.'''

  def set_values(self, scalars):
    '''set_values(self, scalars)

    This method must set the value of Vector to that given by the array
    of scalars.'''

    raise LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no set_values(scalars) method')

  def size(self):
    '''size(self)

    This method must return the number of degrees of freedom in this Vector.'''

    raise LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no size() method')    

  def norm(self):
    '''norm(self)

    This method must return a norm for this vector. It does not matter which
    norm is chosen, so long as it satisfies the usual axioms for a norm.'''

    raise LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no norm() method')    

  def set_random(self):
    '''This method must set the entries of a given vector x to pseudo-random values.'''

    raise LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no set_random() method')        

  def dot_product(self, b):
    '''This method must return the result of dot(self, b).'''

    raise LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no dot_product() method')        


  def as_adj_vector(self):
    '''as_adj_vector(self)

    Returns an adj_vector with this Vector as its data payload.'''

    adj_vec = clib.adj_vector(ptr=python_utils.c_ptr(self))
    return adj_vec

class Matrix(LinAlg):
  '''Base class for adjoint matrix objects.'''

  def as_adj_matrix(self):
    '''as_adj_matrix(self)

    Returns an adj_matrix with this Matrix as its data payload.'''

    adj_mat = clib.adj_matrix(ptr=python_utils.c_ptr(self))
    return adj_mat

def vector(adj_vector):
  '''vector(adj_vector)

  Return the Python Vector object contained in a C adj_vector'''

  return python_utils.c_deref(adj_vector.ptr)

def matrix(adj_matrix):
  '''matrix(adj_matrix)

  Return the Python Matrix object contained in a C adj_matrix'''

  return python_utils.c_deref(adj_matrix.ptr)