import libadjoint.clibadjoint_constants as constants
import ctypes
import libadjoint.exceptions as exceptions

def handle_error(ierr):
  if ierr != 0:
    exception = libadjoint_exceptions.get_exception(ierr)
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
  def __init__(self, name, timestep, iteration=0, auxiliary=False):
    self.var = clib.adj_variable()
    self.name = name
    clib.adj_create_variable(name, timestep, iteration, auxiliary, self.var)
    self.c_object = self.var

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
      c_context = byref(context)
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
      c_context = byref(context)

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
  def __init__(self, var, blocks, targets, rhs_deps=None, rhs_context=None):
    self.equation = clib.adj_equation()
    self.c_object = self.equation

    assert len(blocks) == len(targets)

    clib.adj_create_equation(var.var, len(blocks), list_to_carray(blocks, clib.adj_block), list_to_carray(targets, clib.adj_variable), self.equation)

    if rhs_deps is not None and len(rhs_deps) > 0:
      clib.adj_equation_set_rhs_dependencies(self.equation, len(rhs_deps), list_to_carray(rhs_deps, clib.adj_variable), rhs_context)

  def __del__(self):
    clib.adj_destroy_equation(self.equation)

class Storage(object):
  '''Wrapper class for Vectors that contains additional information about how and where libadjoint 
  stores the vector values. This should never be used directly, instead use MemoryStorage or DiskStorage.'''

  storage_data = clib.adj_storage_data

  def set_compare(self, tol):
    if tol>0.0:
      adj_storage_set_compare(c_ptr(self.storage_data), c_int(1), c_double(tol))
    else:
      adj_storage_set_compare(c_ptr(self.storage_data), c_int(0), c_double(tol))

  def set_overwrite(self, flag=True):
    if flag:
      adj_storage_set_overwrite(c_ptr(self.storage_data), c_int(1))
    else:
      adj_storage_set_overwrite(c_ptr(self.storage_data), c_int(0))

  def set_checkpoint(self, flag=True):
    if flag:
      adj_storage_set_checkpoint(c_ptr(self.storage_data), c_int(1))
    else:
      adj_storage_set_checkpoint(c_ptr(self.storage_data), c_int(0))

class MemoryStorage(Storage):
  '''Wrapper class for Vectors that contains additional information for storing the vector values in memory.'''
  def __init__(self, vec):
    self.adj_storage_data.value = c_ptr(vec)
    adj_storage_memory_incref(vec, c_ptr(self.storage_data))
    incref(vec)

  def __del__(self):
    vec = deref(self.adj_storage_data.value)
    decref(vec)

class DiskStorage(Storage):
  '''Wrapper class for Vectors that contains additional information for storing the vector values in memory.'''
  def __init__(self, vec):
    self.adj_storage_data.value = c_ptr(vec)
    adj_storage_disk_incref(vec, c_ptr(self.storage_data))
    incref(vec)

  def __del__(self):
    vec = deref(self.adj_storage_data.value)
    decref(vec)

class Adjointer(object):
  def __init__(self):
    self.adjointer = clib.adj_adjointer()
    clib.adj_create_adjointer(self.adjointer)
    self.c_object = self.adjointer
    self.__register_data_callbacks__()

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

  def record_variable(self, var, storage):
    '''record_variable(self, var, storage)

    This method records the provided variable according to the settings in storage.'''

    adj_var = c_ptr(var.ptr)
    adj_record_variable(self.adjointer, adj_var, storage.storage_data)

  def to_html(self, filename, viztype):
    try:
      typecode = {"forward": 1, "adjoint": 2, "tlm": 3}[viztype]
    except IndexError:
      raise exceptions.LibadjointErrorInvalidInput, "Argument viztype has to be one of the following: 'forward', 'adjoint', 'tlm'"

    if viztype == 'tlm':
      raise exceptions.LibadjointErrorNotImplemented, "HTML output for TLM is not implemented"

    clib.adj_adjointer_to_html(self.adjointer, filename, typecode)

  def __register_data_callbacks__(self):
    self.__register_data_callback__('ADJ_VEC_DUPLICATE_CB', self.__vec_duplicate_callback__)

  def __register_data_callback__(self, type_name, func):
    type_id = int(constants.adj_constants[type_name])

    try:
      cfunc = ctypes.CFUNCTYPE(None)
      clib.adj_register_data_callback(self.adjointer, ctypes.c_int(type_id), cfunc(func))
    except ctypes.ArgumentError:
      raise libadjoint_exceptions.LibadjointErrorInvalidInputs, 'Wrong function interface in register_data_callback for "%s".' % type_name 

  def __vec_duplicate_callback__(self, adj_vec, adj_vec_ptr):
    vec = c_derefr(adj_vec.ptr)
    new_vec = vec.duplicate()

    # Increase the reference counter of the new object to protect it from deallocation at the end of the callback
    incref(new_vec)
    adj_vec_ptr.ptr = c_ptr(new_vec)


class Vector(object):
  '''Base class for adjoint vector objects. User applications should
  subclass this and provide their own data'''
  def __init__(self):
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
  
  def destroy(self):
    '''destroy(self)

    This method must cause the Vector to be deallocated.'''

    raise LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no destroy() method')
  
  def set_values(self, scalars):
    '''set_values(self, scalars)

    This method must set the value of Vector to that given by the array
    of scalars.'''

    raise LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no set_values(scalars) method')

#  def get_size(
