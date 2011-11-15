import clibadjoint as clib
import clibadjoint_constants
import ctypes
import libadjoint_exceptions

def handle_error(ierr):
  if ierr != 0:
    exception = libadjoint_exceptions.get_exception(ierr)
    errstr  = clib.adj_get_error_string(ierr)
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


class Adjointer(object):
  def __init__(self):
    self.adjointer = clib.adj_adjointer()
    clib.adj_create_adjointer(self.adjointer)
    self.c_object = self.adjointer

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

  def to_html(self, filename, viztype):
    try:
      typecode = {"forward": 1, "adjoint": 2, "tlm": 3}[viztype]
    except IndexError:
      raise exceptions.LibadjointErrorInvalidInput, "Argument viztype has to be one of the following: 'forward', 'adjoint', 'tlm'"

    if viztype == 'tlm':
      raise exceptions.LibadjointErrorNotImplemented, "HTML output for TLM is not implemented"

    clib.adj_adjointer_to_html(self.adjointer, filename, typecode)

  def register_data_callback(self, type_name, func):
    try:
      index = zip(*clib.adj_data_callbacks._fields_)[0].index(type_name)
    except ValueError:
      raise libadjoint_exceptions.LibadjointErrorInvalidInputs, 'Wrong data callback type name in register_data_callback. Valid names are: "%s".' % '", "'.join(str(i) for i in zip(*clib.adj_data_callbacks._fields_)[0])
      return

    cfunctiontype = zip(*clib.adj_data_callbacks._fields_)[1][index]
    type_id = int(clibadjoint_constants.adj_constants['ADJ_'+type_name.upper()+'_CB'])

    try:
      cfunc = cfunctiontype(func)
      clib.adj_register_data_callback(self.adjointer, c_int(type_id), cfunc)
    except ctypes.ArgumentError:
      raise libadjoint_exceptions.LibadjointErrorInvalidInputs, 'Wrong function interface in register_data_callback for "%s".' % type_name 
      return
