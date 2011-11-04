import clibadjoint as clib
import ctypes
import libadjoint_exceptions

def handle_error(ierr):
  if ierr != 0:
    exception = libadjoint_exceptions.get_exception(ierr)
    errstr  = clib.adj_get_error_string(ierr)
    raise exception, errstr

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

    deplisttype = clib.adj_variable * len(dependencies)
    deplist = deplisttype()
    for i in range(len(dependencies)):
      deplist[i] = dependencies[i].var

    clib.adj_create_nonlinear_block(name, len(dependencies), deplist, c_context, self.nblock)

    if coefficient is not None:
      self.set_coefficient(coefficient)

  def __del__(self):
    clib.adj_destroy_nonlinear_block(self.nblock)

  def set_coefficient(self, c):
    clib.adj_nonlinear_block_set_coefficient(self.nblock, c)

class Block(object):
  def __init__(self, name, nblock=None, context=None, coefficient=None, hermitian=False):
    self.block = clib.adj_block()
    c_context = None
    if context is not None:
      c_context = byref(context)

    if nblock is None:
      c_nblock = None
    else:
      c_nblock = nblock.nblock

    clib.adj_create_block(name, c_nblock, c_context, self.block)

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
  def __init__(self, var, blocks, targets):
    self.equation = clib.adj_equation()

    assert len(blocks) == len(targets)

    blocklisttype = clib.adj_block * len(blocks)
    blocklist = blocklisttype()
    for i in range(len(blocks)):
      blocklist[i] = blocks[i].block

    targetlisttype = clib.adj_variable * len(targets)
    targetlist = targetlisttype()
    for i in range(len(targets)):
      targetlist[i] = targets[i].var

    clib.adj_create_equation(var.var, len(blocks), blocklist, targetlist, self.equation)

  def __del__(self):
    clib.adj_destroy_equation(self.equation)


class Adjointer(object):
  def __init__(self):
    self.adjointer = clib.adj_adjointer()
    clib.adj_create_adjointer(self.adjointer)

  def __del__(self):
    clib.adj_destroy_adjointer(self.adjointer)

  def register_equation(self, equation):
    clib.adj_register_equation(self.adjointer, equation.equation)

  def to_html(self, filename, viztype):
    try:
      typecode = {"forward": 1, "adjoint": 2, "tlm": 3}[viztype]
    except IndexError:
      raise exceptions.LibadjointErrorInvalidInput, "Argument viztype has to be one of the following: 'forward', 'adjoint', 'tlm'"

    if viztype == 'tlm':
      raise exceptions.LibadjointErrorNotImplemented, "HTML output for TLM is not implemented"

    clib.adj_adjointer_to_html(self.adjointer, filename, typecode)

