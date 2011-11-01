import clibadjoint as clib
import ctypes
import exceptions

def handle_error(ierr):
  if ierr != 0:
    exception = exceptions.get_exception(ierr)
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
    if name == "name":
      return self.name
    elif name == "timestep":
      timestep = ctypes.c_int()
      clib.adj_variable_get_timestep(self.var, timestep)
      return timestep.value
    elif name == "iteration":
      iteration = ctypes.c_int()
      clib.adj_variable_get_iteration(self.var, iteration)
      return iteration.value
    else:
      raise AttributeError

class Adjointer(object):
  def __init__(self):
    self.adjointer = clib.adj_adjointer()
    clib.adj_create_adjointer(self.adjointer)

  def __del__(self):
    clib.adj_destroy_adjointer(self.adjointer)
