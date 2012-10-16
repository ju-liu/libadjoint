import clibadjoint_constants as constants
import ctypes
import exceptions
import clibadjoint as clib
import python_utils

adj_scalar = ctypes.c_double
references_taken = []

def adj_test_assert(test_pass, msg=None):
  assert isinstance(test_pass, bool)
  clib.adj_test_assert(ctypes.c_int(test_pass), msg)

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
handler_exceptions = ["adj_variable_equal"]
for member in dir(clib):
  # Looping over all of the objects this module offers us ...
  if member.startswith("adj_") and member not in handler_exceptions:
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

  def copy(self):
     ''' Creates a deep copy of the variable. '''
     return Variable(self.var.name, self.var.timestep, self.var.iteration, self.var.auxiliary)

  def iteration_count(self, adjointer):
    '''Returns the number of iterations at the variables timestep'''
    iteration_count = ctypes.c_int()
    clib.adj_iteration_count(adjointer.adjointer, self.var, iteration_count)
    return iteration_count.value

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
    elif name == 'type':
      int_type_map = {1: 'ADJ_FORWARD', 2: 'ADJ_ADJOINT', 3: 'ADJ_TLM'}
      return int_type_map[self.c_object.type]
    else:
      raise AttributeError

  def __setattr__(self, key, val):
    if key == "timestep":
      self.var.timestep = val
    elif key == "iteration":
      self.var.iteration = val
    elif key == "type":
      type_map = {'ADJ_FORWARD': 1, 'ADJ_ADJOINT': 2, 'ADJ_TLM': 3}
      self.var.type = type_map[val]
    else:
      object.__setattr__(self, key, val)

  def __eq__(self, other):

    if not isinstance(other, Variable): 
      return False
    else:
      return (clib.adj_variable_equal(self.var, other.var, 1)==1)

  def to_adjoint(self, functional):
    adj_var = Variable(self.name, self.timestep, self.iteration)
    adj_var.c_object.type = int(constants.adj_constants['ADJ_ADJOINT'])
    adj_var.c_object.functional = str(functional)
    return adj_var

  def to_tlm(self, parameter):
    tlm_var = Variable(self.name, self.timestep, self.iteration)
    tlm_var.c_object.type = int(constants.adj_constants['ADJ_TLM'])
    tlm_var.c_object.functional = str(parameter)
    return tlm_var

  def to_forward(self):
    fwd_var = Variable(self.name, self.timestep, self.iteration)
    return fwd_var

  def equation_nb(self, adjointer):
    ''' Returns the number of the foward equation that solves for this variable '''
    equation_nb = ctypes.c_int()
    clib.adj_find_variable_equation_nb(adjointer.adjointer, self.var, equation_nb)
    return equation_nb.value

class NonlinearBlock(object):
  def __init__(self, name, dependencies, context=None, coefficient=None, test_hermitian=None, test_derivative=None):
    self.name = name
    self.nblock = clib.adj_nonlinear_block()
    c_context = None
    if context is not None:
      c_context = ctypes.byref(context)
    self.c_object = self.nblock

    clib.adj_create_nonlinear_block(name, len(dependencies), list_to_carray(dependencies, clib.adj_variable), c_context, 1.0, self.nblock)

    if coefficient is not None:
      self.set_coefficient(coefficient)

    if test_hermitian is not None and test_hermitian is not False:
      if test_hermitian is True:
        number_of_tests = 10
        tolerance = 1.0e-14
      else:
        number_of_tests = test_hermitian[0]
        tolerance = test_hermitian[1]
      self.set_test_hermitian(number_of_tests, tolerance)

    if test_derivative is not None and test_derivative is not False:
      if test_derivative is True:
        number_of_rounds = 6
      else:
        number_of_rounds = test_derivative
      self.set_test_derivative(number_of_rounds)

  def __del__(self):
    clib.adj_destroy_nonlinear_block(self.nblock)

  def set_coefficient(self, c):
    clib.adj_nonlinear_block_set_coefficient(self.nblock, c)

  def set_test_hermitian(self, number_of_tests, tolerance):
    clib.adj_nonlinear_block_set_test_hermitian(self.nblock, 1, number_of_tests, tolerance)

  def set_test_derivative(self, number_of_rounds):
    if number_of_rounds is True:
      number_of_rounds = 5 # default
    clib.adj_nonlinear_block_set_test_derivative(self.nblock, 1, number_of_rounds)

class Block(object):
  def __init__(self, name, nblock=None, context=None, coefficient=None, hermitian=False, dependencies=None, test_hermitian=None, test_derivative=None):
    self.block = clib.adj_block()
    self.c_object = self.block
    self.name = name

    c_context = None
    if context is not None:
      c_context = ctypes.byref(context)

    if nblock is not None and dependencies is not None:
      raise exceptions.LibadjointErrorInvalidInputs, "Cannot have both nblock and dependencies"

    if nblock is None:
      if dependencies is not None and len(dependencies) > 0:
        nblock = NonlinearBlock(name, dependencies, context=context, test_hermitian=test_hermitian, test_derivative=test_derivative)

    if nblock is not None:
      self.nblock = nblock
      c_nblock = nblock.nblock
    else:
      c_nblock = None

    clib.adj_create_block(name, c_nblock, c_context, 1.0, self.block)

    if coefficient is not None:
      self.set_coefficient(coefficient)

    if hermitian:
      self.set_hermitian(hermitian)

    if test_hermitian is not None and test_hermitian is not False:
      if test_hermitian is True:
        number_of_tests = 10
        tolerance = 1.0e-14
      else:
        number_of_tests = test_hermitian[0]
        tolerance = test_hermitian[1]
      self.set_test_hermitian(number_of_tests, tolerance)

  def __del__(self):
    clib.adj_destroy_block(self.block)

  def set_coefficient(self, c):
    clib.adj_block_set_coefficient(self.block, c)

  def set_hermitian(self, hermitian):
    clib.adj_block_set_hermitian(self.block, hermitian)

  def set_test_hermitian(self, number_of_tests, tolerance):
    clib.adj_block_set_test_hermitian(self.block, 1, number_of_tests, tolerance)

  @staticmethod
  def assemble(dependencies, values, hermitian, coefficient, context):
    '''def assemble(dependencies, values, hermitian, coefficient, context)

    This method must assemble the block and return a tuple (matrix, rhs) where matrix.

      dependencies is a list of Variable objects giving the variables on which the block depends.

      values is a list of Vector objects giving the values of those values.

      hermitian is a boolean value indicating whether the hermitian of the operator is to be constructed.

      coefficient is a coefficient by which the routine must scale the output.

      context is the python object passed to the Block on construction.
    '''

    # The registration code will notice unimplemented methods and fail to register them.
    pass


  @staticmethod
  def action(dependencies, values, hermitian, coefficient, input, context):
    '''def action(dependencies, values, hermitian, coefficient, input, context)

    If hermitian is False, this method must return:
                    coefficient * dot(block, input)
    If hermitian is True, this method must return:
                    coefficient * dot(block*, input)

      dependencies is a list of Variable objects giving the variables on which the block depends.

      values is a list of Vector objects giving the values of those values.

      hermitian is a boolean value indicating whether the hermitian of the operator is to be constructed.

      coefficient is a coefficient by which the routine must scale the output.

      input is a Vector which is is to be the subject of the action.

      context is the python object passed to the Block on construction.
    '''

    # The registration code will notice unimplemented methods and fail to register them.
    pass

  @staticmethod
  def derivative_action(dependencies, values, variable, contraction_vector, hermitian, input, coefficient, context):
    '''def derivative_action(dependencies, values, variable, contraction_vector, hermitian, input, coefficient, context)

    If hermitian is False, this method must return:
                    coefficient * (diff(block, variable).contraction) . input
    If hermitian is True, this method must return:
                    coefficient * (diff(block, variable).contraction)^* . input
    See the libadjoint manual for more details (ADJ_NBLOCK_DERIVATIVE_ACTION_CB).

      dependencies is a list of Variable objects giving the variables on which the block depends.

      values is a list of Vector objects giving the values of those values.

      variable is the Variable with respect to which the block must be differentiated.

      contraction_vector is the Vector with which the derivative operation is to be contracted.

      hermitian is a boolean value indicating whether the hermitian of the operator is to be constructed.

      input is a Vector which is is to be the subject of the action.

      coefficient is a coefficient by which the routine must scale the output.

      context is the python object passed to the Block on construction.
    '''

    # The registration code will notice unimplemented methods and fail to register them.
    pass


class Equation(object):
  def __init__(self, var, blocks, targets, rhs=None):
    self.equation = clib.adj_equation()
    self.c_object = self.equation
    self.blocks=blocks
    self.var=var

    assert len(blocks) == len(targets)

    clib.adj_create_equation(var.var, len(blocks), list_to_carray(blocks, clib.adj_block), list_to_carray(targets, clib.adj_variable), self.equation)

    if (rhs is not None):
      self.rhs=rhs
      rhs.register(self)


  def __del__(self):
    clib.adj_destroy_equation(self.equation)


class Storage(object):
  '''Wrapper class for Vectors that contains additional information about how and where libadjoint 
  stores the vector values. This should never be used directly, instead use MemoryStorage or DiskStorage.'''

  def set_compare(self, tol=None):
    if tol is not None:
      clib.adj_storage_set_compare(self.storage_data, 1, tol)
    else:
      clib.adj_storage_set_compare(self.storage_data, 0, 0.0)

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
  def __init__(self, vec, copy=True, cs=False):
    self.storage_data = clib.adj_storage_data()
    self.c_object = self.storage_data
    if copy:
      clib.adj_storage_memory_copy(vec.as_adj_vector(), self.storage_data)
    else:
      adjvec = vec.as_adj_vector()
      clib.adj_storage_memory_incref(adjvec, self.storage_data)
      references_taken.append(adjvec)

    self.set_checkpoint(cs)

    # Ensure that the storage object always holds a reference to the vec
    self.vec = vec

class DiskStorage(Storage):
  '''Wrapper class for Vectors that contains additional information for storing the vector values in memory.'''
  def __init__(self, vec, cs=False):
    self.storage_data = clib.adj_storage_data()
    self.c_object = self.storage_data
    clib.adj_storage_disk(vec.as_adj_vector(), self.storage_data)

    self.set_checkpoint(cs)

    # Ensure that the storage object always holds a reference to the vec
    self.vec = vec

class Functional(object):
  '''Base class for functionals and their derivatives.'''
  def __init__(self):
    pass

  def __call__(self, adjointer, timestep, dependencies, values):
    '''__call__(self, dependencies, values)

    Evaluate functional given dependencies with values. The result must be a scalar.
    '''

    raise exceptions.LibadjointErrorNotImplemented("No __call__ method provided for functional.")

  def __str__(self):

    return hex(id(self))

  def derivative(self, adjointer, variable, dependencies, values):
    '''derivative(self, variable, dependencies, values)

    Evaluate the derivative of the functional with respect to variable given dependencies with values. The result will be a Vector.
    '''

    raise exceptions.LibadjointErrorNotImplemented("No derivative method provided for functional.")

  def dependencies(self, adjointer, timestep):
    '''dependencies(self, adjointer, timestep)

    Return the list of Variables on which this functional depends at timestep. The adjointer may be queried to find the actual value of time if required.'''

    raise exceptions.LibadjointErrorNotImplemented("No dependencies method provided for functional.")

  def __cfunc_from_derivative__(self):
    '''Return a c-callable function wrapping the derivative method.'''

    def cfunc(adjointer_c, variable_c, ndepends_c, dependencies_c, values_c, name_c, output_c):
      # build the Python objects from the C objects
      adjointer = Adjointer(adjointer_c)
      variable  = Variable(var=variable_c)
      dependencies = [Variable(var=dependencies_c[i]) for i in range(ndepends_c)]
      values = [vector(values_c[i]) for i in range(ndepends_c)]

      # Now call the callback we've been given
      output = self.derivative(adjointer, variable, dependencies, values)

      # Now cast the outputs back to C
      output_c[0].klass = 0
      output_c[0].flags = 0
      output_c[0].ptr = 0
      if not isinstance(output, Vector):
        raise exceptions.LibadjointErrorInvalidInputs("Output from functional derivative must be a Vector.")
      output_c[0].ptr = python_utils.c_ptr(output)
      references_taken.append(output)

    functional_derivative_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_adjointer), clib.adj_variable, ctypes.c_int, ctypes.POINTER(clib.adj_variable), ctypes.POINTER(clib.adj_vector), ctypes.c_char_p, ctypes.POINTER(clib.adj_vector))
    return functional_derivative_type(cfunc)

  def __cfunc_from_functional__(self):
    '''Return a c-callable function wrapping the derivative method.'''

    def cfunc(adjointer_c, timestep, ndepends_c, dependencies_c, values_c, name_c, output_c):
      # build the Python objects from the C objects
      adjointer = Adjointer(adjointer_c)
      dependencies = [Variable(var=dependencies_c[i]) for i in range(ndepends_c)]
      values = [vector(values_c[i]) for i in range(ndepends_c)]

      # Now call the callback we've been given
      output_c[0] = self(adjointer, timestep, dependencies, values)

    functional_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_adjointer), ctypes.c_int, ctypes.c_int, ctypes.POINTER(clib.adj_variable), ctypes.POINTER(clib.adj_vector), ctypes.c_char_p, ctypes.POINTER(ctypes.c_double))
    return functional_type(cfunc)

class Parameter(object):
  '''Base class for parameters and their source terms.'''
  def __init__(self):
    pass

  def __call__(self, adjointer, equation, dependencies, values, variable):
    '''__call__(self, dependencies, values)

    Evaluate dF/dm associated with the equation for variable, given dependencies with values. The result must be a Vector.
    '''

    raise exceptions.LibadjointErrorNotImplemented("No __call__ method provided for parameter.")

  def __str__(self):

    return hex(id(self))

  def dependencies(self, adjointer, variable):
    '''dependencies(self, adjointer, variable)

    Return the list of Variables on which this parameter's source term depends on for the equation associated with variable.'''

    raise exceptions.LibadjointErrorNotImplemented("No dependencies method provided for variable.")

  def __cfunc_from_parameter_source__(self):
    '''Return a c-callable function wrapping the parameter source method.'''

    def cfunc(adjointer_c, equation_c, variable_c, ndepends_c, dependencies_c, values_c, name_c, output_c, has_output_c):
      # build the Python objects from the C objects
      adjointer = Adjointer(adjointer_c)
      equation = equation_c
      variable  = Variable(var=variable_c)
      dependencies = [Variable(var=dependencies_c[i]) for i in range(ndepends_c)]
      values = [vector(values_c[i]) for i in range(ndepends_c)]

      # Now call the callback we've been given
      output = self(adjointer, equation, dependencies, values, variable)

      has_output_c[0] = (output is not None)
      output_c[0].klass = 0
      output_c[0].flags = 0
      output_c[0].ptr = 0
      if output:
        if not isinstance(output, Vector):
          raise exceptions.LibadjointErrorInvalidInputs("Output from parameter source term must be a Vector.")
        output_c[0].ptr = python_utils.c_ptr(output)
        references_taken.append(output)

    parameter_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_adjointer), ctypes.c_int, clib.adj_variable, ctypes.c_int, ctypes.POINTER(clib.adj_variable), ctypes.POINTER(clib.adj_vector), ctypes.c_char_p, ctypes.POINTER(clib.adj_vector), ctypes.POINTER(ctypes.c_int))
    return parameter_type(cfunc)

class RHS(object):
  '''Base class for equation Right Hand Sides and their derivatives.'''
  def __init__(self):
    pass

  def __call__(self, dependencies, values):
    '''__call__(self, dependencies, values)

    Evaluate RHS given the values corresponding to dependencies. The result must be a Vector.
    '''

    raise exceptions.LibadjointErrorNotImplemented("No __call__ method provided for RHS.")

  def __str__(self):

    return hex(id(self))

  def derivative_action(self, dependencies, values, variable, contraction_vector, hermitian):
    '''derivative_action(self, dependencies, values, variable, contraction_vector, hermitian):

    Evaluate the action of the derivative of the RHS with respect to variable given dependencies with values on contraction_vector. The result will be a Vector.
    '''

    raise exceptions.LibadjointErrorNotImplemented("No derivative action method provided for RHS.")

  def derivative_assembly(self, dependencies, values, variable, hermitian):
    '''derivative_assembly(self, dependencies, values, variable, hermitian):

    Evaluate the derivative of the RHS with respect to variable given dependencies with values on contraction_vector. The result will be a Matrix.
    '''

    raise exceptions.LibadjointErrorNotImplemented("No derivative assembly method provided for RHS.")

  def dependencies(self):
    '''dependencies(self)

    Return the list of Variables on which this functional depends.'''

    raise exceptions.LibadjointErrorNotImplemented("No dependencies method provided for RHS.")

  def register(self, equation):
    '''register(self, equation)

    Register this RHS as the RHS of equation.'''

    rhs_deps = self.dependencies()
    if rhs_deps is not None and len(rhs_deps) > 0:
      clib.adj_equation_set_rhs_dependencies(equation.equation, len(rhs_deps), list_to_carray(rhs_deps, clib.adj_variable), None)

    equation.rhs_fn = self.__cfunc_from_rhs__()
    clib.adj_equation_set_rhs_callback(equation.equation, equation.rhs_fn)

    equation.rhs_derivative_action_fn = self.__cfunc_from_derivative_action__()
    clib.adj_equation_set_rhs_derivative_action_callback(equation.equation, equation.rhs_derivative_action_fn)

    equation.rhs_derivative_assembly_fn = self.__cfunc_from_derivative_assembly__()
    clib.adj_equation_set_rhs_derivative_assembly_callback(equation.equation, equation.rhs_derivative_assembly_fn)

  def __cfunc_from_rhs__(self):
    '''Given a rhs function defined using the Pythonic interface, we want to translate that into a function that
    can be called from C. This routine does exactly that.'''

    def cfunc(adjointer_c, variable_c, ndepends_c, dependencies_c, values_c, context_c, output_c, has_output_c):
      # build the Python objects from the C objects
      variable  = Variable(var=variable_c)
      dependencies = [Variable(var=dependencies_c[i]) for i in range(ndepends_c)]
      values = [vector(values_c[i]) for i in range(ndepends_c)]

      # Now call the callback we've been given
      output = self(dependencies, values)

      # Now cast the outputs back to C
      has_output_c[0] = (output is not None)
      output_c[0].klass = 0
      output_c[0].flags = 0
      output_c[0].ptr = 0
      if output:
        if not isinstance(output, Vector):
          raise exceptions.LibadjointErrorInvalidInputs("Output from RHS callback must be None or a Vector.")
        output_c[0].ptr = python_utils.c_ptr(output)

      references_taken.append(output)

    rhs_func_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_adjointer), clib.adj_variable, ctypes.c_int, ctypes.POINTER(clib.adj_variable), ctypes.POINTER(clib.adj_vector), ctypes.POINTER(None),
                                     ctypes.POINTER(clib.adj_vector), ctypes.POINTER(ctypes.c_int))

    return rhs_func_type(cfunc)


  def __cfunc_from_derivative_action__(self):
    '''Return a c-callable function wrapping the derivative_action method.'''

    def cfunc(adjointer_c, source_variable_c, ndepends_c, dependencies_c, values_c, variable_c, contraction_c, hermitian_c, context_c, output_c, has_output_c):
      # build the Python objects from the C objects
      variable  = Variable(var=variable_c)
      dependencies = [Variable(var=dependencies_c[i]) for i in range(ndepends_c)]
      values = [vector(values_c[i]) for i in range(ndepends_c)]
      hermitian = (hermitian_c==1)
      contraction_vector = vector(contraction_c)

      # Now call the callback we've been given
      output = self.derivative_action(dependencies, values, variable, contraction_vector, hermitian)

      # Now cast the outputs back to C
      has_output_c[0] = (output is not None)
      output_c[0].klass = 0
      output_c[0].flags = 0
      output_c[0].ptr = 0
      if output:
        if not isinstance(output, Vector):
          raise exceptions.LibadjointErrorInvalidInputs("Output from RHS derivative_action callback must be None or a Vector.")
        output_c[0].ptr = python_utils.c_ptr(output)

      references_taken.append(output)

    rhs_deriv_action_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_adjointer), clib.adj_variable, ctypes.c_int, ctypes.POINTER(clib.adj_variable),
        ctypes.POINTER(clib.adj_vector), clib.adj_variable, clib.adj_vector, ctypes.c_int, ctypes.POINTER(None), ctypes.POINTER(clib.adj_vector), ctypes.POINTER(ctypes.c_int))
    return rhs_deriv_action_type(cfunc)

  def __cfunc_from_derivative_assembly__(self):
    '''Return a c-callable function wrapping the derivative_assembly method.'''

    def cfunc(adjointer_c, source_variable_c, ndepends_c, dependencies_c, values_c, hermitian_c, context_c, output_c):
      # build the Python objects from the C objects
      variable  = Variable(var=source_variable_c)
      dependencies = [Variable(var=dependencies_c[i]) for i in range(ndepends_c)]
      values = [vector(values_c[i]) for i in range(ndepends_c)]
      hermitian = (hermitian_c==1)

      # Now call the callback we've been given
      output = self.derivative_assembly(dependencies, values, variable, hermitian)

      # Now cast the outputs back to C
      output_c[0].klass = 0
      output_c[0].flags = 0
      output_c[0].ptr = 0
      if not isinstance(output, Matrix):
        raise exceptions.LibadjointErrorInvalidInputs("Output from RHS derivative_assembly callback must be a Matrix.")
      output_c[0].ptr = python_utils.c_ptr(output)

      references_taken.append(output)

    rhs_deriv_assembly_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_adjointer), clib.adj_variable, ctypes.c_int, ctypes.POINTER(clib.adj_variable),
        ctypes.POINTER(clib.adj_vector), ctypes.c_int, ctypes.POINTER(None), ctypes.POINTER(clib.adj_matrix))
    return rhs_deriv_assembly_type(cfunc)

class AdjointerTime(object):
  """class to facilitate recording the simulation time at which timesteps occur"""
  def __init__(self, adjointer):
    self.time_levels = []
    self.finished = False
    self.adjointer = adjointer

  def start(self, time):
    """start(self, time) 
    Set the start time of the simulation."""
    
    if len(self.time_levels)!=0:
      raise exceptions.LibadjointErrorInvalidInputs(
        "time.start() called after simulation started!")

    self.time_levels = [time]
    
  def finish(self):
    """finish()
    Record that the annotation has finished."""
    self.finished = True
    clib.adj_set_finished(self.adjointer.adjointer, 1)
    
  def next(self, time):
    """Increment the timestep counter, and mark the end of the timestep."""
    
    if self.finished:
      raise exceptions.LibadjointErrorInvalidInputs(
        "time.next() called after simulation finished!")
      
    if self.time_levels==[]:
      raise exceptions.LibadjointErrorInvalidInputs(
        "time.next() called before time started!")

    self.time_levels.append(time)
    timestep = len(self.time_levels) - 2
    self.adjointer.set_times(timestep, self.time_levels[-2], self.time_levels[-1])

  def reset(self):
    self.__init__(self.adjointer)
    

class Adjointer(object):
  def __init__(self, adjointer=None):
    self.functions_registered = []
    self.set_function_apis()

    self.equation_timestep = []

    # This gets clobbered during casting.
    self.time = AdjointerTime(self)
    
    if adjointer is None:
      self.adjointer = clib.adj_adjointer()
      clib.adj_create_adjointer(self.adjointer)
      self.adjointer_created = True
      self.c_object = self.adjointer
      self.__register_data_callbacks__()
    else:
      self.adjointer_created = False
      self.adjointer = adjointer
      self.c_object = self.adjointer

  def set_checkpoint_strategy(self, strategy):
    try:
      strategy_id = int(constants.adj_constants['ADJ_CHECKPOINT_REVOLVE_' + strategy.upper()])
    except KeyError:
      raise libadjoint.exceptions.LibadjointErrorInvalidInputs("Unknown checkpointing strategy " + strategy + ". Known strategies: ['offline', 'online', 'multistage'].")
    clib.adj_set_checkpoint_strategy(self.adjointer, strategy_id)

  def set_revolve_options(self, steps, snaps_on_disk, snaps_in_ram, verbose=False):
      clib.adj_set_revolve_options(self.adjointer, steps, snaps_on_disk, snaps_in_ram, verbose)

  def set_revolve_debug_options(self, overwrite, comparison_tolerance):
      clib.adj_set_revolve_debug_options(self.adjointer, overwrite, comparison_tolerance)

  def check_checkpoints(self):
      clib.adj_adjointer_check_checkpoints(self.adjointer)

  def set_function_apis(self):
    self.block_assembly_type = ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.POINTER(clib.adj_variable), ctypes.POINTER(clib.adj_vector), ctypes.c_int, adj_scalar, ctypes.POINTER(None),
                                                ctypes.POINTER(clib.adj_matrix), ctypes.POINTER(clib.adj_vector))
    self.block_action_type = ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.POINTER(clib.adj_variable), ctypes.POINTER(clib.adj_vector), ctypes.c_int, adj_scalar, clib.adj_vector, 
                                              ctypes.POINTER(None),ctypes.POINTER(clib.adj_vector))
    self.nblock_derivative_action_type = ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.POINTER(clib.adj_variable), ctypes.POINTER(clib.adj_vector), clib.adj_variable, clib.adj_vector, ctypes.c_int, clib.adj_vector,
                                                          adj_scalar, ctypes.POINTER(None), ctypes.POINTER(clib.adj_vector))
    self.nblock_action_type = ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.POINTER(clib.adj_variable), ctypes.POINTER(clib.adj_vector), clib.adj_vector, ctypes.POINTER(None), ctypes.POINTER(clib.adj_vector))
    self.vec_duplicate_type = ctypes.CFUNCTYPE(None, clib.adj_vector, ctypes.POINTER(clib.adj_vector))
    self.vec_destroy_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_vector))
    self.vec_axpy_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_vector), adj_scalar, clib.adj_vector)
    self.vec_get_norm_type = ctypes.CFUNCTYPE(None, clib.adj_vector, ctypes.POINTER(adj_scalar))
    self.vec_dot_product_type = ctypes.CFUNCTYPE(None, clib.adj_vector, clib.adj_vector, ctypes.POINTER(adj_scalar))
    self.vec_set_random_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_vector))
    self.vec_set_values_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_vector), ctypes.POINTER(adj_scalar))
    self.vec_get_values_type = ctypes.CFUNCTYPE(None, clib.adj_vector, ctypes.POINTER(ctypes.POINTER(adj_scalar)))
    self.vec_get_size_type = ctypes.CFUNCTYPE(None, clib.adj_vector, ctypes.POINTER(ctypes.c_int))
    self.vec_write_type = ctypes.CFUNCTYPE(None, clib.adj_variable, clib.adj_vector)
    self.vec_read_type = ctypes.CFUNCTYPE(None, clib.adj_variable, ctypes.POINTER(clib.adj_vector))
    self.vec_delete_type = ctypes.CFUNCTYPE(None, clib.adj_variable)
    self.mat_duplicate_type = ctypes.CFUNCTYPE(None, clib.adj_matrix, ctypes.POINTER(clib.adj_matrix))
    self.mat_destroy_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_matrix))
    self.mat_action_type = ctypes.CFUNCTYPE(None, clib.adj_matrix, clib.adj_vector, ctypes.POINTER(clib.adj_vector))
    self.mat_axpy_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_matrix), adj_scalar, clib.adj_matrix)
    self.solve_type = ctypes.CFUNCTYPE(None, clib.adj_variable, clib.adj_matrix, clib.adj_vector, ctypes.POINTER(clib.adj_vector))
    self.functional_type = ctypes.CFUNCTYPE(None, ctypes.POINTER(clib.adj_adjointer), ctypes.c_int, ctypes.c_int, ctypes.POINTER(clib.adj_variable), ctypes.POINTER(clib.adj_vector), 
                                                  ctypes.c_char_p, ctypes.POINTER(clib.adj_vector))

  def reset(self):
    ''' Resets all time information and forgets the annotation. '''
    if self.adjointer_created:
      clib.adj_destroy_adjointer(self.adjointer)
      assert len(references_taken) == 0

    self.functions_registered = []
    self.equation_timestep=[]
    self.adjointer = clib.adj_adjointer()
    clib.adj_create_adjointer(self.adjointer)
    self.adjointer_created = True
    self.c_object = self.adjointer
    self.__register_data_callbacks__()
    self.time.reset()

  def __del__(self):
    if self.adjointer_created:
      clib.adj_destroy_adjointer(self.adjointer)
      if len(references_taken) != 0:
        print "Warning: references still exist!"
        print "References: ", references_taken
      assert len(references_taken) == 0

  def __getattr__(self, name):
    if name == "equation_count":
      equation_count = ctypes.c_int()
      clib.adj_equation_count(self.adjointer, equation_count)
      return equation_count.value
    elif name == "timestep_count":
      timestep_count = ctypes.c_int()
      clib.adj_timestep_count(self.adjointer, timestep_count)
      return timestep_count.value
    if name == "finished":
      finished = ctypes.c_int()
      clib.adj_get_finished(self.adjointer, finished)
      return finished.value
    else:
      raise AttributeError(name)

  def register_equation(self, equation):

    self.equation_timestep.append(equation.var.timestep)

    cs = ctypes.c_int()
    clib.adj_register_equation(self.adjointer, equation.equation, cs)

    if hasattr(equation, 'rhs_fn'):
      self.functions_registered.append(equation.rhs_fn)

    if hasattr(equation, 'rhs_derivative_action_fn'):
      self.functions_registered.append(equation.rhs_derivative_action_fn)

    if hasattr(equation, 'rhs_derivative_assembly_fn'):
      self.functions_registered.append(equation.rhs_derivative_assembly_fn)

    for block in equation.blocks:
      if not (block.assemble is Block.assemble):
        # There is an assemble method to register.
        self.__register_operator_callback__(block.name, "ADJ_BLOCK_ASSEMBLY_CB", block.assemble)

      if not (block.action is Block.action):
        # There is an action method to register.
        self.__register_operator_callback__(block.name, "ADJ_BLOCK_ACTION_CB", block.action)

      if block.block.has_nonlinear_block:
        if not (block.derivative_action is Block.derivative_action):
          self.__register_operator_callback__(block.nblock.name, "ADJ_NBLOCK_DERIVATIVE_ACTION_CB", block.derivative_action)

        if block.action is not Block.action:
          def nblock_action(dependencies, values, input, context):
            return block.action(dependencies, values, False, 1.0, input, context)
          self.__register_operator_callback__(block.nblock.name, "ADJ_NBLOCK_ACTION_CB", nblock_action)

    # Return the checkpoint flag
    return cs.value

  def __register_functional__(self, functional):
    assert(isinstance(functional, Functional))

    cfunc = functional.__cfunc_from_derivative__()
    self.functions_registered.append(cfunc)

    clib.adj_register_functional_derivative_callback(self.adjointer, str(functional), cfunc)

    cfunc = functional.__cfunc_from_functional__()
    self.functions_registered.append(cfunc)

    clib.adj_register_functional_callback(self.adjointer, str(functional), cfunc)

  def __register_parameter__(self, parameter):
    assert(isinstance(parameter, Parameter))

    cfunc = parameter.__cfunc_from_parameter_source__()
    self.functions_registered.append(cfunc)

    clib.adj_register_parameter_source_callback(self.adjointer, str(parameter), cfunc)

  def set_functional_dependencies(self, functional, timestep):

    dependencies = functional.dependencies(self,timestep)

    list_type = clib.adj_variable * len(dependencies)
    c_dependencies=list_type(*[var.var for var in dependencies])

    try:
      clib.adj_timestep_set_functional_dependencies(self.adjointer, timestep, str(functional), len(c_dependencies), c_dependencies)
    except exceptions.LibadjointErrorInvalidInputs:
      # Don't die in the case where these dependencies have already been set.
      pass

  def record_variable(self, var, storage):
    '''record_variable(self, var, storage)

    This method records the provided variable according to the settings in storage.'''

    storage_class = storage.vec.__class__
    raised_exception = False
    try:
      clib.adj_record_variable(self.adjointer, var.var, storage.storage_data)
    except exceptions.LibadjointWarnException, err:
      print err
      raised_exception = True

    # At this point we should also reregister the read and the delete callbacks.
    # Note that the initial callback implementation could not access
    # the user implementation of the read() and delete() functions, and but here we can.
    def __vec_read_callback__(adj_var, adj_vec_ptr):
        var = Variable(var=adj_var)
        y = storage_class.read(var)
        # Increase the reference counter of the new object to protect it from deallocation at the end of the callback
        references_taken.append(y)
        adj_vec_ptr[0].ptr = python_utils.c_ptr(y)
        adj_vec_ptr[0].klass = 0
        adj_vec_ptr[0].flags = 0

    def __vec_delete_callback__(adj_var):
        var = Variable(var=adj_var)
        storage_class.delete(var)

    self.__register_data_callback__('ADJ_VEC_DELETE_CB', __vec_delete_callback__)
    self.__register_data_callback__('ADJ_VEC_READ_CB', __vec_read_callback__)

    return not raised_exception

  def evaluate_functional(self, functional, timestep):
    '''evaluate_functional(self, functional, timestep)

    Evaluate the functional provided at t=timestep.'''

    self.__register_functional__(functional)
    self.set_functional_dependencies(functional, timestep)

    output=clib.c_double()

    clib.adj_evaluate_functional(self.adjointer, timestep, functional.__str__(), output)

    return output.value

  def to_html(self, filename, viztype):
    try:
      typecode = {"forward": 1, "adjoint": 2, "tlm": 3}[viztype]
    except IndexError:
      raise exceptions.LibadjointErrorInvalidInputs, "Argument viztype has to be one of the following: 'forward', 'adjoint', 'tlm'"

    if viztype == 'tlm':
      raise exceptions.LibadjointErrorNotImplemented, "HTML output for TLM is not implemented"

    clib.adj_adjointer_to_html(self.adjointer, filename, typecode)

  def get_forward_equation(self, equation):
    lhs = clib.adj_matrix()
    rhs = clib.adj_vector()
    fwd_var = clib.adj_variable()
    clib.adj_get_forward_equation(self.adjointer, equation, lhs, rhs, fwd_var)
    lhs_py = python_utils.c_deref(lhs.ptr)
    references_taken.remove(lhs_py)
    rhs_py = python_utils.c_deref(rhs.ptr)
    references_taken.remove(rhs_py)

    return (Variable(var=fwd_var), lhs_py, rhs_py)

  def get_forward_solution(self, equation):
    output = clib.adj_vector()
    fwd_var = clib.adj_variable()
    clib.adj_get_forward_solution(self.adjointer, equation, output, fwd_var)
    output_py = python_utils.c_deref(output.ptr)
    references_taken.remove(output_py)

    return (Variable(var=fwd_var), output_py)

  def get_adjoint_equation(self, equation, functional):

    self.__register_functional__(functional)
    self.set_functional_dependencies(functional, self.equation_timestep[equation])

    lhs = clib.adj_matrix()
    rhs = clib.adj_vector()
    adj_var = clib.adj_variable()
    clib.adj_get_adjoint_equation(self.adjointer, equation, str(functional), lhs, rhs, adj_var)
    lhs_py = python_utils.c_deref(lhs.ptr)
    references_taken.remove(lhs_py)
    rhs_py = python_utils.c_deref(rhs.ptr)
    references_taken.remove(rhs_py)

    return (Variable(var=adj_var), lhs_py, rhs_py)

  def get_adjoint_solution(self, equation, functional):

    self.__register_functional__(functional)
    self.set_functional_dependencies(functional, self.equation_timestep[equation])

    output = clib.adj_vector()
    adj_var = clib.adj_variable()
    clib.adj_get_adjoint_solution(self.adjointer, equation, str(functional), output, adj_var)
    output_py = python_utils.c_deref(output.ptr)
    references_taken.remove(output_py)

    return (Variable(var=adj_var), output_py)

  def get_tlm_equation(self, equation, parameter):

    self.__register_parameter__(parameter)
    #self.__set_parameter_dependencies__(parameter, self.equation_timestep[equation])

    lhs = clib.adj_matrix()
    rhs = clib.adj_vector()
    tlm_var = clib.adj_variable()
    clib.adj_get_tlm_equation(self.adjointer, equation, str(parameter), lhs, rhs, tlm_var)
    lhs_py = python_utils.c_deref(lhs.ptr)
    references_taken.remove(lhs_py)
    rhs_py = python_utils.c_deref(rhs.ptr)
    references_taken.remove(rhs_py)

    return (Variable(var=tlm_var), lhs_py, rhs_py)

  def get_tlm_solution(self, equation, parameter):

    self.__register_parameter__(parameter)
    #self.__set_parameter_dependencies__(parameter, self.equation_timestep[equation])

    output = clib.adj_vector()
    adj_var = clib.adj_variable()
    clib.adj_get_tlm_solution(self.adjointer, equation, str(parameter), output, adj_var)
    output_py = python_utils.c_deref(output.ptr)
    references_taken.remove(output_py)

    return (Variable(var=adj_var), output_py)

  def get_forward_variable(self, equation):
    fwd_var = clib.adj_variable()
    clib.adj_get_forward_variable(self.adjointer, equation, fwd_var)
    return Variable(var=fwd_var)

  def timestep_start_equation(self, timestep):
    timestep_start = ctypes.c_int()
    clib.adj_timestep_start_equation(self.adjointer, timestep, timestep_start)
    return timestep_start.value

  def timestep_end_equation(self, timestep):
    timestep_end = ctypes.c_int()
    clib.adj_timestep_end_equation(self.adjointer, timestep, timestep_end)
    return timestep_end.value

  def variable_known(self, variable):
    known = ctypes.c_int()
    clib.adj_variable_known(self.adjointer, variable.var, known)
    return (known.value == 1)

  def forget_adjoint_equation(self, equation):
    clib.adj_forget_adjoint_equation(self.adjointer, equation)

  def forget_forward_equation(self, equation):
    clib.adj_forget_forward_equation(self.adjointer, equation)

  def forget_tlm_equation(self, equation):
    clib.adj_forget_tlm_equation(self.adjointer, equation)

  def forget_adjoint_values(self, equation):
    clib.adj_forget_adjoint_values(self.adjointer, equation)

  def forget_tlm_values(self, equation):
    clib.adj_forget_tlm_values(self.adjointer, equation)

  def set_times(self, timestep, start, end):
    clib.adj_timestep_set_times(self.adjointer, timestep, start, end)

  def get_times(self, timestep):
    start = adj_scalar()
    end = adj_scalar()

    clib.adj_timestep_get_times(self.adjointer, timestep, start, end)
    return (start.value, end.value)

  def __register_data_callbacks__(self):
    self.__register_data_callback__('ADJ_VEC_DUPLICATE_CB', self.__vec_duplicate_callback__)
    self.__register_data_callback__('ADJ_VEC_DESTROY_CB', self.__vec_destroy_callback__)
    self.__register_data_callback__('ADJ_VEC_AXPY_CB', self.__vec_axpy_callback__)
    self.__register_data_callback__('ADJ_VEC_GET_NORM_CB', self.__vec_norm_callback__)
    self.__register_data_callback__('ADJ_VEC_DOT_PRODUCT_CB', self.__vec_dot_callback__)
    self.__register_data_callback__('ADJ_VEC_SET_RANDOM_CB', self.__vec_set_random_callback__)
    self.__register_data_callback__('ADJ_VEC_SET_VALUES_CB', self.__vec_set_values_callback__)
    self.__register_data_callback__('ADJ_VEC_GET_VALUES_CB', self.__vec_get_values_callback__)
    self.__register_data_callback__('ADJ_VEC_GET_SIZE_CB', self.__vec_get_size_callback__)
    self.__register_data_callback__('ADJ_VEC_WRITE_CB', self.__vec_write_callback__)
    self.__register_data_callback__('ADJ_VEC_READ_CB', self.__vec_read_callback__)
    self.__register_data_callback__('ADJ_VEC_DELETE_CB', self.__vec_delete_callback__)
    self.__register_data_callback__('ADJ_MAT_DUPLICATE_CB', self.__mat_duplicate_callback__)
    self.__register_data_callback__('ADJ_MAT_DESTROY_CB', self.__mat_destroy_callback__)
    self.__register_data_callback__('ADJ_MAT_ACTION_CB', self.__mat_action_callback__)
    self.__register_data_callback__('ADJ_MAT_AXPY_CB', self.__mat_axpy_callback__)
    self.__register_data_callback__('ADJ_SOLVE_CB', self.__mat_solve_callback__)

  def __register_data_callback__(self, type_name, func):
    type_id = int(constants.adj_constants[type_name])

    def newfunc(*args, **kwargs):
      try:
        func(*args, **kwargs)
      except:
        import sys
        import traceback
        print
        print "Python traceback: "
        traceback.print_exc()

        print

        # Try to print out a C traceback, too
        import ctypes
        try:
          libc = ctypes.CDLL("libc.so.6")
          datatype = ctypes.c_void_p * 200
          pointers = datatype()
          size = libc.backtrace(pointers, 200)
          print "C traceback: "
          libc.backtrace_symbols_fd(pointers, size, 2)
        except (OSError, AttributeError):
          pass

        sys.exit(1)

    type_to_api = {"ADJ_VEC_DESTROY_CB": self.vec_destroy_type,
                   "ADJ_VEC_DUPLICATE_CB": self.vec_duplicate_type,
                   "ADJ_VEC_AXPY_CB": self.vec_axpy_type,
                   'ADJ_VEC_GET_NORM_CB': self.vec_get_norm_type,
                   'ADJ_VEC_DOT_PRODUCT_CB': self.vec_dot_product_type,
                   'ADJ_VEC_SET_RANDOM_CB': self.vec_set_random_type,
                   'ADJ_VEC_SET_VALUES_CB': self.vec_set_values_type,
                   'ADJ_VEC_GET_VALUES_CB': self.vec_get_values_type,
                   'ADJ_VEC_GET_SIZE_CB': self.vec_get_size_type,
                   'ADJ_VEC_WRITE_CB': self.vec_write_type,
                   'ADJ_VEC_READ_CB': self.vec_read_type,
                   "ADJ_VEC_DELETE_CB": self.vec_delete_type,
                   "ADJ_MAT_DUPLICATE_CB": self.mat_duplicate_type,
                   "ADJ_MAT_DESTROY_CB": self.mat_destroy_type,
                   "ADJ_MAT_ACTION_CB": self.mat_action_type,
                   "ADJ_MAT_AXPY_CB": self.mat_axpy_type,
                   "ADJ_SOLVE_CB": self.solve_type}
    if type_name in type_to_api:
      clib.adj_register_data_callback.argtypes = [ctypes.POINTER(clib.adj_adjointer), ctypes.c_int, type_to_api[type_name]]
      data_function = type_to_api[type_name](newfunc)
      self.functions_registered.append(data_function)
      clib.adj_register_data_callback(self.adjointer, ctypes.c_int(type_id), data_function)

      clib.adj_register_data_callback.argtypes = [ctypes.POINTER(clib.adj_adjointer), ctypes.c_int, ctypes.CFUNCTYPE(None)]
    else:
      raise exceptions.LibadjointErrorNotImplemented("Unknown API for data callback " + type_name)

  def __cfunc_from_block_assembly__(self, bassembly_cb):
    '''Given a block assembly function defined using the Pythonic interface, we want to translate that into a function
    that can be called from C. This routine does exactly that.'''

    def cfunc(ndepends_c, variables_c, dependencies_c, hermitian_c, coefficient_c, context_c, output_c, rhs_c):
      # build the Python objects from the C objects
      variables = [Variable(var=variables_c[i]) for i in range(ndepends_c)]
      dependencies = [vector(dependencies_c[i]) for i in range(ndepends_c)]
      hermitian = (hermitian_c == 1)
      coefficient = coefficient_c
      context = context_c

      # Now call the callback we've been given
      (matrix, rhs) = bassembly_cb(variables, dependencies, hermitian, coefficient, context)

      # Now cast the outputs back to C
      output_c[0].klass = 0
      output_c[0].flags = 0
      output_c[0].ptr = 0
      if not isinstance(matrix, Matrix):
        raise exceptions.LibadjointErrorInvalidInputs("matrix object returned from block assembly callback must be a subclass of Matrix")
      output_c[0].ptr = python_utils.c_ptr(matrix)
      references_taken.append(matrix)

      rhs_c[0].klass = 0
      rhs_c[0].flags = 0
      rhs_c[0].ptr = 0
      if not isinstance(rhs, Vector):
        raise exceptions.LibadjointErrorInvalidInputs("rhs object returned from block assembly callback must be a subclass of Vector")
      rhs_c[0].ptr = python_utils.c_ptr(rhs)
      references_taken.append(rhs)

    return self.block_assembly_type(cfunc)

  def __cfunc_from_block_action__(self, baction_cb):
    '''Given a block action function defined using the Pythonic interface, we want to translate that into a function
    that can be called from C. This routine does exactly that.'''

    def cfunc(ndepends_c, variables_c, dependencies_c, hermitian_c, coefficient_c, input_c, context_c, output_c):
      # build the Python objects from the C objects
      variables = [Variable(var=variables_c[i]) for i in range(ndepends_c)]
      dependencies = [vector(dependencies_c[i]) for i in range(ndepends_c)]
      hermitian = (hermitian_c == 1)
      coefficient = coefficient_c
      input = vector(input_c)
      context = context_c

      # Now call the callback we've been given
      output = baction_cb(variables, dependencies, hermitian, coefficient, input, context)

      output_c[0].klass = 0
      output_c[0].flags = 0
      output_c[0].ptr = 0
      if not isinstance(output, Vector):
        raise exceptions.LibadjointErrorInvalidInputs("object returned from block action callback must be a subclass of Vector")
      output_c[0].ptr = python_utils.c_ptr(output)
      references_taken.append(output)


    return self.block_action_type(cfunc)

  def __cfunc_from_nblock_derivative_action__(self, nbaction_cb):
    '''Given a nonlinear block derivative action function defined using the Pythonic interface, we want to translate that into a function
    that can be called from C. This routine does exactly that.'''

    def cfunc(ndepends_c, dependencies_c, values_c, variable_c, contraction_c, hermitian_c, input_c, coefficient_c, context_c, output_c):
      # build the Python objects from the C objects
      dependencies = [Variable(var=dependencies_c[i]) for i in range(ndepends_c)]
      values = [vector(values_c[i]) for i in range(ndepends_c)]
      variable = Variable(var=variable_c)
      contraction = vector(contraction_c)
      hermitian = (hermitian_c == 1)
      input = vector(input_c)
      coefficient = coefficient_c
      context = context_c

      # Now call the callback we've been given
      output = nbaction_cb(dependencies, values, variable, contraction, hermitian, input, coefficient, context)

      output_c[0].klass = 0
      output_c[0].flags = 0
      output_c[0].ptr = 0
      if not isinstance(output, Vector):
        raise exceptions.LibadjointErrorInvalidInputs("object returned from nonlinear block derivative action callback must be a subclass of Vector")
      output_c[0].ptr = python_utils.c_ptr(output)
      references_taken.append(output)

    return self.nblock_derivative_action_type(cfunc)

  def __cfunc_from_nblock_action__(self, nbaction_cb):
    '''Given a nonlinear block action function defined using the Pythonic interface, we want to translate that into a function
    that can be called from C. This routine does exactly that.'''

    def cfunc(ndepends_c, dependencies_c, values_c, input_c, context_c, output_c):
      # build the Python objects from the C objects
      dependencies = [Variable(var=dependencies_c[i]) for i in range(ndepends_c)]
      values = [vector(values_c[i]) for i in range(ndepends_c)]
      input = vector(input_c)
      context = context_c

      # Now call the callback we've been given
      output = nbaction_cb(dependencies, values, input, context)

      output_c[0].klass = 0
      output_c[0].flags = 0
      output_c[0].ptr = 0
      if not isinstance(output, Vector):
        raise exceptions.LibadjointErrorInvalidInputs("object returned from nonlinear block derivative action callback must be a subclass of Vector")
      output_c[0].ptr = python_utils.c_ptr(output)
      references_taken.append(output)

    return self.nblock_action_type(cfunc)

  def __register_operator_callback__(self, name, type_name, func):
    type_to_api = {"ADJ_BLOCK_ASSEMBLY_CB": (self.block_assembly_type, self.__cfunc_from_block_assembly__),
                   "ADJ_BLOCK_ACTION_CB": (self.block_action_type, self.__cfunc_from_block_action__),
                   "ADJ_NBLOCK_DERIVATIVE_ACTION_CB": (self.nblock_derivative_action_type, self.__cfunc_from_nblock_derivative_action__),
                   "ADJ_NBLOCK_ACTION_CB": (self.nblock_action_type, self.__cfunc_from_nblock_action__)}
    if type_name in type_to_api:
      clib.adj_register_operator_callback.argtypes = [ctypes.POINTER(clib.adj_adjointer), ctypes.c_int, ctypes.c_char_p, type_to_api[type_name][0]]
      fn=type_to_api[type_name][1](func)
      self.functions_registered.append(fn)
      clib.adj_register_operator_callback(self.adjointer, int(constants.adj_constants[type_name]), name, fn)
      clib.adj_register_operator_callback.argtypes = [ctypes.POINTER(clib.adj_adjointer), ctypes.c_int, ctypes.c_char_p, ctypes.CFUNCTYPE(None)]
    else:
      raise exceptions.LibadjointErrorNotImplemented("Unknown API for data callback " + type_name)

  @staticmethod
  def __vec_duplicate_callback__(adj_vec, adj_vec_ptr):
    vec = vector(adj_vec)
    new_vec = vec.duplicate()

    # Increase the reference counter of the new object to protect it from deallocation at the end of the callback
    references_taken.append(new_vec)
    adj_vec_ptr[0].ptr = python_utils.c_ptr(new_vec)
    adj_vec_ptr[0].klass = 0
    adj_vec_ptr[0].flags = 0

  @staticmethod
  def __vec_destroy_callback__(adj_vec_ptr):
    vec = vector(adj_vec_ptr[0])

    # Do the corresponding decref of the object, so that the Python GC can pick it up
    try:
      references_taken.remove(vec)
    except:
      print vec.data
      print references_taken
      raise

  @staticmethod
  def __vec_axpy_callback__(adj_vec_ptr, alpha, adj_vec):
    y = vector(adj_vec_ptr[0])
    x = vector(adj_vec)
    assert isinstance(y, Vector)
    assert isinstance(x, Vector)
    y.axpy(alpha, x)

  @staticmethod
  def __vec_norm_callback__(adj_vec, x):
    y = vector(adj_vec)

    x[0] = y.norm()

  @staticmethod
  def __vec_dot_callback__(adj_vec_x, adj_vec_y, dot):
    x = vector(adj_vec_x)
    y = vector(adj_vec_y)

    dot[0] = x.dot_product(y)

  @staticmethod
  def __vec_set_random_callback__(adj_vec_ptr):
    y = vector(adj_vec_ptr[0])
    y.set_random()

  @staticmethod
  def __vec_set_values_callback__(adj_vec_ptr, values):
    import numpy
    y = vector(adj_vec_ptr[0])

    sz = y.size()
    nparray = numpy.zeros(sz)
    for i in range(sz):
      nparray[i] = values[i]
    y.set_values(numpy.array(nparray))

  @staticmethod
  def __vec_get_values_callback__(adj_vec, values_ptr):
    import numpy
    y = vector(adj_vec)

    sz = y.size()
    nparray = numpy.zeros(sz)
    y.get_values(nparray)

    for i in range(sz):
      values_ptr[0][i] = nparray[i]

  @staticmethod
  def __vec_get_size_callback__(adj_vec, sz):
    y = vector(adj_vec)
    sz[0] = y.size()

  @staticmethod
  def __vec_write_callback__(adj_var, adj_vec):
    var = Variable(var=adj_var)
    vec = vector(adj_vec)
    vec.write(var)

  @staticmethod
  def __vec_read_callback__(adj_var, adj_vec_ptr):
    raise exceptions.LibadjointErrorInvalidInputs(
        'Internal error: called vec_read callback before recording any variables.')

  @staticmethod
  def __vec_delete_callback__(adj_var):
    raise exceptions.LibadjointErrorInvalidInputs(
        'Internal error: called vec_delete callback before recording any variables.')

    print "To be implemented"
    assert(False)
    vec = vector(adj_vec_ptr[0])

    # Do the corresponding decref of the object, so that the Python GC can pick it up
    try:
      references_taken.remove(vec)
    except:
      print vec.data
      print references_taken
      raise

  @staticmethod
  def __mat_duplicate_callback__(adj_mat, adj_mat_ptr):
    mat = matrix(adj_mat)
    new_mat = mat.duplicate()

    # Increase the reference counter of the new object to protect it from deallocation at the end of the callback
    references_taken.append(new_mat)
    adj_mat_ptr.ptr = python_utils.c_ptr(new_mat)
    adj_mat_ptr.klass = 0
    adj_mat_ptr.flags = 0

  @staticmethod
  def __mat_destroy_callback__(adj_mat_ptr):
    mat = matrix(adj_mat_ptr[0])

    # Do the corresponding decref of the object, so that the Python GC can pick it up
    references_taken.remove(mat)

  @staticmethod
  def __mat_action_callback__(adj_mat, x_vec, y_ptr):
    y = vector(y_ptr[0])
    x = vector(x_vec)
    mat = matrix(adj_mat)
    mat.action(x, y)

  @staticmethod
  def __mat_axpy_callback__(adj_mat_ptr, alpha, adj_mat):
    y = matrix(adj_mat_ptr[0])
    x = matrix(adj_mat)
    y.axpy(alpha, x)

  @staticmethod
  def __mat_solve_callback__(adj_var, adj_mat, adj_rhs, adj_soln_ptr):
    A = matrix(adj_mat)
    b = vector(adj_rhs)

    x = A.solve(Variable(var=adj_var), b)
    references_taken.append(x)

    adj_soln_ptr[0].ptr = python_utils.c_ptr(x)
    adj_soln_ptr[0].klass = 0
    adj_soln_ptr[0].flags = 0

  def compute_gst(self, ic, ic_norm, final, final_norm, nrv):
    '''Computes the singular value decomposition of the propagator,
    possibly with some norms for the initial and final condition spaces.
    Pass ic_norm=None and final_norm=None to use the vector l2 norm.

    The propagator is the operator that maps
    (perturbations in the initial condition)
    to
    (perturbations in the final state)
    in a linear manner. Essentially, the propagator is the inverse
    of the tangent linear model.

    The singular value decomposition of the propagator is the basic
    tool in generalised stability and predictability analysis; see
    ``Atmospheric Modelling, Data Assimilation and Predictability''
    by E. Kalnay, chapter 6.

    ic -- an adj_variable corresponding to the initial condition
    ic_norm -- an adj_matrix with a norm for the initial condition.
               must be symmetric positive-definite
    final -- an adj_variable corresponding to the final condition
    final_norm -- an adj_matrix with a norm for the final condition.
                  must be symmetric positive-definite
    nrv -- number of requested singular vectors.'''

    handle = clib.adj_gst()
    ncv = ctypes.c_int()

    # Some hideous hackery to keep references around
    orig_final_norm = final_norm
    if final_norm is not None:
      final_norm = final_norm.as_adj_matrix()

    orig_ic_norm = ic_norm
    if ic_norm is not None:
      ic_norm = ic_norm.as_adj_matrix()

    clib.adj_compute_gst(self.adjointer, ic.var, ic_norm, final.var, final_norm, nrv, handle, ncv)

    gst = GSTHandle(handle, ncv)

    # We need to keep a reference to the adj_matrix, in case it gets deallocated
    gst.final_norm = final_norm
    gst.orig_final_norm = orig_final_norm
    gst.ic_norm = ic_norm
    gst.orig_ic_norm = orig_ic_norm
    return gst

  def get_variable_value(self, var):
    vec = clib.adj_vector()
    clib.adj_get_variable_value(self.adjointer, var.var, vec)
    return vector(vec)

class LinAlg(object):
  '''Base class for adjoint vector or matrix objects. In libadjoint,
  the operations performed on these are quite similar, so the common ones
  are factored out here. User applications should subclass Vector and Matrix
  instead of this directly.'''

  def __init__(self):
    raise exceptions.LibadjointErrorNotImplemented("Shouldn't ever instantiate a LinAlg object directly")

  def axpy(self, alpha, x):
    '''axpy(self, alpha, x)

    This method must update the Vector with self=self+alpha*x where
    alpha is a scalar and x is a Vector'''

    raise exceptions.LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no axpy(alpha,x) method')

class Vector(LinAlg):
  '''Base class for adjoint vector objects. User applications should
  subclass this and provide their own data and methods.'''

  def duplicate(self):
    '''duplicate(self)

    This method must return a newly allocated duplicate of its
    parent. The value of every entry of the duplicate must be zero.'''
    raise exceptions.LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no copy() method')

  def set_values(self, scalars):
    '''set_values(self, scalars)

    This method must set the value of Vector to that given by the array
    of scalars.'''

    raise exceptions.LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no set_values(scalars) method')

  def set_values(self, scalars):
    '''set_values(self, scalars)

    This method must set the value of scalars to that given by the local degrees of freedom
    of the Vector.'''

    raise exceptions.LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no get_values(scalars) method')

  def size(self):
    '''size(self)

    This method must return the number of degrees of freedom in this Vector.'''

    raise exceptions.LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no size() method')    

  def norm(self):
    '''norm(self)

    This method must return a norm for this vector. It does not matter which
    norm is chosen, so long as it satisfies the usual axioms for a norm.'''

    raise exceptions.LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no norm() method')    

  def set_random(self):
    '''This method must set the entries of a given vector x to pseudo-random values.'''

    raise exceptions.LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no set_random() method')        

  def dot_product(self, b):
    '''This method must return the result of dot(self, b).'''

    raise exceptions.LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no dot_product() method')        

  @staticmethod
  def read(var):
    '''This method must return a vector containing the values of var from disk.'''

    raise exceptions.LibadjointErrorNeedCallback(
      'Class Vector has no read() method')        

  @staticmethod
  def delete(var):
    '''This method must delete the vector containing the values of var from disk.''' 

    raise exceptions.LibadjointErrorNeedCallback(
      'Class Vector has no delete() method')        

  def write(self, var):
    '''This method must write this vector as values of var to disk.'''

    raise exceptions.LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no write() method')        

  def as_adj_vector(self):
    '''as_adj_vector(self)

    Returns an adj_vector with this Vector as its data payload.'''

    adj_vec = clib.adj_vector(ptr=python_utils.c_ptr(self))
    return adj_vec

class Matrix(LinAlg):
  '''Base class for adjoint matrix objects.'''

  def solve(self, var, b):
    '''solve(self, var, b)

    This method must solve the system self*x = b and return the answer. The adj_variable corresponding to the variable
    to be solved for is given by var.'''

    raise exceptions.LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no solve() method')        

  def action(self, x, y):
    '''action(self, x, y)

    This method must compute the action y = self*x.'''

    raise exceptions.LibadjointErrorNeedCallback(
      'Class '+self.__class__.__name__+' has no action() method')        

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

class GSTHandle(object):
  '''An object that wraps the result of a generalised stability theory computation.
     Request the computed results with get_gst.'''
  def __init__(self, handle, ncv):
    self.handle = handle
    self.ncv = ncv.value
    self.allocated = True

  def __del__(self):
    self.destroy()

  def get_gst(self, i, return_vectors=False, return_residual=False):
    if return_vectors:
      u = clib.adj_vector()
      v = clib.adj_vector()
    else:
      u = None
      v = None

    if return_residual:
      error = adj_scalar()
    else:
      error = None

    sigma = adj_scalar()

    clib.adj_get_gst(self.handle, i, sigma, u, v, error)

    retval = [sigma.value]
    if return_vectors:
      u_vec = vector(u)
      v_vec = vector(v)
      references_taken.remove(u_vec)
      references_taken.remove(v_vec)
      retval += [u_vec, v_vec]
    if return_residual:
      retval += [error.value]

    if len(retval) == 1:
      retval = retval[0]

    return retval

  def destroy(self):
    if self.allocated:
      clib.adj_destroy_gst(self.handle)
      self.allocated = False
