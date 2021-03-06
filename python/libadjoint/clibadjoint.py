from ctypes import *

import six

# _library = CDLL("/Users/minrk/conda/envs/fenics-py2/lib/libadjoint.dylib")
# so = 'dylib' if darwin else 'so'
from os.path import abspath, dirname, join
import sys

# prefix/lib/python2.7/site-packages/libadjoint/
_install_prefix = abspath(join(dirname(__file__), '..', '..', '..', '..'))

so = 'dylib' if sys.platform == 'darwin' else 'so'
# first, try loading abspath with prefix:
_libname = 'libadjoint.' + so
_library = None
for prefix in (_install_prefix, sys.prefix):
    try:
        _library = CDLL(join(prefix, 'lib', _libname))
    except OSError:
        pass

# Finally, fallback on default loader path:
if _library is None:
    _library = CDLL(_libname)

class STRING(object):
    @classmethod
    def from_param(cls, value):
        if isinstance(value, six.text_type):
            value = value.encode('utf8')
        return value

STRING_POINTER = POINTER(c_char_p)

class adj_adjointer(Structure):
    pass
adj_reset_revolve = _library.adj_reset_revolve
adj_reset_revolve.restype = c_int
adj_reset_revolve.argtypes = [POINTER(adj_adjointer)]
adj_advance_to_adjoint_run_revolve = _library.adj_advance_to_adjoint_run_revolve
adj_advance_to_adjoint_run_revolve.restype = c_int
adj_advance_to_adjoint_run_revolve.argtypes = [POINTER(adj_adjointer)]
adj_create_adjointer = _library.adj_create_adjointer
adj_create_adjointer.restype = c_int
adj_create_adjointer.argtypes = [POINTER(adj_adjointer)]
adj_destroy_adjointer = _library.adj_destroy_adjointer
adj_destroy_adjointer.restype = c_int
adj_destroy_adjointer.argtypes = [POINTER(adj_adjointer)]
adj_deactivate_adjointer = _library.adj_deactivate_adjointer
adj_deactivate_adjointer.restype = c_int
adj_deactivate_adjointer.argtypes = [POINTER(adj_adjointer)]
adj_get_checkpoint_strategy = _library.adj_get_checkpoint_strategy
adj_get_checkpoint_strategy.restype = c_int
adj_get_checkpoint_strategy.argtypes = [POINTER(adj_adjointer), POINTER(c_int)]
adj_set_checkpoint_strategy = _library.adj_set_checkpoint_strategy
adj_set_checkpoint_strategy.restype = c_int
adj_set_checkpoint_strategy.argtypes = [POINTER(adj_adjointer), c_int]
adj_set_revolve_options = _library.adj_set_revolve_options
adj_set_revolve_options.restype = c_int
adj_set_revolve_options.argtypes = [POINTER(adj_adjointer), c_int, c_int, c_int, c_int]
adj_set_revolve_debug_options = _library.adj_set_revolve_debug_options
adj_set_revolve_debug_options.restype = c_int
adj_set_revolve_debug_options.argtypes = [POINTER(adj_adjointer), c_int, c_double]
adj_equation_count = _library.adj_equation_count
adj_equation_count.restype = c_int
adj_equation_count.argtypes = [POINTER(adj_adjointer), POINTER(c_int)]
class adj_equation(Structure):
    pass
class adj_variable(Structure):
    pass
adj_variable._fields_ = [
    ('name', c_char * 4080),
    ('timestep', c_int),
    ('iteration', c_int),
    ('type', c_int),
    ('auxiliary', c_int),
    ('functional', c_char * 4080),
]
class adj_block(Structure):
    pass
class adj_vector(Structure):
    pass
adj_vector._fields_ = [
    ('ptr', c_void_p),
    ('klass', c_int),
    ('flags', c_int),
]
class adj_matrix(Structure):
    pass
adj_equation._fields_ = [
    ('variable', adj_variable),
    ('nblocks', c_int),
    ('blocks', POINTER(adj_block)),
    ('targets', POINTER(adj_variable)),
    ('nrhsdeps', c_int),
    ('rhsdeps', POINTER(adj_variable)),
    ('rhs_context', c_void_p),
    ('rhs_callback', CFUNCTYPE(None, c_void_p, adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), c_void_p, POINTER(adj_vector), POINTER(c_int))),
    ('rhs_deriv_action_callback', CFUNCTYPE(None, c_void_p, adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), adj_variable, adj_vector, c_int, c_void_p, POINTER(adj_vector), POINTER(c_int))),
    ('rhs_second_deriv_action_callback', CFUNCTYPE(None, c_void_p, adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), adj_variable, adj_vector, adj_variable, c_int, adj_vector, c_void_p, POINTER(adj_vector), POINTER(c_int))),
    ('rhs_deriv_assembly_callback', CFUNCTYPE(None, c_void_p, adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), c_int, c_void_p, POINTER(adj_matrix))),
    ('memory_checkpoint', c_int),
    ('disk_checkpoint', c_int),
]
adj_register_equation = _library.adj_register_equation
adj_register_equation.restype = c_int
adj_register_equation.argtypes = [POINTER(adj_adjointer), adj_equation, POINTER(c_int)]
class adj_storage_data(Structure):
    pass
adj_storage_data._fields_ = [
    ('compare', c_int),
    ('comparison_tolerance', c_double),
    ('overwrite', c_int),
    ('value', adj_vector),
    ('storage_memory_type', c_int),
    ('storage_memory_has_value', c_int),
    ('storage_memory_is_checkpoint', c_int),
    ('storage_disk_has_value', c_int),
    ('storage_disk_is_checkpoint', c_int),
]
adj_record_variable = _library.adj_record_variable
adj_record_variable.restype = c_int
adj_record_variable.argtypes = [POINTER(adj_adjointer), adj_variable, adj_storage_data]
adj_register_operator_callback = _library.adj_register_operator_callback
adj_register_operator_callback.restype = c_int
adj_register_operator_callback.argtypes = [POINTER(adj_adjointer), c_int, STRING, CFUNCTYPE(None)]
adj_register_data_callback = _library.adj_register_data_callback
adj_register_data_callback.restype = c_int
adj_register_data_callback.argtypes = [POINTER(adj_adjointer), c_int, CFUNCTYPE(None)]
adj_register_functional_callback = _library.adj_register_functional_callback
adj_register_functional_callback.restype = c_int
adj_register_functional_callback.argtypes = [POINTER(adj_adjointer), STRING, CFUNCTYPE(None, POINTER(adj_adjointer), c_int, c_int, POINTER(adj_variable), POINTER(adj_vector), c_char_p, POINTER(c_double))]
adj_register_functional_derivative_callback = _library.adj_register_functional_derivative_callback
adj_register_functional_derivative_callback.restype = c_int
adj_register_functional_derivative_callback.argtypes = [POINTER(adj_adjointer), STRING, CFUNCTYPE(None, POINTER(adj_adjointer), adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), c_char_p, POINTER(adj_vector))]
adj_register_functional_second_derivative_callback = _library.adj_register_functional_second_derivative_callback
adj_register_functional_second_derivative_callback.restype = c_int
adj_register_functional_second_derivative_callback.argtypes = [POINTER(adj_adjointer), STRING, CFUNCTYPE(None, POINTER(adj_adjointer), adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), adj_vector, c_char_p, POINTER(adj_vector))]
adj_register_parameter_source_callback = _library.adj_register_parameter_source_callback
adj_register_parameter_source_callback.restype = c_int
adj_register_parameter_source_callback.argtypes = [POINTER(adj_adjointer), STRING, CFUNCTYPE(None, POINTER(adj_adjointer), c_int, adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), c_char_p, POINTER(adj_vector), POINTER(c_int))]
adj_forget_adjoint_equation = _library.adj_forget_adjoint_equation
adj_forget_adjoint_equation.restype = c_int
adj_forget_adjoint_equation.argtypes = [POINTER(adj_adjointer), c_int]
adj_forget_forward_equation = _library.adj_forget_forward_equation
adj_forget_forward_equation.restype = c_int
adj_forget_forward_equation.argtypes = [POINTER(adj_adjointer), c_int]
adj_forget_tlm_equation = _library.adj_forget_tlm_equation
adj_forget_tlm_equation.restype = c_int
adj_forget_tlm_equation.argtypes = [POINTER(adj_adjointer), c_int]
adj_forget_adjoint_values = _library.adj_forget_adjoint_values
adj_forget_adjoint_values.restype = c_int
adj_forget_adjoint_values.argtypes = [POINTER(adj_adjointer), c_int]
adj_forget_tlm_values = _library.adj_forget_tlm_values
adj_forget_tlm_values.restype = c_int
adj_forget_tlm_values.argtypes = [POINTER(adj_adjointer), c_int]
adj_timestep_count = _library.adj_timestep_count
adj_timestep_count.restype = c_int
adj_timestep_count.argtypes = [POINTER(adj_adjointer), POINTER(c_int)]
adj_iteration_count = _library.adj_iteration_count
adj_iteration_count.restype = c_int
adj_iteration_count.argtypes = [POINTER(adj_adjointer), adj_variable, POINTER(c_int)]
adj_timestep_start_equation = _library.adj_timestep_start_equation
adj_timestep_start_equation.restype = c_int
adj_timestep_start_equation.argtypes = [POINTER(adj_adjointer), c_int, POINTER(c_int)]
adj_timestep_end_equation = _library.adj_timestep_end_equation
adj_timestep_end_equation.restype = c_int
adj_timestep_end_equation.argtypes = [POINTER(adj_adjointer), c_int, POINTER(c_int)]
adj_timestep_set_times = _library.adj_timestep_set_times
adj_timestep_set_times.restype = c_int
adj_timestep_set_times.argtypes = [POINTER(adj_adjointer), c_int, c_double, c_double]
adj_timestep_get_times = _library.adj_timestep_get_times
adj_timestep_get_times.restype = c_int
adj_timestep_get_times.argtypes = [POINTER(adj_adjointer), c_int, POINTER(c_double), POINTER(c_double)]
adj_timestep_set_functional_dependencies = _library.adj_timestep_set_functional_dependencies
adj_timestep_set_functional_dependencies.restype = c_int
adj_timestep_set_functional_dependencies.argtypes = [POINTER(adj_adjointer), c_int, STRING, c_int, POINTER(adj_variable)]
adj_storage_memory_copy = _library.adj_storage_memory_copy
adj_storage_memory_copy.restype = c_int
adj_storage_memory_copy.argtypes = [adj_vector, POINTER(adj_storage_data)]
adj_storage_memory_incref = _library.adj_storage_memory_incref
adj_storage_memory_incref.restype = c_int
adj_storage_memory_incref.argtypes = [adj_vector, POINTER(adj_storage_data)]
adj_storage_disk = _library.adj_storage_disk
adj_storage_disk.restype = c_int
adj_storage_disk.argtypes = [adj_vector, POINTER(adj_storage_data)]
adj_storage_set_compare = _library.adj_storage_set_compare
adj_storage_set_compare.restype = c_int
adj_storage_set_compare.argtypes = [POINTER(adj_storage_data), c_int, c_double]
adj_storage_set_overwrite = _library.adj_storage_set_overwrite
adj_storage_set_overwrite.restype = c_int
adj_storage_set_overwrite.argtypes = [POINTER(adj_storage_data), c_int]
adj_storage_set_checkpoint = _library.adj_storage_set_checkpoint
adj_storage_set_checkpoint.restype = c_int
adj_storage_set_checkpoint.argtypes = [POINTER(adj_storage_data), c_int]
adj_variable_known = _library.adj_variable_known
adj_variable_known.restype = c_int
adj_variable_known.argtypes = [POINTER(adj_adjointer), adj_variable, POINTER(c_int)]
adj_get_variable_value = _library.adj_get_variable_value
adj_get_variable_value.restype = c_int
adj_get_variable_value.argtypes = [POINTER(adj_adjointer), adj_variable, POINTER(adj_vector)]
adj_set_finished = _library.adj_set_finished
adj_set_finished.restype = c_int
adj_set_finished.argtypes = [POINTER(adj_adjointer), c_int]
adj_get_finished = _library.adj_get_finished
adj_get_finished.restype = c_int
adj_get_finished.argtypes = [POINTER(adj_adjointer), POINTER(c_int)]
adj_get_forward_variable = _library.adj_get_forward_variable
adj_get_forward_variable.restype = c_int
adj_get_forward_variable.argtypes = [POINTER(adj_adjointer), c_int, POINTER(adj_variable)]
adj_adjointer_to_html = _library.adj_adjointer_to_html
adj_adjointer_to_html.restype = c_int
adj_adjointer_to_html.argtypes = [POINTER(adj_adjointer), STRING, c_int]
adj_get_adjoint_equation = _library.adj_get_adjoint_equation
adj_get_adjoint_equation.restype = c_int
adj_get_adjoint_equation.argtypes = [POINTER(adj_adjointer), c_int, STRING, POINTER(adj_matrix), POINTER(adj_vector), POINTER(adj_variable)]
adj_get_adjoint_solution = _library.adj_get_adjoint_solution
adj_get_adjoint_solution.restype = c_int
adj_get_adjoint_solution.argtypes = [POINTER(adj_adjointer), c_int, STRING, POINTER(adj_vector), POINTER(adj_variable)]
adj_get_forward_equation = _library.adj_get_forward_equation
adj_get_forward_equation.restype = c_int
adj_get_forward_equation.argtypes = [POINTER(adj_adjointer), c_int, POINTER(adj_matrix), POINTER(adj_vector), POINTER(adj_variable)]
adj_get_forward_solution = _library.adj_get_forward_solution
adj_get_forward_solution.restype = c_int
adj_get_forward_solution.argtypes = [POINTER(adj_adjointer), c_int, POINTER(adj_vector), POINTER(adj_variable)]
adj_get_tlm_equation = _library.adj_get_tlm_equation
adj_get_tlm_equation.restype = c_int
adj_get_tlm_equation.argtypes = [POINTER(adj_adjointer), c_int, STRING, POINTER(adj_matrix), POINTER(adj_vector), POINTER(adj_variable)]
adj_get_tlm_solution = _library.adj_get_tlm_solution
adj_get_tlm_solution.restype = c_int
adj_get_tlm_solution.argtypes = [POINTER(adj_adjointer), c_int, STRING, POINTER(adj_vector), POINTER(adj_variable)]
adj_get_soa_equation = _library.adj_get_soa_equation
adj_get_soa_equation.restype = c_int
adj_get_soa_equation.argtypes = [POINTER(adj_adjointer), c_int, STRING, STRING, POINTER(adj_matrix), POINTER(adj_vector), POINTER(adj_variable)]
adj_get_soa_solution = _library.adj_get_soa_solution
adj_get_soa_solution.restype = c_int
adj_get_soa_solution.argtypes = [POINTER(adj_adjointer), c_int, STRING, STRING, POINTER(adj_vector), POINTER(adj_variable)]
adj_matrix._fields_ = [
    ('ptr', c_void_p),
    ('klass', c_int),
    ('flags', c_int),
]
class adj_nonlinear_block(Structure):
    pass
adj_nonlinear_block._fields_ = [
    ('name', c_char * 4080),
    ('coefficient', c_double),
    ('ndepends', c_int),
    ('depends', POINTER(adj_variable)),
    ('context', c_void_p),
    ('test_deriv_hermitian', c_int),
    ('number_of_tests', c_int),
    ('tolerance', c_double),
    ('test_derivative', c_int),
    ('number_of_rounds', c_int),
]
adj_block._fields_ = [
    ('name', c_char * 4080),
    ('has_nonlinear_block', c_int),
    ('nonlinear_block', adj_nonlinear_block),
    ('context', c_void_p),
    ('hermitian', c_int),
    ('coefficient', c_double),
    ('test_hermitian', c_int),
    ('number_of_tests', c_int),
    ('tolerance', c_double),
]
class adj_term(Structure):
    pass
adj_term._fields_ = [
    ('nblocks', c_int),
    ('blocks', POINTER(adj_block)),
    ('targets', POINTER(adj_variable)),
]
class adj_variable_data(Structure):
    pass
adj_variable_data._fields_ = [
    ('equation', c_int),
    ('type', c_int),
    ('ntargeting_equations', c_int),
    ('targeting_equations', POINTER(c_int)),
    ('ndepending_equations', c_int),
    ('depending_equations', POINTER(c_int)),
    ('nrhs_equations', c_int),
    ('rhs_equations', POINTER(c_int)),
    ('ndepending_timesteps', c_int),
    ('depending_timesteps', POINTER(c_int)),
    ('nadjoint_equations', c_int),
    ('adjoint_equations', POINTER(c_int)),
    ('storage', adj_storage_data),
    ('next', POINTER(adj_variable_data)),
]
class adj_data_callbacks(Structure):
    pass
adj_data_callbacks._fields_ = [
    ('vec_duplicate', CFUNCTYPE(None, adj_vector, POINTER(adj_vector))),
    ('vec_axpy', CFUNCTYPE(None, POINTER(adj_vector), c_double, adj_vector)),
    ('vec_destroy', CFUNCTYPE(None, POINTER(adj_vector))),
    ('vec_set_values', CFUNCTYPE(None, POINTER(adj_vector), POINTER(c_double))),
    ('vec_get_values', CFUNCTYPE(None, adj_vector, POINTER(POINTER(c_double)))),
    ('vec_get_size', CFUNCTYPE(None, adj_vector, POINTER(c_int))),
    ('vec_divide', CFUNCTYPE(None, POINTER(adj_vector), adj_vector)),
    ('vec_get_norm', CFUNCTYPE(None, adj_vector, POINTER(c_double))),
    ('vec_dot_product', CFUNCTYPE(None, adj_vector, adj_vector, POINTER(c_double))),
    ('vec_set_random', CFUNCTYPE(None, POINTER(adj_vector))),
    ('vec_write', CFUNCTYPE(None, adj_variable, adj_vector)),
    ('vec_read', CFUNCTYPE(None, adj_variable, POINTER(adj_vector))),
    ('vec_delete', CFUNCTYPE(None, adj_variable)),
    ('mat_duplicate', CFUNCTYPE(None, adj_matrix, POINTER(adj_matrix))),
    ('mat_axpy', CFUNCTYPE(None, POINTER(adj_matrix), c_double, adj_matrix)),
    ('mat_destroy', CFUNCTYPE(None, POINTER(adj_matrix))),
    ('mat_action', CFUNCTYPE(None, adj_matrix, adj_vector, POINTER(adj_vector))),
    ('solve', CFUNCTYPE(None, adj_variable, adj_matrix, adj_vector, POINTER(adj_vector))),
]
class adj_op_callback(Structure):
    pass
adj_op_callback._fields_ = [
    ('name', c_char * 4080),
    ('callback', CFUNCTYPE(None)),
    ('next', POINTER(adj_op_callback)),
]
class adj_op_callback_list(Structure):
    pass
adj_op_callback_list._fields_ = [
    ('firstnode', POINTER(adj_op_callback)),
    ('lastnode', POINTER(adj_op_callback)),
]
class adj_func_callback(Structure):
    pass
adj_func_callback._fields_ = [
    ('name', c_char * 4080),
    ('callback', CFUNCTYPE(None, c_void_p, c_int, c_int, POINTER(adj_variable), POINTER(adj_vector), c_char_p, POINTER(c_double))),
    ('next', POINTER(adj_func_callback)),
]
class adj_func_deriv_callback(Structure):
    pass
adj_func_deriv_callback._fields_ = [
    ('name', c_char * 4080),
    ('callback', CFUNCTYPE(None, c_void_p, adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), c_char_p, POINTER(adj_vector))),
    ('next', POINTER(adj_func_deriv_callback)),
]
class adj_func_second_deriv_callback(Structure):
    pass
adj_func_second_deriv_callback._fields_ = [
    ('name', c_char * 4080),
    ('callback', CFUNCTYPE(None, c_void_p, adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), adj_vector, c_char_p, POINTER(adj_vector))),
    ('next', POINTER(adj_func_second_deriv_callback)),
]
class adj_parameter_source_callback(Structure):
    pass
adj_parameter_source_callback._fields_ = [
    ('name', c_char * 4080),
    ('callback', CFUNCTYPE(None, c_void_p, c_int, adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), c_char_p, POINTER(adj_vector), POINTER(c_int))),
    ('next', POINTER(adj_parameter_source_callback)),
]
class adj_func_callback_list(Structure):
    pass
adj_func_callback_list._fields_ = [
    ('firstnode', POINTER(adj_func_callback)),
    ('lastnode', POINTER(adj_func_callback)),
]
class adj_func_deriv_callback_list(Structure):
    pass
adj_func_deriv_callback_list._fields_ = [
    ('firstnode', POINTER(adj_func_deriv_callback)),
    ('lastnode', POINTER(adj_func_deriv_callback)),
]
class adj_func_second_deriv_callback_list(Structure):
    pass
adj_func_second_deriv_callback_list._fields_ = [
    ('firstnode', POINTER(adj_func_second_deriv_callback)),
    ('lastnode', POINTER(adj_func_second_deriv_callback)),
]
class adj_parameter_source_callback_list(Structure):
    pass
adj_parameter_source_callback_list._fields_ = [
    ('firstnode', POINTER(adj_parameter_source_callback)),
    ('lastnode', POINTER(adj_parameter_source_callback)),
]
class adj_nonlinear_block_derivative(Structure):
    pass
adj_nonlinear_block_derivative._fields_ = [
    ('nonlinear_block', adj_nonlinear_block),
    ('variable', adj_variable),
    ('contraction', adj_vector),
    ('hermitian', c_int),
    ('outer', c_int),
]
class adj_nonlinear_block_second_derivative(Structure):
    pass
adj_nonlinear_block_second_derivative._fields_ = [
    ('nonlinear_block', adj_nonlinear_block),
    ('inner_variable', adj_variable),
    ('inner_contraction', adj_vector),
    ('outer_variable', adj_variable),
    ('outer_contraction', adj_vector),
    ('hermitian', c_int),
    ('block_action', adj_vector),
]
class adj_variable_hash(Structure):
    pass
class adj_hash_handle(Structure):
    pass
class UT_hash_table(Structure):
    pass
adj_hash_handle._fields_ = [
    ('tbl', POINTER(UT_hash_table)),
    ('prev', c_void_p),
    ('next', c_void_p),
    ('hh_prev', POINTER(adj_hash_handle)),
    ('hh_next', POINTER(adj_hash_handle)),
    ('key', c_void_p),
    ('keylen', c_uint),
    ('hashv', c_uint),
]
adj_variable_hash._fields_ = [
    ('variable', adj_variable),
    ('data', POINTER(adj_variable_data)),
    ('hh', adj_hash_handle),
]
class adj_dictionary_entry(Structure):
    pass
adj_dictionary_entry._fields_ = [
    ('key', c_char * 32768),
    ('value', c_char * 32768),
    ('hh', adj_hash_handle),
]
class adj_dictionary(Structure):
    pass
adj_dictionary._fields_ = [
    ('dict', POINTER(adj_dictionary_entry)),
]
class adj_functional_data(Structure):
    pass
adj_functional_data._fields_ = [
    ('name', c_char * 4080),
    ('ndepends', c_int),
    ('dependencies', POINTER(adj_variable)),
    ('next', POINTER(adj_functional_data)),
]
class adj_timestep_data(Structure):
    pass
adj_timestep_data._fields_ = [
    ('start_equation', c_int),
    ('start_time', c_double),
    ('end_time', c_double),
    ('functional_data_start', POINTER(adj_functional_data)),
    ('functional_data_end', POINTER(adj_functional_data)),
]
class adj_revolve_data(Structure):
    pass
class CRevolve(Structure):
    pass
CRevolve._fields_ = [
    ('ptr', c_void_p),
]

# values for enumeration 'CACTION'
CACTION_ADVANCE = 0
CACTION_TAKESHOT = 1
CACTION_RESTORE = 2
CACTION_FIRSTRUN = 3
CACTION_YOUTURN = 4
CACTION_TERMINATE = 5
CACTION_ERROR = 6
CACTION = c_int # enum
adj_revolve_data._fields_ = [
    ('revolve', CRevolve),
    ('snaps', c_int),
    ('snaps_in_ram', c_int),
    ('steps', c_int),
    ('current_action', CACTION),
    ('current_timestep', c_int),
    ('verbose', c_int),
    ('overwrite', c_int),
    ('comparison_tolerance', c_double),
]
adj_adjointer._fields_ = [
    ('equations', POINTER(adj_equation)),
    ('nequations', c_int),
    ('equations_sz', c_int),
    ('ntimesteps', c_int),
    ('timestep_data', POINTER(adj_timestep_data)),
    ('revolve_data', adj_revolve_data),
    ('varhash', POINTER(adj_variable_hash)),
    ('options', c_int * 3),
    ('callbacks', adj_data_callbacks),
    ('nonlinear_action_list', adj_op_callback_list),
    ('nonlinear_derivative_action_list', adj_op_callback_list),
    ('nonlinear_derivative_assembly_list', adj_op_callback_list),
    ('block_action_list', adj_op_callback_list),
    ('block_assembly_list', adj_op_callback_list),
    ('nonlinear_second_derivative_action_list', adj_op_callback_list),
    ('nonlinear_derivative_outer_action_list', adj_op_callback_list),
    ('functional_list', adj_func_callback_list),
    ('functional_derivative_list', adj_func_deriv_callback_list),
    ('functional_second_derivative_list', adj_func_second_deriv_callback_list),
    ('parameter_source_list', adj_parameter_source_callback_list),
    ('finished', c_int),
]
adj_create_variable = _library.adj_create_variable
adj_create_variable.restype = c_int
adj_create_variable.argtypes = [STRING, c_int, c_int, c_int, POINTER(adj_variable)]
adj_variable_get_name = _library.adj_variable_get_name
adj_variable_get_name.restype = c_int
adj_variable_get_name.argtypes = [adj_variable, STRING_POINTER]
adj_variable_get_timestep = _library.adj_variable_get_timestep
adj_variable_get_timestep.restype = c_int
adj_variable_get_timestep.argtypes = [adj_variable, POINTER(c_int)]
adj_variable_get_iteration = _library.adj_variable_get_iteration
adj_variable_get_iteration.restype = c_int
adj_variable_get_iteration.argtypes = [adj_variable, POINTER(c_int)]
adj_variable_get_type = _library.adj_variable_get_type
adj_variable_get_type.restype = c_int
adj_variable_get_type.argtypes = [adj_variable, POINTER(c_int)]
adj_variable_set_auxiliary = _library.adj_variable_set_auxiliary
adj_variable_set_auxiliary.restype = c_int
adj_variable_set_auxiliary.argtypes = [POINTER(adj_variable), c_int]
size_t = c_ulong
adj_variable_str = _library.adj_variable_str
adj_variable_str.restype = c_int
adj_variable_str.argtypes = [adj_variable, STRING, size_t]
adj_create_nonlinear_block = _library.adj_create_nonlinear_block
adj_create_nonlinear_block.restype = c_int
adj_create_nonlinear_block.argtypes = [STRING, c_int, POINTER(adj_variable), c_void_p, c_double, POINTER(adj_nonlinear_block)]
adj_destroy_nonlinear_block = _library.adj_destroy_nonlinear_block
adj_destroy_nonlinear_block.restype = c_int
adj_destroy_nonlinear_block.argtypes = [POINTER(adj_nonlinear_block)]
adj_nonlinear_block_set_coefficient = _library.adj_nonlinear_block_set_coefficient
adj_nonlinear_block_set_coefficient.restype = c_int
adj_nonlinear_block_set_coefficient.argtypes = [POINTER(adj_nonlinear_block), c_double]
adj_create_block = _library.adj_create_block
adj_create_block.restype = c_int
adj_create_block.argtypes = [STRING, POINTER(adj_nonlinear_block), c_void_p, c_double, POINTER(adj_block)]
adj_destroy_block = _library.adj_destroy_block
adj_destroy_block.restype = c_int
adj_destroy_block.argtypes = [POINTER(adj_block)]
adj_block_set_coefficient = _library.adj_block_set_coefficient
adj_block_set_coefficient.restype = c_int
adj_block_set_coefficient.argtypes = [POINTER(adj_block), c_double]
adj_block_set_hermitian = _library.adj_block_set_hermitian
adj_block_set_hermitian.restype = c_int
adj_block_set_hermitian.argtypes = [POINTER(adj_block), c_int]
adj_block_set_test_hermitian = _library.adj_block_set_test_hermitian
adj_block_set_test_hermitian.restype = c_int
adj_block_set_test_hermitian.argtypes = [POINTER(adj_block), c_int, c_int, c_double]
adj_nonlinear_block_set_test_hermitian = _library.adj_nonlinear_block_set_test_hermitian
adj_nonlinear_block_set_test_hermitian.restype = c_int
adj_nonlinear_block_set_test_hermitian.argtypes = [POINTER(adj_nonlinear_block), c_int, c_int, c_double]
adj_nonlinear_block_set_test_derivative = _library.adj_nonlinear_block_set_test_derivative
adj_nonlinear_block_set_test_derivative.restype = c_int
adj_nonlinear_block_set_test_derivative.argtypes = [POINTER(adj_nonlinear_block), c_int, c_int]
adj_create_equation = _library.adj_create_equation
adj_create_equation.restype = c_int
adj_create_equation.argtypes = [adj_variable, c_int, POINTER(adj_block), POINTER(adj_variable), POINTER(adj_equation)]
adj_equation_set_rhs_dependencies = _library.adj_equation_set_rhs_dependencies
adj_equation_set_rhs_dependencies.restype = c_int
adj_equation_set_rhs_dependencies.argtypes = [POINTER(adj_equation), c_int, POINTER(adj_variable), c_void_p]
adj_destroy_equation = _library.adj_destroy_equation
adj_destroy_equation.restype = c_int
adj_destroy_equation.argtypes = [POINTER(adj_equation)]
adj_create_term = _library.adj_create_term
adj_create_term.restype = c_int
adj_create_term.argtypes = [c_int, POINTER(adj_block), POINTER(adj_variable), POINTER(adj_term)]
adj_add_terms = _library.adj_add_terms
adj_add_terms.restype = c_int
adj_add_terms.argtypes = [adj_term, adj_term, POINTER(adj_term)]
adj_destroy_term = _library.adj_destroy_term
adj_destroy_term.restype = c_int
adj_destroy_term.argtypes = [POINTER(adj_term)]
adj_add_term_to_equation = _library.adj_add_term_to_equation
adj_add_term_to_equation.restype = c_int
adj_add_term_to_equation.argtypes = [adj_term, POINTER(adj_equation)]
adj_equation_set_rhs_callback = _library.adj_equation_set_rhs_callback
adj_equation_set_rhs_callback.restype = c_int
adj_equation_set_rhs_callback.argtypes = [POINTER(adj_equation), CFUNCTYPE(None, POINTER(adj_adjointer), adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), c_void_p, POINTER(adj_vector), POINTER(c_int))]
adj_equation_set_rhs_derivative_action_callback = _library.adj_equation_set_rhs_derivative_action_callback
adj_equation_set_rhs_derivative_action_callback.restype = c_int
adj_equation_set_rhs_derivative_action_callback.argtypes = [POINTER(adj_equation), CFUNCTYPE(None, POINTER(adj_adjointer), adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), adj_variable, adj_vector, c_int, c_void_p, POINTER(adj_vector), POINTER(c_int))]
adj_equation_set_rhs_second_derivative_action_callback = _library.adj_equation_set_rhs_second_derivative_action_callback
adj_equation_set_rhs_second_derivative_action_callback.restype = c_int
adj_equation_set_rhs_second_derivative_action_callback.argtypes = [POINTER(adj_equation), CFUNCTYPE(None, POINTER(adj_adjointer), adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), adj_variable, adj_vector, adj_variable, c_int, adj_vector, c_void_p, POINTER(adj_vector), POINTER(c_int))]
adj_equation_set_rhs_derivative_assembly_callback = _library.adj_equation_set_rhs_derivative_assembly_callback
adj_equation_set_rhs_derivative_assembly_callback.restype = c_int
adj_equation_set_rhs_derivative_assembly_callback.argtypes = [POINTER(adj_equation), CFUNCTYPE(None, POINTER(adj_adjointer), adj_variable, c_int, POINTER(adj_variable), POINTER(adj_vector), c_int, c_void_p, POINTER(adj_matrix))]
adj_variable_equal = _library.adj_variable_equal
adj_variable_equal.restype = c_int
adj_variable_equal.argtypes = [POINTER(adj_variable), POINTER(adj_variable), c_int]
adj_adjointer_check_consistency = _library.adj_adjointer_check_consistency
adj_adjointer_check_consistency.restype = c_int
adj_adjointer_check_consistency.argtypes = [POINTER(adj_adjointer)]
adj_adjointer_check_checkpoints = _library.adj_adjointer_check_checkpoints
adj_adjointer_check_checkpoints.restype = c_int
adj_adjointer_check_checkpoints.argtypes = [POINTER(adj_adjointer)]
adj_dict_init = _library.adj_dict_init
adj_dict_init.restype = c_int
adj_dict_init.argtypes = [POINTER(adj_dictionary)]
adj_dict_set = _library.adj_dict_set
adj_dict_set.restype = c_int
adj_dict_set.argtypes = [POINTER(adj_dictionary), STRING, STRING]
adj_dict_find = _library.adj_dict_find
adj_dict_find.restype = c_int
adj_dict_find.argtypes = [POINTER(adj_dictionary), STRING, STRING_POINTER]
adj_dict_print = _library.adj_dict_print
adj_dict_print.restype = None
adj_dict_print.argtypes = [POINTER(adj_dictionary)]
adj_dict_destroy = _library.adj_dict_destroy
adj_dict_destroy.restype = c_int
adj_dict_destroy.argtypes = [POINTER(adj_dictionary)]
class adj_eps_options(Structure):
    pass
adj_eps_options._fields_ = [
    ('input', adj_vector),
    ('output', adj_vector),
    ('method', c_char_p),
    ('type', c_int),
    ('which', c_int),
    ('monitor', c_int),
    ('neigenpairs', c_int),
]
class adj_eps(Structure):
    pass
adj_eps._fields_ = [
    ('eps_handle', c_void_p),
    ('eps_data', c_void_p),
]
adj_compute_eps = _library.adj_compute_eps
adj_compute_eps.restype = c_int
adj_compute_eps.argtypes = [POINTER(adj_adjointer), adj_matrix, adj_eps_options, POINTER(adj_eps), POINTER(c_int)]
adj_get_eps = _library.adj_get_eps
adj_get_eps.restype = c_int
adj_get_eps.argtypes = [POINTER(adj_eps), c_int, POINTER(c_double), POINTER(c_double), POINTER(adj_vector), POINTER(adj_vector)]
adj_destroy_eps = _library.adj_destroy_eps
adj_destroy_eps.restype = c_int
adj_destroy_eps.argtypes = [POINTER(adj_eps)]
adj_error_msg = (c_char * 1024).in_dll(_library, 'adj_error_msg')
adj_chkierr_private = _library.adj_chkierr_private
adj_chkierr_private.restype = None
adj_chkierr_private.argtypes = [c_int, STRING, c_int]
adj_chkierr_auto_private = _library.adj_chkierr_auto_private
adj_chkierr_auto_private.restype = c_int
adj_chkierr_auto_private.argtypes = [c_int, STRING, c_int]
adj_set_error_checking = _library.adj_set_error_checking
adj_set_error_checking.restype = c_int
adj_set_error_checking.argtypes = [c_int]
adj_evaluate_functional_derivative = _library.adj_evaluate_functional_derivative
adj_evaluate_functional_derivative.restype = c_int
adj_evaluate_functional_derivative.argtypes = [POINTER(adj_adjointer), adj_variable, STRING, POINTER(adj_vector), POINTER(c_int)]
adj_evaluate_functional = _library.adj_evaluate_functional
adj_evaluate_functional.restype = c_int
adj_evaluate_functional.argtypes = [POINTER(adj_adjointer), c_int, STRING, POINTER(c_double)]
class adj_gst(Structure):
    pass
adj_gst._fields_ = [
    ('eps_handle', c_void_p),
    ('gst_data', c_void_p),
]
adj_compute_gst = _library.adj_compute_gst
adj_compute_gst.restype = c_int
adj_compute_gst.argtypes = [POINTER(adj_adjointer), adj_variable, POINTER(adj_matrix), adj_variable, POINTER(adj_matrix), c_int, POINTER(adj_gst), POINTER(c_int), c_int]
adj_get_gst = _library.adj_get_gst
adj_get_gst.restype = c_int
adj_get_gst.argtypes = [POINTER(adj_gst), c_int, POINTER(c_double), POINTER(adj_vector), POINTER(adj_vector), POINTER(c_double)]
adj_destroy_gst = _library.adj_destroy_gst
adj_destroy_gst.restype = c_int
adj_destroy_gst.argtypes = [POINTER(adj_gst)]
adj_test_assert = _library.adj_test_assert
adj_test_assert.restype = None
adj_test_assert.argtypes = [c_int, STRING]
adj_sizeof_adjointer = _library.adj_sizeof_adjointer
adj_sizeof_adjointer.restype = c_int
adj_sizeof_adjointer.argtypes = []
adj_find_variable_equation_nb = _library.adj_find_variable_equation_nb
adj_find_variable_equation_nb.restype = c_int
adj_find_variable_equation_nb.argtypes = [POINTER(adj_adjointer), POINTER(adj_variable), POINTER(c_int)]
class UT_hash_bucket(Structure):
    pass
ptrdiff_t = c_long
uint32_t = c_uint32
UT_hash_table._fields_ = [
    ('buckets', POINTER(UT_hash_bucket)),
    ('num_buckets', c_uint),
    ('log2_num_buckets', c_uint),
    ('num_items', c_uint),
    ('tail', POINTER(adj_hash_handle)),
    ('hho', ptrdiff_t),
    ('ideal_chain_maxlen', c_uint),
    ('nonideal_items', c_uint),
    ('ineff_expands', c_uint),
    ('noexpand', c_uint),
    ('signature', uint32_t),
]
UT_hash_bucket._fields_ = [
    ('hh_head', POINTER(adj_hash_handle)),
    ('count', c_uint),
    ('expand_mult', c_uint),
]
__all__ = ['adj_variable_known', 'CACTION_YOUTURN',
           'adj_adjointer_check_checkpoints',
           'adj_get_forward_equation',
           'adj_equation_set_rhs_second_derivative_action_callback',
           'adj_block', 'adj_parameter_source_callback',
           'adj_storage_memory_incref', 'adj_destroy_gst',
           'adj_advance_to_adjoint_run_revolve', 'adj_get_finished',
           'adj_timestep_get_times', 'adj_get_adjoint_solution',
           'adj_register_parameter_source_callback',
           'adj_func_deriv_callback_list', 'adj_op_callback_list',
           'adj_get_adjoint_equation', 'CACTION',
           'adj_adjointer_to_html', 'adj_set_finished',
           'adj_get_variable_value', 'adj_term',
           'adj_block_set_coefficient',
           'adj_nonlinear_block_set_test_derivative',
           'adj_find_variable_equation_nb', 'adj_dict_destroy',
           'adj_create_equation', 'CACTION_RESTORE',
           'adj_set_revolve_options', 'adj_timestep_set_times',
           'adj_register_equation', 'adj_record_variable',
           'adj_nonlinear_block_set_test_hermitian',
           'adj_set_checkpoint_strategy', 'adj_adjointer',
           'adj_create_term', 'adj_test_assert', 'UT_hash_bucket',
           'adj_storage_memory_copy', 'size_t', 'adj_reset_revolve',
           'adj_get_tlm_equation', 'adj_add_term_to_equation',
           'adj_variable', 'adj_chkierr_auto_private',
           'CACTION_ADVANCE', 'adj_eps', 'adj_storage_disk',
           'adj_hash_handle', 'adj_equation',
           'adj_register_functional_callback',
           'adj_forget_tlm_equation', 'adj_storage_data',
           'adj_create_variable',
           'adj_equation_set_rhs_derivative_action_callback',
           'adj_set_revolve_debug_options', 'adj_functional_data',
           'adj_data_callbacks', 'adj_dictionary_entry',
           'adj_iteration_count', 'adj_add_terms',
           'adj_forget_adjoint_equation',
           'adj_register_data_callback', 'adj_set_error_checking',
           'adj_func_second_deriv_callback',
           'adj_storage_set_compare', 'adj_op_callback',
           'adj_block_set_test_hermitian',
           'adj_storage_set_overwrite', 'adj_get_forward_solution',
           'adj_nonlinear_block_set_coefficient', 'CACTION_ERROR',
           'ptrdiff_t', 'adj_nonlinear_block_derivative',
           'adj_storage_set_checkpoint', 'adj_create_nonlinear_block',
           'adj_dict_find', 'adj_variable_set_auxiliary',
           'adj_destroy_nonlinear_block', 'adj_variable_hash',
           'adj_func_callback_list', 'adj_get_soa_equation',
           'adj_parameter_source_callback_list',
           'adj_nonlinear_block_second_derivative', 'adj_compute_eps',
           'adj_destroy_block', 'adj_dictionary', 'CACTION_TERMINATE',
           'adj_timestep_set_functional_dependencies',
           'adj_timestep_data', 'adj_get_gst', 'adj_dict_set',
           'CRevolve', 'adj_compute_gst',
           'adj_register_functional_second_derivative_callback',
           'adj_equation_set_rhs_callback', 'adj_nonlinear_block',
           'CACTION_TAKESHOT', 'adj_timestep_start_equation',
           'adj_variable_get_iteration', 'UT_hash_table',
           'adj_func_deriv_callback', 'adj_chkierr_private',
           'adj_equation_count', 'adj_destroy_equation',
           'adj_create_adjointer', 'adj_get_soa_solution',
           'adj_vector', 'adj_variable_equal', 'adj_revolve_data',
           'adj_block_set_hermitian', 'adj_gst',
           'adj_get_forward_variable', 'adj_variable_get_name',
           'adj_forget_forward_equation', 'CACTION_FIRSTRUN',
           'adj_error_msg', 'adj_timestep_count',
           'adj_variable_get_type',
           'adj_equation_set_rhs_dependencies',
           'adj_destroy_adjointer', 'adj_dict_init',
           'adj_create_block', 'adj_matrix',
           'adj_deactivate_adjointer', 'adj_forget_tlm_values',
           'adj_variable_data', 'adj_func_callback',
           'adj_destroy_term', 'adj_evaluate_functional',
           'adj_eps_options', 'adj_get_tlm_solution',
           'adj_adjointer_check_consistency', 'adj_dict_print',
           'adj_func_second_deriv_callback_list',
           'adj_evaluate_functional_derivative',
           'adj_forget_adjoint_values', 'adj_variable_get_timestep',
           'adj_equation_set_rhs_derivative_assembly_callback',
           'adj_register_operator_callback',
           'adj_get_checkpoint_strategy', 'adj_get_eps',
           'adj_destroy_eps', 'adj_variable_str',
           'adj_register_functional_derivative_callback', 'uint32_t',
           'adj_sizeof_adjointer', 'adj_timestep_end_equation']
