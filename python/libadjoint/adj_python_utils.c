#include <python2.7/Python.h>

static PyObject * incref(PyObject *self, PyObject *args)
{
  /* Expose the INCREF macro to Python so as to enable libadjoint data callbacks to protect their targets from deallocation when control passes to C. */ 
  PyObject * object;
  
  if (!PyArg_ParseTuple(args, "O", &object))
    return NULL;

  Py_INCREF(object);

  Py_RETURN_NONE;
};

static PyObject * decref(PyObject *self, PyObject *args)
{
  /* Expose the DECREF macro to Python so as to enable libadjoint data callbacks to protect their targets from deallocation when control passes to C. */ 
  PyObject * object;
  
  if (!PyArg_ParseTuple(args, "O", &object))
    return NULL;

  Py_DECREF(object);

  Py_RETURN_NONE;
};

static PyObject * c_ptr(PyObject *self, PyObject *args)
{
  /* Return a python integer containing the pointer to the object provided.*/
  PyObject * object;
  
  if (!PyArg_ParseTuple(args, "O", &object))
    return NULL;

  /* Note that we purposefully do not dereference object. */ 
  return PyInt_FromLong((long)object);
};

static PyObject * c_deref(PyObject *self, PyObject *args)
{
  /* Given a python integer containing the pointer to a python object, return the object*/
  PyObject * pointer;
  
  
  if (!PyArg_ParseTuple(args, "O", &pointer))
    return NULL;
  
  /* The value of the python integer is the pointer */
  PyObject* pyobj = (PyObject*) PyInt_AsLong(pointer);
  Py_XINCREF(pyobj);
  return pyobj;
};

static PyMethodDef AdjMethods[] = {
/*    {"incref",  incref, METH_VARARGS,
     "Increment the reference count of a python object."},
    {"decref",  decref, METH_VARARGS,
     "Decrement the reference count of a python object."}, */
    {"c_ptr",  c_ptr, METH_VARARGS,
     "Find the memory address of a python object."},
    {"c_deref",  c_deref, METH_VARARGS,
     "Cast a memory address to a python object."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initpython_utils(void)
{
    PyObject *m;

    m = Py_InitModule("python_utils", AdjMethods);
    if (m == NULL)
        return;

}
