#include <python2.7/Python.h>

static PyObject * incref(PyObject *self, PyObject *args)
{
  /* Expose the INCREF macro to Python so as to enable libadjoint data callbacks to protect their targets from deallocation when control passes to C. */ 
  PyObject * object;
  
  if (!PyArg_ParseTuple(args, "O", &object))
    return NULL;

  Py_INCREF(object);

};

static PyObject * decref(PyObject *self, PyObject *args)
{
  /* Expose the DECREF macro to Python so as to enable libadjoint data callbacks to protect their targets from deallocation when control passes to C. */ 
  PyObject * object;
  
  if (!PyArg_ParseTuple(args, "O", &object))
    return NULL;

  Py_DECREF(object);

};

static PyMethodDef SpamMethods[] = {
    {"incref",  spam_system, METH_VARARGS,
     "Increment the reference count of a python object."},
    {"decref",  spam_system, METH_VARARGS,
     "Decrement the reference count of a python object."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initadj_python_utils(void)
{
    PyObject *m;

    m = Py_InitModule("adj_python_utils", AdjMethods);
    if (m == NULL)
        return;

}
