#!/usr/bin/env python

from libadjoint.libadjoint_numpy import *

A=libadjoint.Adjointer()

v=Vector(numpy.random.rand(10))

var=libadjoint.Variable('foo', 0)

b=libadjoint.Block("Identity")

def id_assemble(variables, dependencies, hermitian, coefficient, context):
    return (Matrix(coefficient*numpy.eye(10)), Vector(coefficient*numpy.zeros(10)))

b.assemble=id_assemble

def rhs_cb(adjointer, variable, dependencies, values, context):
    return v

e=libadjoint.Equation(var, [b], [var], rhs_cb=rhs_cb)

A.register_equation(e)

A.record_variable(var, libadjoint.MemoryStorage(v))

var1=libadjoint.Variable('foo', 1)

def rhs_cb1(adjointer, variable, dependencies, values, context):
    return dependencies[0]
 
e1=libadjoint.Equation(var1, [b], [var1], rhs_cb=rhs_cb1, rhs_deps=[var])

A.register_equation(e1)

(var, soln0) = A.get_forward_solution(0)
libadjoint.adj_test_assert(all(soln0.vec[:] == v.vec[:]), "First solution should be v")

(var, soln1) = A.get_forward_solution(1)
libadjoint.adj_test_assert(all(soln0.vec[:] == soln1.vec[:]), "Second solution should also be v")
