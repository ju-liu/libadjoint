#!/usr/bin/env python

from libadjoint.libadjoint_numpy import *

A=libadjoint.Adjointer()

v=Vector(numpy.random.rand(10))

var=libadjoint.Variable('foo', 0)

b=libadjoint.Block("Identity")

def id_assemble(variables, dependencies, hermitian, coefficient, context):
    return (Matrix(coefficient*numpy.eye(10)), Vector(coefficient*numpy.zeros(10)))

b.assemble=id_assemble

e=libadjoint.Equation(var, [b], [var], rhs=RHS(v))

A.register_equation(e)

A.record_variable(var, libadjoint.MemoryStorage(v))

(var, soln0) = A.get_forward_solution(0)
libadjoint.adj_test_assert(all(soln0.vec[:] == v.vec[:]), "First solution should be v")
