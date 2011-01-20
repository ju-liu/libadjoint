#!/usr/bin/env sage

# This is a sage script (http://www.sagemath.org)
# that, given a matrix that is a function of some variables,
# computes its derivative with respect to those variables
# contracted with some other vector.

def differentiate(m, x, c):
  d = m.new_matrix()

  for j in range(len(x)):
    d[:, j] = diff(m, x[j]) * c

  return d

if __name__ == "__main__":
  n = 3
  x = []
  m = matrix(SR, n, n)

  for i in range(n):
    x.append(var('x%d' % i))
    m[i,i] = x[i]

  c = vector([float(i+1) for i in range(n)])

  print "Matrix to differentiate: "
  print m
  
  print "Contraction vector: "
  print c
  d = differentiate(m, x, c)
  # d should be c on the diagonal
  print "dM/dx * c: "
  print d
