#!/usr/bin/python

f = open('include/libadjoint/adj_constants.h', 'r')
constants = {}

for line in f:
  if "#define " in line.lower():
    try:
      m, k, v = line.split()
      constants[k] = v
    except ValueError:
      pass
f.close()

code = 'adj_constants = %s' % repr(constants)

fc = open('lib/clibadjoint_constants.py', 'w')
fc.write(code)
fc.close()
