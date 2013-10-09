#!/usr/bin/env python

import sys
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'
else:
    path = ''
print path

f = open(path + 'include/libadjoint/adj_constants.h', 'r')
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

fc = open(path + 'python/libadjoint/clibadjoint_constants.py', 'w')
fc.write(code)
fc.close()
