#!/usr/bin/env python

import sys
import os.path

source_path = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), os.pardir)
print source_path

if len(sys.argv) > 1:
    target_path = sys.argv[1] + '/'
else:
    target_path = source_path
print target_path

f = open(os.path.join(source_path, 'include', 'libadjoint', 'adj_constants.h'), 'r')
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

fc = open(os.path.join(target_path, 'clibadjoint_constants.py'), 'w')
fc.write(code)
fc.close()
