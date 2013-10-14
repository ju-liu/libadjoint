#!/usr/bin/env python

import sys
if len(sys.argv) > 1:
    source_path = sys.argv[1] + '/'
else:
    source_path = ''
print source_path

if len(sys.argv) > 2:
    target_path = sys.argv[2] + '/'
else:
    target_path = source_path
print target_path

f = open(source_path + 'include/libadjoint/adj_constants.h', 'r')
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

fc = open(target_path + '/clibadjoint_constants.py', 'w')
fc.write(code)
fc.close()
