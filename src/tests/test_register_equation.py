#!/usr/bin/python
import libadjoint

def my_vec_duplicate(a, b):
  print "Vec duplicate"

adjointer = libadjoint.Adjointer()
adjointer.register_data_callback('vec_duplicate', my_vec_duplicate)

