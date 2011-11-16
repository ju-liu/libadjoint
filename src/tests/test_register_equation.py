#!/usr/bin/python
import libadjoint

def my_vec_duplicate(a, b):
  print "Vec duplicate"

adjointer = libadjoint.Adjointer()
adjointer.register_data_callback('ADJ_VEC_DUPLICATE_CB', my_vec_duplicate)

