#!/usr/bin/env python
import sys
import os
import os.path

sys.path.append(os.path.join(os.path.dirname(sys.argv[0]), os.pardir))

from ctypeslib.h2xml import main

if __name__ == "__main__":
    sys.exit(main())
