from distutils.core import setup, Extension

python_utils=Extension("libadjoint.python_utils",
                       sources=["libadjoint/adj_python_utils.c"]
                       )

setup (name = 'libadjoint',
       version = '2016.1.0',
       description = 'libadjoint python bindings',
       author = 'The libadjoint team',
       author_email = 'patrick.farrell@maths.ox.ac.uk',
       packages = ['libadjoint'],
       package_dir = {'libadjoint': 'libadjoint'},
       ext_modules = [python_utils])
