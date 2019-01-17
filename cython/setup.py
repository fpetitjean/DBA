from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize("DBA.pyx",compiler_directives={'boundscheck':False,'wraparound':False}),
    include_dirs=[numpy.get_include()]
)
