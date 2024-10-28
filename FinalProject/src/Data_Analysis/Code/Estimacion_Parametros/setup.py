from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize("LB_AdvectionDiffusion.pyx", compiler_directives={'language_level': "3"}),
    include_dirs=[numpy.get_include()]
)

#python setup.py build_ext --inplace
