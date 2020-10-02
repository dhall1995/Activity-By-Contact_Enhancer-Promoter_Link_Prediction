from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize("link_utils.pyx"))
setup(ext_modules=cythonize("dtrack_utils.pyx"))
