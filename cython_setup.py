# This script can be run using the following command:
#
#    python cython_setup.py build_ext --inplace
#

from setuptools import setup, Extension
from Cython.Build import cythonize
from numpy import get_include as get_numpy_include
from os import path

exts = [Extension("dtrack_utils", [path.join("utils/cython", "dtrack_utils.pyx")],
                  include_dirs=[get_numpy_include()]),
	Extension("link_utils", [path.join("utils/cython","link_utils.pyx")],
		  include_dirs=[get_numpy_include()])
]

setup(
    name="Activity by Contact link prediction",
    version="0.0.1",
    description="A tool to calculate ",
    packages=['dtrack_utils','link_utils'],
    ext_modules = cythonize(exts)
)

