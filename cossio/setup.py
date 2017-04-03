#
"""
Run using :
    python setup.py build_ext --inplace
"""
from distutils.core import setup
from Cython.Build import cythonize

setup(
        ext_modules = cythonize("*.pyx"),
)
#setup(
#            ext_modules = cythonize("cossio.pyx",
#                    compiler_directives={'profile': True})
#            )

