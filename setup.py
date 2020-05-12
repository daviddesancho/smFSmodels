#
"""
Run using :
    python setup.py build_ext --build-lib smfsmodels
"""
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extension = [
        Extension("smfsmodels/cossio", ["smfsmodels/cossio.pyx"]
            )
] 

setup(
        name='smfs',
        url='https://github.com/daviddesancho/smFSmodels',
        author='David De Sancho',
        author_email='daviddesancho.at.gmail.com',
        license='MIT',
        ext_modules = cythonize(["smfsmodels/*.pyx"])
)
