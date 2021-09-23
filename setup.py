#
"""
Run using :
    python setup.py build_ext --build-lib smfsmodels
"""
import sys
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extension = [
        Extension("smfsmodels/cossio", ["smfsmodels/cossio.pyx"]
            ),
        Extension("smfsmodels/schutte", ["smfsmodels/schutte.pyx"]
            )
] 

setup(
        name='smfs',
        url='https://github.com/daviddesancho/smFSmodels',
        author='David De Sancho',
        author_email='daviddesancho.at.gmail.com',
        license='MIT',
#        ext_modules = cythonize(["smfsmodels/cossio.pyx", "smfsmodels/schutte.pyx"])
        ext_modules = cythonize(["smfsmodels/*.pyx"], \
                compiler_directives={'language_level' : sys.version_info[0]})
)
