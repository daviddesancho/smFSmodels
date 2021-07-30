#
"""
Run using :
    python setup.py build_ext --build-lib smfsmodels
"""
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extension = [
        Extension("schutte/schutte", ["schutte/schutte.pyx"]
            )
] 

setup(
        name='schutte',
        url='https://github.com/imitxelena003',
        author='Ion Mitxelena',
        author_email='ifrentzua.at.gmail.com',
        license='MIT',
        ext_modules = cythonize(["schutte/*.pyx"])
)
