#
"""
Run using :
    python setup.py build_ext --build-lib smfs
"""
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

extension = [
        Extension("smfs/cossio", ["smfs/cossio.pyx"]
            )
] 

#,
#        ("smfs/cossio_ramp",
#            ["smfs/cossio_ramp.pyx"]
#    )
#]

setup(
        name='smfs',
        url='https://github.com/daviddesancho/smFS',
        author='David De Sancho',
        author_email='daviddesancho.at.gmail.com',
        license='MIT',
#        ext_modules = cythonize(extension)
        ext_modules = cythonize(["smfs/*.pyx"])
)
