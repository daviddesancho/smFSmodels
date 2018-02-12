#
"""
Run using :
    python setup.py build_ext --build-lib smfs

"""
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize

extension = [Extension
        ("smfs/cossio",
            ["smfs/cossio.pyx"]
    ),
]

setup(
        name='smfs',
        url='https://github.com/daviddesancho/smFS',
        author='David De Sancho',
        author_email='daviddesancho.at.gmail.com',
        license='MIT',
        packages=find_packages(),
        ext_modules = cythonize(extension),
)
