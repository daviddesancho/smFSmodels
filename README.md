# smFS 
This module is an implementation of methods to analyze and interpret single molecule force 
spectroscopy (smFS) data. It is largely based on the calculations reported by Cossio, Hummer
 and Szabo in [PNAS (2015)](http://dx.doi.org/10.1073/pnas.1519633112).
It includes a module for the estimation of the kinetics from the results.
The code is distributed under the MIT License. Use at your own risk.


Installation
------------
Simply download the code and leave it in your working directory. To generate the 
extension from the `pyx` file simply run 

     python setup.py build_ext --build-lib smfsmodels

This should result in a cossio.so file that you can then import.

Example
-------
Jupyter notebooks with working examples are included.
