# Cossio
This module is an implementation of the calculations reported by Cossio, Hummer and
Szabo in [PNAS (2015)](http://dx.doi.org/10.1073/pnas.1519633112) and it is 
distributed under the MIT License. It includes a module for the estimation of the 
kinetics from the results. Use at your own risk.


Installation
------------
Simply download the code and leave it in your working directory. To generate the 
extension from the `pyx` file simply run 

    python setup.py build_ext --inplace

This should result in a cossio.so file that you can then import.

Example
-------
A Jupyter notebook with a working example reproducing some of the results of the paper is 
included.
