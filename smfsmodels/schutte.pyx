#!/usr/bin/env python
import numpy as np
from scipy.stats import norm

"""
Implementation of the Schutte potential.
Schutte et al, J. Chem. Phys. (2011)

"""

# Globals
cdef float gamma = 1. # friction coefficient
cdef float beta = 4.
cdef float Diff = 1./(beta*gamma)

def f(float x):
    """
    Functional form of the potential

    x : float
        Value of molecular coordinate x.

    """
    if x < 0.0:
        return (1. - x*x)*(1. - x*x)
    elif x > 8.0:
        return (1. - (x - 8.)*(x - 8.))*(1. - (x - 8.)*(x - 8.))
    else:
        return (4./5. + 1./5. * np.cos(x*np.pi))

def df(float x):
    """
    Derivative of the potential

    x : float
        Value of molecular coordinate x.

    """
    if x < 0.0:
        return (4.*(x*x - 1.)*x)
    elif x > 8.:
        return (4.*((x - 8.)*(x - 8.) - 1)*(x - 8.))
    else:
        return (-np.pi/5. * np.sin(x*np.pi))

def delta_x_eff(float x, float dt, float rg):
    """
    Displacement in x

    x : float
        Value of molecular coordinate x.

    rg : float
        Normally distributed random number.

    dt : float
        Time step.

    """
    return (-beta*Diff*dt*df(x) + np.sqrt(2*Diff*dt)*rg)
    #return -gamma*beta*dt*df(x) + np.sqrt(2.*gamma*dt)*rg 

def run_brownian(float x0=5., float dt=5e-4, int numsteps=100000, int fwrite=1):
    """
    Brownian dynamics runner (see Cossio, Hummer, Szabo, PNAS (2015))
    for a potential from Schutte et al. JCP 2011.
    

    Parameters
    ----------
    x0 : float
        Initial position on molecular coordinate.

    dt : float
        Timestep for BD integraiton.

    numsteps : int
        Length of the run.

    fwrite : int
        Frequency of writing output.

    """

    cdef int kk, k
    cdef float x

    rgaussx = norm.rvs(size=fwrite, loc=0.0, scale=1.0)
    
    x = x0
    k = 0
    kk = 0
    t = 0.
    xk = [x]
    time = [0.]
    while True:
        x += delta_x_eff(x, dt, rgaussx[kk])
        t += dt
        k +=1
        kk +=1
        if k%fwrite == 0:
            time.append(t)
            xk.append(x)
            rgaussx = norm.rvs(size=fwrite, loc=0.0, scale=1.0)
            kk = 0
            if k >= numsteps:
                break
    return time, xk

def fmod(float x):
    """
    Functional form of the potential

    x : float
        Value of molecular coordinate x.

    """
    if x < 0.0:
        return (1.-x*x)*(1.-x*x)
    elif x > 8.0:
        return (1.-(x-8.)*(x-8.))*(1.-(x-8.)*(x-8.))
    else:
        return 2.7/5. + 2.3/5. * np.cos(x*np.pi)

def dfmod(float x):
    """
    Derivative of the potential

    x : float
        Value of molecular coordinate x.

    """
    if x < 0.0:
        return 4*(x*x-1)*x
    elif x > 8.:
        return 4.*((x-8.)*(x-8.)-1)*(x-8.)
    else:
        return -2.3*np.pi/5. * np.sin(x*np.pi)
