#!/usr/bin/env python
import numpy as np
from scipy.stats import norm

# Globals
cdef float xddagger = 4 # in nm
cdef float beta = 1.0

def Gx(float x, float barrier=5.):
    """
    Functional form of the molecular free energy surface
    
    """
    return barrier*f(x)

def V(float q, float x, float kl):
    """
    Harmonic contribution from the instrument

    """
    return 0.5*kl*(q - x)**2

def Gqx(float q, float x, float barrier=5.0, float kl=0.):
    """
    Functional form of the 2D free energy surface

    """
    return Gx(x, barrier) + V(q, x, kl)

def f(float x):
    """
    Functional form of the bistable potential

    """
    cdef float xxd

    xxd = x/xddagger
    if abs(xxd) < 0.5:
        return -2*xxd**2  
    else:
        return 2*(np.abs(xxd)-1)**2 - 1

def df(float x):
    """
    Derivative of the bistable potential

    """
    cdef float xxd

    xxd = x/xddagger
    if abs(xxd) < 0.5:
        return -4.*x/xddagger**2  
    else:
        return 4.*x/(xddagger**2)*(1. - 1./np.abs(xxd))
        #return 4*x*(np.abs(x)-1)/np.abs(x) 
#    if abs(x) < 0.5:
#        return -4*x
#    else:
#        return 4*x*(np.abs(x)-1)/np.abs(x) 
 
def dVdx(float q, float x, float kl):
    """
    Derivative of instrumental contribution with respect to x

    """
    return -kl*(q - x)

def dGqxdq(float q, float x, float kl=0.):
    """
    Derivative of the 2D free energy surface with respect to q

    """
    return kl*(q - x)

def dGqxdx(float q, float x, float barrier=5., float kl=0.):
    """
    Derivative of the 2D free energy surface with respect to x

    """
    return barrier*df(x) + dVdx(q, x, kl)

def delta_q_eff(float q, float x, float dt, float kl, float rg,\
        float Dq, float sqrt2Dqdt):
    """
    Displacement in q

    """
    return -dt*Dq*beta*dGqxdq(q,x,kl) + sqrt2Dqdt*rg

def delta_x_eff(float q, float x, float dt, float kl, float rg, \
        float Dx, float sqrt2Dxdt):
    """
    Displacement in x

    """
    return -dt*Dx*beta*dGqxdx(q,x,kl) + sqrt2Dxdt*rg 

def run_brownian(float x0=5., float q0=0., float Dx=0., float Dq=0.,\
        float kl=0., float dt=1.e-3, int numsteps=100000):
    """
    Brownian dynamics runner for anisotropic diffusion model.
    Cossio, Hummer, Szabo, PNAS (2015)

    Parameters
    ----------
    x0 : float
        Initial position on molecular coordinate.

    q0 : float
        Initial position on measuring coordinate.

    Dx : float
        Diffusion coefficient for molecule.

    Dq : float
        Diffusion coefficient for instrument.

    kl : float
        Spring constant of the instrument.

    dt : float
        Timestep for BD integraiton.

    numsteps :int
        Length of the run.
    """

    cdef int k
    cdef float x, q
    cdef float sqrt2Dxdt, sqrt2Dqdt

    sqrt2Dxdt = np.sqrt(2.*Dx*dt)
    sqrt2Dqdt = np.sqrt(2.*Dq*dt)

    rgaussx = norm.rvs(size=numsteps, loc=0.0, scale=1.0)
    rgaussq = norm.rvs(size=numsteps, loc=0.0, scale=1.0)
    
    x, q = [x0, q0]
    k = 0
    xk = [x]
    qk = [q]
    while True:
        x += delta_x_eff(q, x, dt, kl, rgaussx[k], Dx, sqrt2Dxdt)
        q += delta_q_eff(q, x, dt, kl, rgaussq[k], Dq, sqrt2Dqdt)
        k +=1
        if k%100 == 0:
            xk.append(x)
            qk.append(q)
        if k == numsteps:
            break
    return xk, qk
