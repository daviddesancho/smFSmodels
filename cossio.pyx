#!/usr/bin/env python
import numpy as np
from scipy.stats import norm

# Globals
cdef float xddagger = 4 # in nm
cdef float beta = 1.0

def Gx(float x, float barrier=5.):
    """
    Functional form of the molecular free energy surface

    x : float
        Value of molecular coordinate x.

    barrier : float
        Value of free energy barrier on molecular coordinate.
    
    """
    return barrier*f(x)

def V(float q, float x, float kl=0.):
    """
    Harmonic contribution from the instrument

    q : float
        Value of measured coordinate q.

    x : float
        Value of molecular coordinate x.

    kl : float
        Spring constant of harmonic potential.

    """
    return 0.5*kl*(q - x)**2

def Gqx(float q, float x, float barrier=5.0, float kl=0.):
    """
    Functional form of the 2D free energy surface

    q : float
        Value of measured coordinate q.

    x : float
        Value of molecular coordinate x.

    barrier : float
        Value of free energy barrier on molecular coordinate.
    
    kl : float
        Spring constant of harmonic potential.

    """
    return Gx(x, barrier) + V(q, x, kl)

def f(float x):
    """
    Functional form of the bistable potential

    x : float
        Value of molecular coordinate x.

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

    x : float
        Value of molecular coordinate x.

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
 
def dVdx(float q, float x, float kl=0.):
    """
    Derivative of instrumental contribution with respect to x

    q : float
        Value of measured coordinate q.

    x : float
        Value of molecular coordinate x.

    kl : float
        Spring constant of harmonic potential.

    """
    return -kl*(q - x)

def dGqxdq(float q, float x, float kl=0.):
    """
    Derivative of the 2D free energy surface with respect to q

    q : float
        Value of measured coordinate q.

    x : float
        Value of molecular coordinate x.

    kl : float
        Spring constant of harmonic potential.

    """
    return kl*(q - x)

def dGqxdx(float q, float x, float barrier=5., float kl=0.):
    """
    Derivative of the 2D free energy surface with respect to x

    q : float
        Value of measured coordinate q.

    x : float
        Value of molecular coordinate x.

    barrier : float
        Value of free energy barrier on molecular coordinate.
    
    kl : float
        Spring constant of harmonic potential.

    """
    return barrier*df(x) + dVdx(q, x, kl)

def delta_q_eff(float q, float x, float bDqdt, float sqrt2Dqdt, \
        float rg, float kl):
    """
    Displacement in q

    q : float
        Value of measured coordinate q.

    x : float
        Value of molecular coordinate x.

    kl : float
        Spring constant of harmonic potential.

    rg : float
        Normally distributed random number.

    bDxDt : float
        beta*Dx*dt

    sqrt2Dqdt : float
        sqrt(2*Dq*dt)

    """
    return -bDqdt*dGqxdq(q,x, kl) + sqrt2Dqdt*rg

def delta_x_eff(float q, float x, float barrier, \
        float bDxdt, float sqrt2Dxdt, float rg, float kl):
    """
    Displacement in x

    q : float
        Value of measured coordinate q.

    x : float
        Value of molecular coordinate x.

    barrier : float
        Value of free energy barrier on molecular coordinate.
    
    kl : float
        Spring constant of harmonic potential.

    rg : float
        Normally distributed random number.

    bDxdt : float
        beta*Dx*dt

    sqrt2Dxdt : float
        sqrt(2*Dx*dt)

    """
    return -bDxdt*dGqxdx(q,x,barrier, kl) + sqrt2Dxdt*rg 

def run_brownian(float x0=5., float q0=0., float dt=5e-4, float barrier=5., \
        float kl=0., float Dx=0., float Dq=0., int numsteps=100000, int fwrite=1):
    """
    Brownian dynamics runner for anisotropic diffusion model.
    Cossio, Hummer, Szabo, PNAS (2015)

    Parameters
    ----------
    x0 : float
        Initial position on molecular coordinate.

    q0 : float
        Initial position on measuring coordinate.

    barrier : float
        Value of free energy barrier on molecular coordinate.
    
    Dx : float
        Diffusion coefficient for molecule.

    Dq : float
        Diffusion coefficient for instrument.

    kl : float
        Spring constant of the instrument.

    dt : float
        Timestep for BD integraiton.

    numsteps : int
        Length of the run.

    fwrite : int
        Frequency of writing output.

    """

    cdef int k, kk
    cdef float x, q
    cdef float sqrt2Dxdt, sqrt2Dqdt
    cdef float bDqdt, bDxdt    
    
    dt = 5.e-4
    bDqdt = beta*Dq*dt
    bDxdt = beta*Dx*dt
    sqrt2Dxdt = np.sqrt(2.*Dx*dt)
    sqrt2Dqdt = np.sqrt(2.*Dq*dt)

    rgaussx = norm.rvs(size=fwrite, loc=0.0, scale=1.0)
    rgaussq = norm.rvs(size=fwrite, loc=0.0, scale=1.0)
    
    x, q = [x0, q0]
    k = 0
    kk = 0
    t = 0.
    xk = [x]
    qk = [q]
    time = [0.]
    while True:
        x += delta_x_eff(q, x, barrier,  bDxdt, sqrt2Dxdt, rgaussx[kk], kl)
        q += delta_q_eff(q, x,  bDqdt, sqrt2Dqdt, rgaussq[kk], kl)
        t += dt
        k +=1
        kk +=1
        if k%fwrite == 0:
            time.append(t)
            xk.append(x)
            qk.append(q)
            rgaussx = norm.rvs(size=fwrite, loc=0.0, scale=1.0)
            rgaussq = norm.rvs(size=fwrite, loc=0.0, scale=1.0)
            kk = 0
            if k >= numsteps:
                break
    return time, xk, qk
