import sys
import h5py
import numpy as np

pathtocossio = "/Users/daviddesancho/Research/code/smFSmodels"
sys.path.append(pathtocossio)
from smfsmodels import cossio


def cossio_runner(inp):
    np.random.seed()
    b = inp[0]
    kl = inp[1]
    sc = inp[2]
    numsteps = inp[3]

    dt = 5e-4
    Dx = 1
    Dq = sc*Dx
    x, q = [5., 5.]

    tt, xk, qk = cossio.run_brownian(x0=x, barrier=b, kl=kl, \
                            Dx=Dx, Dq=Dq, numsteps=numsteps, \
                                     fwrite=int(1./dt))
    data = np.column_stack((tt,xk,qk))

    h5file = "data/cossio_barrier%g_kl%g_Dx%g_Dq%g.h5"%(b, kl, Dx, Dq)
    with h5py.File(h5file, "w") as hf:
        hf.create_dataset("data", data=data)

    return h5file
