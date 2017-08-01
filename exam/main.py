import math
import cmath
import sys
import numpy as np
import matplotlib.pyplot as plt

from conf import *

def sim(Omega, sigma, x0):
    def V(x):
        return Omega*Omega / 2 * x*x

    def psi_t0(x):
        ret = math.pow(math.pi * sigma*sigma, -0.25)
        ret *= math.exp(-(x-x0)*(x-x0)/(2*sigma*sigma))
        return ret

    psi = np.zeros(L, dtype=complex)
    for l in range(L):
        psi[l] = psi_t0(xMin + l*Delta)

    expV = np.zeros((L,L), dtype=complex)
    for l in range(L):
        expV[l][l] = cmath.exp(-1j * tau * (1/(Delta*Delta) + V(xMin + l*Delta)))

    c = math.cos(tau / (4*Delta*Delta))
    s = math.sin(tau / (4*Delta*Delta))

    expK1 = np.eye(L, dtype=complex) * c
    expK1[L-1][L-1] = 1
    expK2 = np.eye(L, dtype=complex) * c
    expK2[1][1] = 1

    for l in range(L):
        if l+1 < L:
            if l % 2 == 0:
                expK1[l+1][l] = expK1[l][l+1] = 1j * s
            else:
                expK2[l+1][l] = expK2[l][l+1] = 1j * s

    U = expK1 @ expK2 @ expV @ expK2 @ expK1

    out = np.zeros((m+1, L), dtype=complex)

    for i in range(m+1):
        if math.floor(i/m*100) != math.floor((i-1)/m*100):
            print(math.floor(i/m*100), "%")
        t = i * tau

        psi = np.dot(U, psi)
        out[i] = psi

    return out

# on sheet: (1,1,0), (1,1,1), (1,2,0), (2,1,1), (2,2,2)
Omega, sigma, x0 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
set_name = "Omega={}, sigma={}, x0={}".format(Omega, sigma, x0)
print(set_name)
data = sim(Omega, sigma, x0)
np.save(set_name, data)