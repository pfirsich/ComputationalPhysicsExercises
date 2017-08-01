import sys
import math
import cmath

import numpy as np
import matplotlib.pyplot as plt

from conf import *

N = 200
x = np.linspace(xMin, xMax, num=N)
omega = 1

def hermite(n, x):
    if n == 0:
        return 1
    elif n == 1:
        return 2*x
    else:
        return 2*x*hermite(n-1, x) - 2*(n-1)*hermite(n-2, x)

def psi_0(x, x0, sigma):
    arg = (x - x0) / sigma
    return math.pow(math.pi*sigma*sigma, -0.25) * np.exp(-arg*arg/2)

def basestate(n, x):
    ret = math.pow(omega/math.pi, 1/4) / math.sqrt(math.pow(2, n) * math.factorial(n))
    return ret * np.exp(-1/2*omega*x*x) * hermite(n, math.sqrt(omega)*x)

def dot(a, b):
    return np.sum(np.conjugate(a)*b) * xRange/N

D = 30
i = np.arange(D)
v = np.zeros(D)
for j in range(D):
    v[j] = dot(psi_0(x, 0, 2), basestate(j, x))

plt.plot(i, v)
plt.show()
