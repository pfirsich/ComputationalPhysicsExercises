import sys
import math
import cmath

import numpy as np
import matplotlib.pyplot as plt

from conf import *

omega = 1

# r = 1
# phi = float(sys.argv[1])*math.pi
# if len(sys.argv) > 2:
#     r = float(sys.argv[2])
# alpha = cmath.rect(r, phi)
# print("alpha = {r} * exp({phi} * pi) = {comp}".format(r=r, phi=phi/math.pi, comp=alpha))

alpha = complex(sys.argv[1])
print("alpha =", alpha)

def hermite(n, x):
    if n == 0:
        return 1
    elif n == 1:
        return 2*x
    else:
        return 2*x*hermite(n-1, x) - 2*(n-1)*hermite(n-2, x)

def gaussian(x, mu, sig):
    e = (x - mu) / sig
    return 1/sig/math.sqrt(2*math.pi) * np.exp(-1/2 * e*e)

def cstate(x):
    ret = 0
    alphapow = 1
    for n in range(0, 20, 2):
        ret += alphapow / math.factorial(n) / math.sqrt(math.pow(2, n)) * hermite(n, math.sqrt(omega)*x)
        alphapow *= alpha
    return ret * math.exp(-abs(alpha)*abs(alpha) / 2) * math.exp(-1/2*omega*x*x) * math.pow(omega/math.pi, 1/4)

N = 200
x = np.linspace(xMin, xMax, num=N)
val = np.zeros(N, dtype=complex)

for i in range(N):
    val[i] = cstate(xMin + xRange / (N-1) * i)

absval = np.absolute(val)
curve = absval*absval

# calculate mean and stddev
mean = np.sum(x * curve) * xRange/N
var = np.sum(x*x*curve) * xRange/N - mean*mean
print("Mean: {}, Var: {} (stddev: {})".format(mean, var, math.sqrt(var)))
print("Expect: mean: {}".format(alpha.real * math.sqrt(2/omega)))

plt.plot(x, curve)
plt.plot(x, gaussian(x, mean, math.sqrt(var)), label="test")
plt.legend()
plt.show()