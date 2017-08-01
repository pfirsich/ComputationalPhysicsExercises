import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from conf import *
x = np.arange(L) * Delta + xMin

def pd(psi):
    v = np.absolute(psi)
    return v*v

def expect(psi, f):
    assert f.shape == psi.shape
    return np.dot(f, pd(psi)) * Delta

def mean_pos(psi):
    return expect(psi, x)

def var_pos(psi):
    mp = mean_pos(psi)
    return expect(psi, x*x) - mp*mp

data = np.load(sys.argv[1])
name = sys.argv[1][:-4]

t = np.arange(m+1) * tau
mean = np.zeros(m+1)
var = np.zeros(m+1)
for i in range(m+1):
    psi = data[i]
    # snapshots at t=0,2,4,6,8,10
    # => i = 0, 8000, 16000, ..., 40000
    if i % 8000 == 0:
         plt.plot(x, pd(psi), label="t={}".format(i * tau))

    mean[i] = mean_pos(psi)
    var[i] = var_pos(psi)

plt.xlabel("x")
plt.ylabel("P(x,t)")
plt.grid(True)
plt.legend()
plt.savefig(name + "-probdist.pdf")
plt.show()
plt.close()

plt.xlabel("t")
plt.ylabel("Average")
plt.grid(True)
plt.plot(t, mean, label="$<x(t)>$")
plt.plot(t, var, label="$<x(t)^2> - <x(t)>^2$")
plt.legend()
plt.savefig(name + "-pos.pdf")
plt.show()