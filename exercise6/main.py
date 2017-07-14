import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

sigma = 3
x0 = 20
q = 1

Delta = 0.1
L = 1001
tau = 0.001
m = 50000

rescale = True
#rescale = False
V_barrier = 2
#V_barrier = 0

def V(x):
    if x >= 50.0 and x <= 50.5:
        return V_barrier
    else:
        return 0

def psi_t0(x):
    ret = math.pow(2 * math.pi * sigma*sigma, -0.25)
    ret *= cmath.exp(1j * q * (x-x0))
    ret *= math.exp(-(x-x0)*(x-x0)/(4*sigma*sigma))
    return ret

def pd(psi):
    v = np.absolute(psi)
    return v*v

psi = np.zeros(L, dtype=complex)
for l in range(L):
    psi[l] = psi_t0(l*Delta)

expV = np.zeros((L,L), dtype=complex)
for l in range(L):
    expV[l][l] = cmath.exp(-1j * tau * V(l*Delta))

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

x = np.arange(L) * Delta
plt.xlabel("x")
plt.ylabel("P(x,t)")
plt.grid(True)

ax = plt.gca()
ylim = 0.14
plt.ylim(0, ylim)
if V_barrier > 0:
    ax.fill_between([50.0, 50.5], 0, ylim, facecolor='grey', alpha=0.5)

right_index = math.floor(50.5 / Delta + 0.5)
right_normalization = 0
snapshots = []
for i in range(m+1):
    if math.floor(i/m*100) != math.floor((i-1)/m*100):
        print(math.floor(i/m*100), "%")
    t = i * tau

    # snapshot
    if i in [0, m/10, m/5*4, m/50*45, m]:
        total_prob = np.sum(pd(psi)) * Delta
        right_prob = np.sum(pd(psi[right_index:])) * Delta
        if right_prob > right_normalization:
            right_normalization = right_prob
        snapshots.append((t, pd(psi)))

    psi = np.dot(U, psi)

for t, pdpsi in snapshots:
    if rescale:
        pdpsi[right_index:] /= right_normalization
    plt.plot(x, pdpsi, label="t={}".format(t))

plt.legend()
plt.savefig("tdse_V={}_rescale_right={}.pdf".format(V_barrier, rescale))
plt.show()