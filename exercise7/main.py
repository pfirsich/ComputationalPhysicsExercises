import math
import numpy as np
import matplotlib.pyplot as plt

invT1 = 0
invT2 = 0
M0 = [0, 1, 0]
phi = 0 # in units of pi

######

f0 = 4
f1 = 1/4
B0 = 2*math.pi*f0
h = 2*math.pi*f1
gamma = 1
omega0 = B0

t_min = 0
t_max = 4
t_range = t_max - t_min
steps = 1000
tau = t_range / steps

def _B(t):
    arg = omega0*t + math.pi*phi
    return np.asarray([h*math.cos(arg), -h*math.sin(arg), B0])

def _expC():
    expT1 = math.exp(-tau/2 * invT1)
    expT2 = math.exp(-tau/2 * invT2)
    return np.diag([expT2, expT2, expT1])

def _expB(t, B):
    omegaSqr = B[0]*B[0] + B[1]*B[1] + B[2]*B[2]
    omega = math.sqrt(omegaSqr)
    BB = np.outer(B, B)
    c = math.cos(omega*t)
    s = math.sin(omega*t)

    ret = BB * (1 - c)
    ret += np.eye(3) * omega * omega * c
    v = omega * s * B
    # essentially: M_ij = eps_ijk * v_k (eps being the levi-civita symbol)
    ret += np.asarray([[0,v[2],-v[1]], [-v[2],0,v[0]], [v[1],-v[0],0]])

    ret /= omegaSqr

    return ret

def expB(t):
    C = _expC()
    return C @ _expB(tau, _B(t + tau/2)) @ C

t = [t_min]
M = [np.asarray(M0)]
for i in range(steps):
    M.append(np.dot(expB(t[-1]), M[-1]))
    t.append(t[-1] + tau)
M = np.asarray(M)

plt.plot(t, M[:,0], label="$M^x$")
plt.plot(t, M[:,1], label="$M^y$")
plt.plot(t, M[:,2], label="$M^z$")
plt.xlabel("t")
plt.ylabel("Magnetization")
plt.xlim(t_min, t_max)
plt.ylim(-1,1)
plt.legend()
plt.grid()
plt.savefig("invT1={},invT2={},M0={}{}{},phi={}.pdf".format(invT1, invT2, M0[0], M0[1], M0[2], phi))
plt.show()
