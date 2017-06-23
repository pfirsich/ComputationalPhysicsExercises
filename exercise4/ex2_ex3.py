import sys
import math
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    INTEGRATION_TYPE = sys.argv[1]
else:
    INTEGRATION_TYPE = "verlet"

if INTEGRATION_TYPE not in ["verlet", "eulerca", "eulercb"]:
    print("Invalid integration type!")
    quit()

def a(x):
    return -x # k = 1

def integrate_euler_cromer_a(x, v, dt):
    vn = v + dt * a(x)
    xn = x + dt * vn
    return xn, vn

def integrate_euler_cromer_b(x, v, dt):
    xn = x + v * dt
    vn = v + a(xn) * dt
    return xn, vn

dt = 0.01
N = 1000
maxT = N * dt
t = np.linspace(0, maxT, num=N)
x, v = np.zeros(N), np.ones(N)

if INTEGRATION_TYPE == "verlet":
    # one step euler (-cromer A)
    v[1] = v[0] + a(x[0]) * dt
    x[1] = x[0] + v[1] * dt
    lastAccell = a(x[1])
    startIndex = 2
else:
    startIndex = 1

for i in range(startIndex, N):
    if INTEGRATION_TYPE == "verlet":
        # velocity verlet
        x[i] = x[i-1] + v[i-1]*dt + 0.5 * lastAccell * dt*dt
        accell = a(x[i])
        v[i] = v[i-1] + 0.5*(lastAccell + accell)*dt
        lastAccell = accell
    elif INTEGRATION_TYPE == "eulerca":
        x[i], v[i] = integrate_euler_cromer_a(x[i-1], v[i-1], dt)
    elif INTEGRATION_TYPE == "eulercb":
        x[i], v[i] = integrate_euler_cromer_b(x[i-1], v[i-1], dt)
    else:
        print("Invalid integration type!")
        quit()

plt.plot(t, x, label="x(t)".format(dt))
plt.plot(t, v, label="v(t)".format(dt))

plt.plot(t, np.sin(t), linestyle="--", label="exact")

Ee = np.full(N, 0.5)
plt.plot(t, Ee, linestyle="--", label="E exact")
plt.plot(t, np.power(v, 2)/2 + np.power(x, 2)/2)

plt.ylim(-1, 1)
plt.legend()
plt.xlabel("t")
plt.ylabel("x(t), v(t), E(t)")
plt.savefig("ex2_ex3_{}.pdf".format(INTEGRATION_TYPE))
plt.show()