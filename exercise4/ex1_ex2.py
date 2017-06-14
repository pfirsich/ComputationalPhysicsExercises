import sys
import math
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    INTEGRATION_TYPE = sys.argv[1]
else:
    INTEGRATION_TYPE = "euler"

if INTEGRATION_TYPE not in ["euler", "eulerca", "eulercb"]:
    print("Invalid integration type!")
    quit()

def a(x):
    return -x # k = 1

def integrate_euler(x, v, dt):
    xn = x + v * dt
    vn = v + a(x) * dt
    return xn, vn

def integrate_euler_cromer_a(x, v, dt):
    vn = v + dt * a(x)
    xn = x + dt * vn
    return xn, vn

def integrate_euler_cromer_b(x, v, dt):
    xn = x + v * dt
    vn = v + a(xn) * dt
    return xn, vn

maxT = 1000
for dt in [0.1, 0.01, 0.001]:
    N = math.floor(1000/dt)
    t = np.linspace(0, N*dt, num=N)
    x, v = np.zeros(N), np.ones(N)
    for i in range(1, N):
        if INTEGRATION_TYPE == "euler":
            x[i], v[i] = integrate_euler(x[i-1], v[i-1], dt)
        elif INTEGRATION_TYPE == "eulerca":
            x[i], v[i] = integrate_euler_cromer_a(x[i-1], v[i-1], dt)
        elif INTEGRATION_TYPE == "eulercb":
            x[i], v[i] = integrate_euler_cromer_b(x[i-1], v[i-1], dt)
        else:
            print("Invalid integration type!")
            quit()

    plt.plot(t, x, label="dt = {}".format(dt))

t = np.linspace(0, maxT, 100000)
plt.plot(t, np.sin(t), linestyle="--", label="exact", linewidth=2)

plt.ylim(-2, 2)
plt.legend()
plt.xlabel("t")
plt.ylabel("x(t)")
plt.savefig("ex1_ex2_{}.pdf".format(INTEGRATION_TYPE))
plt.show()