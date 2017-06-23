import sys
import math
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    BOUNDARY_TYPE = int(sys.argv[1])
else:
    BOUNDARY_TYPE = 2
print("Boundary type:", BOUNDARY_TYPE)

def a(x, i):
    if i == 0:
        return -(x[0] - x[1])
    elif i == x.size - 1:
        return -(x[i] - x[i-1])
    else:
        return -(2*x[i] - x[i-1] - x[i+1])

x_results = []

for n_i, N in enumerate([4, 16, 128]):
    for dt in [0.1]:
    #for dt in [0.1, 0.01]:
        maxT = 30
        Nt = math.floor(maxT / dt)
        t = np.linspace(0, maxT, num=Nt)

        v = np.zeros((Nt, N))
        if BOUNDARY_TYPE == 1:
            x = np.zeros((Nt, N))
            x[0][math.floor(N/2-1)] = 1
        elif BOUNDARY_TYPE == 2:
            x = np.zeros((Nt, N))
            for i in range(N):
                x[0][i] = np.sin(np.pi * (i+1) / (N + 1))
        elif BOUNDARY_TYPE == 3:
            x = np.zeros((Nt, N))
            for i in range(N):
                x[0][i] = np.sin(np.pi * N/2 * (i+1) / (N + 1))
        else:
            print("invalid BOUNDARY_TYPE!")
            quit()

        # one euler step, prepare lastAccell
        lastAccell = np.zeros(N)
        for n in range(N):
            v[1][n] = v[0][n] + a(x[0], n) * dt
            x[1][n] = x[0][n] + v[1][n] * dt
            lastAccell[n] = a(x[1], n)

        # verlet integration
        for i in range(2, Nt):
            for n in range(N):
                x[i][n] = x[i-1][n] + v[i-1][n]*dt + 0.5 * lastAccell[n] * dt*dt

            for n in range(N):
                accell = a(x[i], n)
                v[i][n] = v[i-1][n] + 0.5*(lastAccell[n] + accell)*dt
                lastAccell[n] = accell

    # plot for last dt
    toPlot = [0, 1, 2, 3]
    for i in toPlot:
        plt.plot(t, x[:,i], label="x_{}(t)".format(i))

    # shift, leave out last column
    energy = np.sum(np.power(x[:,:-1] - np.roll(x, -1, axis=1)[:,:-1], 2), axis=1)/2 + \
        np.sum(np.power(v, 2), axis=1)/2
    plt.plot(t, energy, label="Energy")

    plt.legend()
    plt.xlabel("t")
    plt.ylabel("x(t), E(t)")
    plt.grid()
    plt.xlim(0, 10)
    plt.savefig("ex4_N={}_dt={}_boundaries={}_curves.pdf".format(N, dt, BOUNDARY_TYPE))
    #plt.show()
    plt.close()

    x_results.append((N,x))

# Raster plots
f, ax = plt.subplots(3,1)
for i, (N,x) in enumerate(x_results):
    ax[i].imshow(x, cmap="gray", interpolation="nearest", aspect=N/Nt)
    ax[i].set_xlabel("oscillator")
    ax[i].set_ylabel("time")

plt.savefig("ex4_dt={}_boundaries={}_raster.pdf".format(dt, BOUNDARY_TYPE))
plt.show()
plt.close()
