import math

import numpy as np
import matplotlib.pyplot as plt

LAMBDA = 1
NUM_GRID_POINTS_PER_LAMBDA = 50
DELTA = LAMBDA / NUM_GRID_POINTS_PER_LAMBDA # spatial resolution, 0.02
TAU = 0.9 * DELTA # temporal resolution
#TAU = 1.05 * DELTA # temporal resolution
X = 100*LAMBDA # Length of simulation box in spatial units
L = math.floor(X / DELTA) # Length of simulation box in index units, 5000
F = 1
NUM_TIMESTEPS = 10000
GLASS_THICKNESS = 2*LAMBDA
#GLASS_THICKNESS = 50*LAMBDA
N_GLASS = 1.46
SOURCE_POS = 20*LAMBDA
SOURCE_POS_I = math.floor(SOURCE_POS / DELTA) # 1000
BOUNDARY_THICKNESS = 6*LAMBDA

# turn on and off slowly to get a wave packet
def source(t):
    e = (t-30)/10
    return math.sin(2*math.pi*F*t) * math.exp(-e*e)

def sigma(x):
    if (x >= 0 and x <= BOUNDARY_THICKNESS) or (x >= X-BOUNDARY_THICKNESS and x <= X):
        return 1
    else:
        return 0

def eps(x):
    if x >= X/2 and x < X/2 + GLASS_THICKNESS:
        return N_GLASS*N_GLASS
    else:
        return 1

def mu(x):
    return 1

C = np.zeros(L+1)
D = np.zeros(L+1)

for l in range(L+1):
    x = l*DELTA
    temp = sigma(x)*TAU / (2*eps(x))
    C[l] = (1 - temp) / (1 + temp)
    D[l] = TAU / eps(x) / (1 + temp)

A = np.zeros(L)
B = np.zeros(L)

for l in range(L):
    x = (l + 0.5) * DELTA
    temp = sigma(x)*TAU / (2*mu(x))
    A[l] = (1 - temp) / (1 + temp)
    B[l] = TAU / mu(x) / (1 + temp)

E = np.zeros(L+1) # E_z
H = np.zeros(L) # H_y

PLOT_T = [1000, 2500, 5000, 6000]

plt.rc('text', usetex=True)

for t in range(NUM_TIMESTEPS):
    E[1:L-1] = D[1:L-1] * (H[1:L-1] - H[0:L-2]) / DELTA + C[1:L-1] * E[1:L-1]
    E[SOURCE_POS_I] -= D[SOURCE_POS_I] * source(t*TAU)
    H[0:L-1] = B[0:L-1] * (E[1:L] - E[0:L-1]) / DELTA + A[0:L-1] * H[0:L-1]

    if t in PLOT_T:
        E_max = np.max(np.abs(E[:L/2]))
        print("E_max(t=" + str(t) + ") =", E_max)

        ax = plt.gca()
        plt.plot(E)
        plt.title("t = " + str(t))
        yRange = 0.016
        plt.ylim(-yRange, yRange)
        ax.fill_between([0, BOUNDARY_THICKNESS/DELTA], -yRange, yRange, facecolor='grey', alpha=0.5)
        ax.fill_between([(X-BOUNDARY_THICKNESS)/DELTA, X/DELTA], -yRange, yRange, facecolor='grey', alpha=0.5)
        ax.fill_between([X/2/DELTA, (X/2+GLASS_THICKNESS)/DELTA], -yRange, yRange, facecolor='green', alpha=0.5)
        plt.grid(True)
        plt.xlabel("$x / \Delta$")
        plt.ylabel("$E$")
        plt.savefig("t=" + str(t) + ",glass=" + str(math.floor(GLASS_THICKNESS/LAMBDA)) + ".pdf")
        #plt.show()
        plt.close()