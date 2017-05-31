import math
import random

import numpy as np
import matplotlib.pyplot as plt

PERIODIC_BOUNDARIES = False

#N = 10
#N = 100
N = 1000

N_SAMPLES = 1000
N_SAMPLES = 10000

def energy(spins):
    # shift elements
    shifted = np.roll(spins, -1)
    if not PERIODIC_BOUNDARIES:
        shifted[-1] = 0
    return -np.dot(spins, shifted)

def U_theory(N, beta):
    return -(N - 1)/N * np.tanh(beta)

def C_theory(N, beta):
    x = beta / np.cosh(beta)
    return (N - 1)/N * x*x

def mmc(beta):
    U = 0
    C = 0
    spins = np.random.randint(2, size=N) * 2 - 1 # N random numbers, either -1 or 1
    e = energy(spins)
    N_WAIT = int(N_SAMPLES / 10)
    for sample in range(N_SAMPLES + N_WAIT):
        for i in range(N):
            j = random.randint(0, N-1)

            #e_old = rest + spins[j] * (spins[j-1] + spins[j+1])
            # e_new = rest - spins[j] * (spins[j-1] + spins[j+1])
            # delta_e = e_new - e_old = 2 * spins[j] * (spins[j-1] + spins[j+1])

            # assuming spins[j] is flipped
            delta_e = 0
            if j > 0: delta_e += spins[j-1]
            if j < N-1: delta_e += spins[j+1]
            delta_e *= 2.0 * spins[j]

            q = math.exp(-delta_e * beta)
            if q > random.random():
                # flip the spin
                spins[j] = -spins[j]
                e += delta_e

        if sample >= N_WAIT:
            U += e
            C += e*e
    U = U / N_SAMPLES
    C = beta*beta * (C / N_SAMPLES - U*U)

    return U, C

print("T           beta         U/N          C/N         U_t          C_t         err_U       err_C")
for T in np.arange(4.0, 0.19, -0.2):
    k_B = 1
    beta = 1 / (k_B * T)

    U, C = mmc(beta)

    U_t = U_theory(N, beta)
    C_t = C_theory(N, beta)

    print("{0:.3E}   {1:.3E}    {2:.3E}   {3:.3E}   {4:.3E}   {5:.3E}   {6:.3E}   {7:.3E}".format(
        T, beta, U/N, C/N, U_t, C_t, np.abs((U_t - U/N)/U_t), np.abs((C_t - C/N)/C_t)))
