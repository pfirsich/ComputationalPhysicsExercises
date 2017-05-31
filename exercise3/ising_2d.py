import math
import random
from multiprocessing import Pool
from functools import partial
import json

import numpy as np

PERIODIC_BOUNDARIES = False

N_SAMPLES = 1000
#N_SAMPLES = 10000

def energy(spins):
    # shift elements on y axis <=> spins[i][j] gets the value of spins[i+1][j]
    shiftedy = np.roll(spins, -1, 0)
    # spins[i][j] gets the value of spins[i][j+1]
    shiftedx = np.roll(spins, -1, 1)

    if not PERIODIC_BOUNDARIES:
        shiftedy[-1,:] = 0
        shiftedx[:,-1] = 0

    # multiply matrices component wise
    return - np.sum(np.multiply(spins, shiftedy)) - np.sum(np.multiply(spins, shiftedx))

def mmc(T, N):
    print("T =", T)
    k_B = 1
    beta = 1 / (k_B * T)

    U, C, M = 0, 0, 0

    # N random numbers, either -1 or 1
    spins = np.random.randint(2, size=(N,N)) * 2 - 1
    e = energy(spins)
    N_WAIT = int(N_SAMPLES / 10)
    N_MEASUREMENTS = 0
    for sample in range(N_SAMPLES + N_WAIT):
        for n in range(N*N):
            j1, j2 = random.randint(0, N-1), random.randint(0, N-1)

            # assuming spins[j1,j2] is flipped
            delta_e = 0
            if j1 > 0: delta_e += spins[j1-1,j2]
            if j2 > 0: delta_e += spins[j1,j2-1]
            if j1 < N-1: delta_e += spins[j1+1,j2]
            if j2 < N-1: delta_e += spins[j1,j2+1]
            delta_e *= 2.0 * spins[j1,j2]

            q = math.exp(-delta_e * beta)
            if q > random.random():
                # flip the spin
                spins[j1,j2] = -spins[j1,j2]
                e += delta_e

        if sample >= N_WAIT:
            N_MEASUREMENTS += 1
            U += e
            C += e*e
            M += abs(np.sum(spins))

    U = U / N_MEASUREMENTS
    C = beta*beta * (C / N_MEASUREMENTS - U*U)
    M = M / N_MEASUREMENTS

    return U, C, M

if __name__ == '__main__':
    jsonData = []

    Ts = np.arange(0.2, 4.1, 0.2)
    #for i_N, N in enumerate([10, 50, 100]):
    for i_N, N in enumerate([10, 50, 100]):
        print("N =", N)
        _mmc = partial(mmc, N=N)

        U, C, M = [], [], []
        with Pool(processes=8) as pool:
            for vU, vC, vM in pool.map(_mmc, Ts):
                U.append(vU)
                C.append(vC)
                M.append(vM)

        jsonData.append({'N': N, 'T': list(Ts), 'U': U, 'C': C, 'M': M})

    with open("ising_2d_data.json", "w") as outFile:
        json.dump(jsonData, outFile)