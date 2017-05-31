import json
import math

import numpy as np
import matplotlib.pyplot as plt

def M_theory(beta):
    beta_C = 1 / (2 / math.log(1 + math.sqrt(2)))
    if beta < beta_C:
        return 0
    else:
        x = math.sinh(2*beta)
        return math.pow(1 - 1/(x*x*x*x), 1/8)

plt.rc('text', usetex=True)
f, axarr = plt.subplots(3, 1, figsize=(15,10))

for i in range(3):
    axarr[i].set_xlabel(r"$T$")
axarr[0].set_ylabel(r"$U/N^2$")
axarr[1].set_ylabel(r"$C/N^2$")
axarr[2].set_ylabel(r"$M/N^2$")

with open("ising_2d_data.json") as inFile:
    jsonData = json.load(inFile)

for data_set in jsonData:
    N = data_set["N"]
    Ts = np.array(data_set["T"])
    U = np.array(data_set["U"])
    C = np.array(data_set["C"])
    M = np.array(data_set["M"])

    label = "N = {}".format(N)
    axarr[0].plot(Ts, U/N/N, label=label)
    axarr[1].plot(Ts, C/N/N, label=label)
    axarr[2].plot(Ts, M/N/N, label=label)

axarr[2].plot(Ts, np.vectorize(M_theory)(1/Ts), label="Theory")

for i in range(3):
    axarr[i].legend()

plt.savefig("ising_2d.pdf")
plt.show()