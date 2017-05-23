import numpy as np
import matplotlib.pyplot as plt

N_PART = 10000
N_JUMPS = 1000

part_pos = np.zeros(N_PART, dtype=np.int32)

# Matricle number 309231
np.random.seed(9231)

plt.rc('text', usetex=True)

var = np.zeros(N_JUMPS)
for i_jump in range(N_JUMPS):
    for i_part in range(N_PART):
        if np.random.random() < 0.5:
            part_pos[i_part] += 1
        else:
            part_pos[i_part] -= 1
    var[i_jump] = np.var(part_pos)

    if (i_jump+1) % 250 == 0:
        plt.hist(part_pos, bins=50, alpha=0.5,
            label="{} jumps".format(i_jump+1),
            histtype="stepfilled", range=(-150, 150))
plt.legend()
plt.xlabel("Position")
plt.ylabel("Number of particles")
plt.savefig("walk_hist.pdf")
plt.show()

plt.plot(var)
plt.xlabel("Number of steps $N$")
plt.ylabel("Variance of position $<x^2> - <x>^2$")
plt.savefig("walk_var.pdf")
plt.show()
