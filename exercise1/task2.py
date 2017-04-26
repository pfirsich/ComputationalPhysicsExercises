import numpy as np
import matplotlib.pyplot as plt

def cheby(x, N):
	ret = np.zeros((N+1, len(x)))
	ret[0] = np.ones(len(x))
	ret[1] = np.copy(x)
	for n in range(2,N+1):
		ret[n] = 2*x*ret[n-1] - ret[n-2]
	return ret
		
x = np.linspace(-1, 1, 100)
T = cheby(x, 4)
for i in range(len(T)):
	plt.plot(x, T[i], label="T_" + str(i))
plt.legend()
plt.grid(True)
plt.xlabel("x")
plt.ylabel("T_n(x)")
plt.show()
