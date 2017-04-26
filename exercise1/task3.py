import numpy as np

i = np.arange(1, dtype=np.uint8)
for j in range(0, 300):
	i[0] = i[0] + 1
	print(i[0])
