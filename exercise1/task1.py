import numpy as np 

# 1.
# 309231
np.random.seed(9231)

# Rescale random values from [0,1]
A = np.random.rand(6,6) * 10.0 - 5.0
print(A)

# 2.
max_value = A.max()
max_index = np.where(A == max_value)
print(max_index, max_value)

# 3.
row_max = A.max(axis=0)
col_max = A.max(axis=1)
print(row_max, col_max)

# I hope row_max * col_max is fine (instead of the other way round)
print(np.dot(row_max, col_max))

# 4.
B = np.random.rand(6,6) * 10.0 - 5.0
C = np.dot(A, B)
D = np.dot(B, A)

print(B)
print(C)
print(D)

