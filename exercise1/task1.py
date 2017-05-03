import numpy as np

# Sorry for not commenting these properly
def np_mat2tex(m):
    return " \\\\ \n".join(" & ".join("{:.4f}".format(elem)
        for elem in row) for row in m.tolist())

def np_vec2tex(m, delim=" \\\\ "):
    return delim.join("{:.4f}".format(elem)
        for elem in m.tolist())

# 1.
# 309231
np.random.seed(9231)

# Rescale random values from [0,1]
A = np.random.rand(6,6) * 10.0 - 5.0
print("A = ", np_mat2tex(A), "\n")

# 2.
max_value = A.max()
max_index = np.where(A == max_value)
print(max_index, max_value, "\n")

# 3.
row_max = A.max(axis=0)
col_max = A.max(axis=1)
print("r = ", np_vec2tex(row_max, " & "))
print("c = ", np_vec2tex(col_max), "\n")

# I hope row_max * col_max is fine (instead of the other way round)
print(np.dot(row_max, col_max))

# 4.
B = np.random.rand(6,6) * 10.0 - 5.0
C = np.dot(A, B)
D = np.dot(B, A)

print("B = ", np_mat2tex(B))
print("C = ", np_mat2tex(C))
print("D = ", np_mat2tex(D))

