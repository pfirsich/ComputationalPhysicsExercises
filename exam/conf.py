Delta = 0.025
tau = 0.00025
m = 40000
xMin, xMax = -15, 15
xRange = xMax - xMin
L = (xRange / Delta) + 1
assert abs(L - int(L)) < 1e-8
L = int(L)