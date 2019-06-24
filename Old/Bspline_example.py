'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps
import scipy.interpolate as spi

# plot cubic cardinal B-spline (knots 0, 1, 2, 3, 4)
# plot cubic cardinal B-spline (knots 0, 1, 2, 3, 4)
p = 4
xx = np.linspace(0, p+1, 100)
yy = sps.bspline(xx - (p+1)/2, p)
plt.plot(xx, yy)
plt.show()

# plot cubic non-uniform spline (m=5 DOFs)
xi = [0, 1, 3, 4, 6, 7, 8, 10, 11]
c = [2, -1, 1, 0, 1]
s = spi.BSpline(xi, c, p)
m = len(c)
xx = np.linspace(xi[p], xi[m])
yy = s(xx)
plt.plot(xx, yy)
plt.show()

'''


import numpy as np
import matplotlib.pyplot as plt
import math

x, y = np.genfromtxt('data', unpack=True, skip_header=1)
# find lots of points on the piecewise linear curve defined by x and y
M = 1000
t = np.linspace(0, len(x), M)
x = np.interp(t, np.arange(len(x)), x)
y = np.interp(t, np.arange(len(y)), y)
tol = 1.5
i, idx = 0, [0]
while i < len(x):
    total_dist = 0
    for j in range(i+1, len(x)):
        total_dist += math.sqrt((x[j]-x[j-1])**2 + (y[j]-y[j-1])**2)
        if total_dist > tol:
            idx.append(j)
            break
    i = j+1

xn = x[idx]
yn = y[idx]
fig, ax = plt.subplots()
ax.plot(x, y, '-')
ax.scatter(xn, yn, s=50)
ax.set_aspect('equal')
plt.show()