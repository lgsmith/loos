import numpy as np
import matplotlib.pyplot as plt

# set this equal to argv[1] if you want this to have a crude CLI
filename = 'testmat.asc'
mat = np.loadtxt(filename)
# matrix is diagonal, so get upper triangular view
tri = np.triu(mat, 1)
# density, bins = np.histogram(tri, density=True)
# binw = bins[1] - bins[0]
# cdf = np.cumsum(density) * binw

fig, ax = plt.subplots()
density, bins, patches = ax.hist(tri.flat, density=True)
binws = np.diff(bins)
cdf = np.cumsum(density * binws)
cdf_xvals = bins[:-1] + binws/2
ax.plot(cdf_xvals, cdf)
plt.show()
