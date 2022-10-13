
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import os


def phenotype(pA, pB):
    qA = 1 - pA
    qB = 1 - pB
    mut = qA ** 2 + qB ** 2 - (qA ** 2 * qB ** 2)
    WT = 1 - mut
    return WT


freqA = np.linspace(0, 1, num=101)
freqB = np.linspace(0, 1, num=101)
X, Y = np.meshgrid(freqA, freqB)
zs = phenotype(X, Y)
Z = zs.reshape(X.shape)

# Normalize to [0,1]
Z = (Z - Z.min()) / (Z.max() - Z.min())

colors = cm.Greys(Z)
# rcount, ccount, _ = colors.shape

fig = plt.figure(figsize=(8, 6))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=4, cstride=4,
                       facecolors=colors, shade=False)
surf.set_edgecolor((0, 0, 0, 0.8))

m = cm.ScalarMappable(cmap=cm.Greys)
m.set_array(Z)
plt.colorbar(m, fraction=0.03, pad=0.04)

plt.show()
