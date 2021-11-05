import path3d as p3d
import numpy as np
import matplotlib.pyplot as plt

# from mpl_toolkits.mplot3d import Axes3D

# points de passage
xyzPoints = p3d.looping_points()

# sauvetage des points dans un fichier
np.savetxt('looping_points.txt', xyzPoints.T, fmt='%10.5f')

# chemin et vecteurs
sPath, xyzPath, TPath, CPath = p3d.path(xyzPoints)

# points jalons à afficher sur le graphique
num = 30  # nombre de jalons
length = sPath[-1]
sMarks = np.linspace(0, length, num)
xyzMarks = np.empty((3, num))  # coordonnées
TMarks = np.empty((3, num))  # vecteur tangent
CMarks = np.empty((3, num))  # vecteur de courbure

for i in range(num):
    xyz = p3d.ainterp(sMarks[i], sPath, xyzPath)
    T = p3d.ainterp(sMarks[i], sPath, TPath)
    C = p3d.ainterp(sMarks[i], sPath, CPath)

    xyzMarks[:, i] = xyz
    CMarks[:, i] = C
    TMarks[:, i] = T

# graphique 3D : points, chemin et vecteurs
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_box_aspect(np.ptp(xyzPath, axis=1))
ax.plot(xyzPoints[0], xyzPoints[1], xyzPoints[2], 'bo', label='points')
ax.plot(xyzPath[0], xyzPath[1], xyzPath[2], 'k-', lw=0.5, label='path')
scale = 0.5 * length / num
ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],
          scale * TMarks[0], scale * TMarks[1], scale * TMarks[2],
          color='r', linewidth=0.5, label='T')
ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],
          scale * CMarks[0], scale * CMarks[1], scale * CMarks[2],
          color='g', linewidth=0.5, label='C')
ax.legend()
plt.show()
