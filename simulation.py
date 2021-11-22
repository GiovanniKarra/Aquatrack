import path3d as p3d
import numpy as np
import matplotlib.pyplot as plt


# Paramètres physiques:
g = 9.81  # accélération de gravitation [m/s**2]
b = 0.01  # écart des rails [m]
r = 0.0079375  # rayon de la bille [m]
h = np.sqrt(r ** 2 - b ** 2 / 4)  # hauteur du centre de la bille sur les rails [m]
m = 0.016  # masse de la bille

# e1 = 0.001270193299685568  # coefficient de frottement linéaire [m/(m/s)]
# e1 = 0.000995637359264194
# e1 = 0.00024985386042669975
e1 = 0.01

# chemin de la bille (et autres paramètres)
# xyzPoints = np.loadtxt("looping_points.txt", unpack=True)
xyzPoints = np.loadtxt("points_de_passage.txt", unpack=True)
sPath, xyzPath, tPath, cPath = p3d.path(xyzPoints)

# paramètres pour la simulation:
tEnd = 100  # durée de la simulation [s]
dt = 0.3333333333333  # pas de la simulation [s]

steps = int(tEnd / dt)  # nombre de pas de la simulation
tSim = np.zeros(steps + 1)  # temps: array[steps+1] * [s]
sSim = np.zeros(steps + 1)  # distance curviligne: array[steps+1] * [m]
VsSim = np.zeros(steps + 1)  # vitesse tangentielle: array[steps+1] * [m/s]
zSim = np.zeros(steps + 1)

# valeurs initiales:
tSim[0] = 0
sSim[0] = 0
VsSim[0] = 0

xyzMarks = np.empty((3, steps))
tMarks = np.empty((3, steps))
cMarks = np.empty((3, steps))

# boucle de simulation:
for i in range(steps):
    path = p3d.path_at(sSim[i], (sPath, xyzPath, tPath, cPath))

    tan = path[1]
    norm = path[2]

    gs = (-g * tan[2]) * tan
    gn = np.array((0, 0, -g)) - gs
    Gn = VsSim[i] ** 2 * norm - gn

    As = (np.linalg.norm(gs) - e1 * (VsSim[i]/h) * np.linalg.norm(Gn)) / (1 + 2/5 * (r ** 2 / h ** 2))

    VsSim[i + 1] = VsSim[i] + As * dt
    sSim[i + 1] = sSim[i] + VsSim[i + 1] * dt
    tSim[i + 1] = tSim[i] + dt

    xyz = p3d.ainterp(sSim[i], sPath, xyzPath)
    xyzMarks[:, i] = xyz
    cMarks[:, i] = Gn
    tMarks[:, i] = gs

    zSim[i] = xyz[2]

M = 1 + 2 / 5 * r ** 2 / h ** 2  # coefficient d'inertie [1]

EpSim = g * zSim  # énergie potentielle spécifique [m**2/s**2]
EkSim = 0.5 * M * VsSim ** 2  # énergie cinétique spécifique [m**2/s**2]

# plot énergies
plt.figure()
plt.plot(tSim, EpSim, 'b-', label='Ep/m')
plt.plot(tSim, EkSim, 'r-', label='Ek/m')
plt.plot(tSim, EkSim + EpSim, 'k-', label='E/m')
plt.legend()
plt.ylabel('Energy/mass [J/kg]')
plt.xlabel('t [s]')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_box_aspect(np.ptp(xyzPath, axis=1))
ax.plot(xyzPoints[0], xyzPoints[1], xyzPoints[2], 'bo', label='points')
ax.plot(xyzPath[0], xyzPath[1], xyzPath[2], 'k-', lw=0.5, label='path')
scale = 2 * sPath[-1] / steps
ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],
          scale * tMarks[0], scale * tMarks[1], scale * tMarks[2],
          color='r', linewidth=0.5, label='gs')
ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],
          scale * cMarks[0], scale * cMarks[1], scale * cMarks[2],
          color='g', linewidth=0.5, label='Gn')
ax.legend()
plt.show()
