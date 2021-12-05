import path3d as p3d
import numpy as np
import matplotlib.pyplot as plt


"""


SIMULATION


"""


# Paramètres physiques:
g = 9.81  # accélération de gravitation [m/s**2]
b = 0.01  # écart des rails [m]
r = 0.0079375  # rayon de la bille [m]
h = np.sqrt(r ** 2 - b ** 2 / 4)  # hauteur du centre de la bille sur les rails [m]
m = 0.016  # masse de la bille

e1 = 0.000575  # coefficient de frottement linéaire [m/(m/s)]

# chemin de la bille (et autres paramètres)
# xyzPoints = np.loadtxt("looping_points.txt", unpack=True)
xyzPoints = np.loadtxt("points_de_passage.txt", unpack=True)
sPath, xyzPath, tPath, cPath = p3d.path(xyzPoints)

# paramètres pour la simulation:
tEnd = 6  # durée de la simulation [s]
dt = 0.1  # pas de la simulation [s]

steps = int(tEnd / dt)  # nombre de pas de la simulation
tSim = np.zeros(steps + 1)  # temps: array[steps+1] * [s]
sSim = np.zeros(steps + 1)  # distance curviligne: array[steps+1] * [m]
VsSim = np.zeros(steps + 1)  # vitesse tangentielle: array[steps+1] * [m/s]
AsSim = np.zeros(steps + 1)  # acceleration curviligne: array[steps] * [m/s^2]
zSim = np.zeros(steps + 1)  # coordonnée z des points de la simulation

# valeurs initiales:
tSim[0] = 0
sSim[0] = 0
VsSim[0] = 0

# matries 3 x steps qui stockent les composantes des coordonnées,
# vecteurs tangents et vecteurs normales des points de la simulation
xyzMarks = np.empty((3, steps))
tMarks = np.empty((3, steps))
cMarks = np.empty((3, steps))

M = 1 + 2 / 5 * r ** 2 / h ** 2  # coefficient d'inertie [1]

# boucle de simulation:
for i in range(steps):
    # tuple sous la forme (coordonnées, vecteur tangent, vecteur normal)
    # en fonction de l'abscisse curviligne du circuit
    path = p3d.path_at(sSim[i], (sPath, xyzPath, tPath, cPath))

    tan = path[1]  # vecteur tangent
    norm = path[2]  # vecteur normal

    gs = -g * tan[2]  # norme de l'accélération gravitationnelle tangentielle
    gs_vector = gs * tan  # vecteur g_s
    gn = np.array((0, 0, -g)) - gs_vector  # vecteur de l'accélération gravitationnelle normale
    Gn = VsSim[i] ** 2 * norm - gn

    As = (gs - e1 * (VsSim[i] / h) * np.linalg.norm(Gn)) / M  # accélération curviligne

    AsSim[i] = As
    VsSim[i + 1] = VsSim[i] + As * dt  # on varie la vitesse curviligne suivante selon l'accélération
    sSim[i + 1] = sSim[i] + VsSim[i + 1] * dt  # on varie la position curviligne suivante selon la vitesse
    tSim[i + 1] = tSim[i] + dt  # on varie le temps t suivant selon dt

    # variables utilisés pour le graphique (coordonnées, vecteurs, etc.)
    xyz = p3d.ainterp(sSim[i], sPath, xyzPath)
    xyzMarks[:, i] = xyz
    cMarks[:, i] = Gn
    tMarks[:, i] = gs_vector

    zSim[i] = xyz[2]

    # si la bille dépasse le circuit, on met la même valeur aux itérations restantes
    if sSim[i] >= sPath[-1]:
        print(tSim[i])
        sSim[i:] = sSim[i]
        VsSim[i:] = VsSim[i]
        tSim[i:] = tSim[i]

EpSim = g * (zSim + 0.070)  # énergie potentielle spécifique [m**2/s**2]
EkSim = 0.5 * M * VsSim ** 2  # énergie cinétique spécifique [m**2/s**2]


"""


DESSIN DES DIFFÉRENTS GRAPHIQUES


"""


# plot vitesses, positions, accelerations
plt.figure()
plt.subplot(311)
plt.plot(tSim, sSim, 'b-')
plt.ylabel('abscisse curviligne [m]')
plt.xlabel('t [s]')
plt.subplot(312)
plt.plot(tSim, VsSim, 'b-')
plt.ylabel('vitesse tangentielle [m/s]')
plt.xlabel('t [s]')
plt.subplot(313)
plt.plot(tSim, AsSim, 'b-')
plt.ylabel('acceleration tangentielle [m/s^2]')
plt.xlabel('t [s]')
plt.show()

# plot énergies
plt.figure()
plt.plot(tSim, EpSim, 'b-', label='Ep/m')
plt.plot(tSim, EkSim, 'r-', label='Ek/m')
plt.plot(tSim, EkSim + EpSim, 'k-', label='E/m')
plt.legend()
plt.ylabel('Energy/mass [J/kg]')
plt.xlabel('t [s]')
plt.show()

# plot circuit
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_box_aspect(np.ptp(xyzPath, axis=1))
ax.plot(xyzPoints[0], xyzPoints[1], xyzPoints[2], 'bo', label='points', markersize=1)
ax.plot(xyzPath[0], xyzPath[1], xyzPath[2], 'k-', lw=0.5, label='path')
scale = 0.2 * sPath[-1] / steps
ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],
          scale * tMarks[0], scale * tMarks[1], scale * tMarks[2],
          color='r', linewidth=0.5, label='gs')
ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],
          scale * cMarks[0], scale * cMarks[1], scale * cMarks[2],
          color='g', linewidth=0.5, label='Gn')
ax.legend()
plt.show()
