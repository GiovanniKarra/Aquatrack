import path3d as p3d             # sert à construire le circuit numériquement
import numpy as np               # sert à faire des calcules plus efficacement, surtout avec les vecteurs et matrices
import matplotlib.pyplot as plt  # sert à faire des graphiques avec nos données


"""


SIMULATION


"""


# Paramètres physiques:
g = 9.81                          # accélération de gravitation [m/s**2]
b = 0.01                          # écart des rails [m]
r = 0.0079375                     # rayon de la bille [m]
h = np.sqrt(r ** 2 - b ** 2 / 4)  # hauteur du centre de la bille sur les rails [m]

e1 = 0.000575  # coefficient de frottement linéaire [m/(m/s)]

# charge les points de passage en tant que matrice (3 x nombre de points)
xyzPoints = np.loadtxt("points_de_passage.txt", unpack=True)

# réalise une interpolation des points de passage pour donner...
# sPath : array des abscisses curvilignes
# xyzPath : coordonnées en 3 dimensions
# tPath : vecteur unitaire tangent à chaque point
# cPath : vecteur unitaire normal à chaque point
sPath, xyzPath, tPath, cPath = p3d.path(xyzPoints)

# paramètres pour la simulation:
tEnd = 6  # durée de la simulation [s]
dt = 0.1  # pas de la simulation [s]

steps = int(tEnd / dt)       # nombre de pas de la simulation
tSim = np.zeros(steps + 1)   # temps: array[steps+1] * [s]
sSim = np.zeros(steps + 1)   # distance curviligne: array[steps+1] * [m]
VsSim = np.zeros(steps + 1)  # vitesse tangentielle: array[steps+1] * [m/s]
AsSim = np.zeros(steps + 1)  # acceleration curviligne: array[steps+1] * [m/s**2]
zSim = np.zeros(steps + 1)   # coordonnée z des points de la simulation

# valeurs initiales:
tSim[0] = 0
sSim[0] = 0
VsSim[0] = 0

# matrices (3 x steps) qui stockent les composantes des coordonnées,
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

    tan = path[1]                          # vecteur tangent
    norm = path[2]                         # vecteur normal

    gs = -g * tan[2]                       # norme de l'accélération gravitationnelle tangentielle
    gs_vector = gs * tan                   # vecteur g_s
    gn = np.array((0, 0, -g)) - gs_vector  # vecteur de l'accélération gravitationnelle normale
    Gn = VsSim[i] ** 2 * norm - gn         # accélération normale

    As = (gs - e1 * (VsSim[i] / h) * np.linalg.norm(Gn)) / M  # accélération curviligne

    AsSim[i] = As
    VsSim[i + 1] = VsSim[i] + As * dt          # on varie la vitesse curviligne suivante selon l'accélération
    sSim[i + 1] = sSim[i] + VsSim[i + 1] * dt  # on varie la position curviligne suivante selon la vitesse
    tSim[i + 1] = tSim[i] + dt                 # on varie le temps t suivant selon dt

    # variables utilisées pour le graphique
    xyz = p3d.ainterp(sSim[i], sPath, xyzPath)  # coordonnées du point simulé
    xyzMarks[:, i] = xyz                        # stockage des coordonnées
    cMarks[:, i] = Gn                           # stockage du vecteur de l'accélération normale
    tMarks[:, i] = As * tan                     # stockage du vecteur de l'accélération curviligne

    zSim[i] = xyz[2]                            # stockage des coordonnées z pour le calcul de l'énergie potentielle

EpSim = g * (zSim + 0.070)    # énergie potentielle spécifique [m**2/s**2]
EkSim = 0.5 * M * VsSim ** 2  # énergie cinétique spécifique [m**2/s**2]


"""


DESSIN DES DIFFÉRENTS GRAPHIQUES


"""


# plot vitesses, positions, accelerations

plt.figure()                            # création d'un plot

plt.subplot(311)                        #
plt.plot(tSim, sSim, 'b-')              # DESSIN DU GRAPHIQUE
plt.ylabel('abscisse curviligne [m]')   # DES POSITIONS CURVILIGNES
plt.xlabel('t [s]')                     #

plt.subplot(312)                        #
plt.plot(tSim, VsSim, 'b-')             # DESSIN DU GRAPHIQUE
plt.ylabel('vitesse [m/s]')             # DES VITESSES
plt.xlabel('t [s]')                     #

plt.subplot(313)                        #
plt.plot(tSim, AsSim, 'b-')             # DESSIN DU GRAPHIQUE
plt.ylabel('acceleration [m/s^2]')      # DES ACCÉLÉRATIONS
plt.xlabel('t [s]')                     #

plt.show()                              # affichage du plot sur l'écran


# plot énergies

plt.figure()                                      # création d'un plot

plt.plot(tSim, EpSim, 'b-', label='Ep/m')         #
plt.plot(tSim, EkSim, 'r-', label='Ek/m')         # DESSIN DU GRAPHIQUE DE L'ÉNERGIE MÉCANIQUE
plt.plot(tSim, EkSim + EpSim, 'k-', label='E/m')  #

plt.legend()                                      #
plt.ylabel('Energy/mass [J/kg]')                  # CRÉATION ET AFFICHAGE DE LA LÉGENDE
plt.xlabel('t [s]')                               #

plt.show()                                        # affichage du plot sur l'écran


# plot circuit

fig = plt.figure()                                                  # création d'un plot

ax = fig.add_subplot(projection='3d')                               #
ax.set_box_aspect(np.ptp(xyzPath, axis=1))                          #
ax.plot(xyzPoints[0], xyzPoints[1], xyzPoints[2],                   #
        'bo', label='points', markersize=1)                         #
ax.plot(xyzPath[0], xyzPath[1], xyzPath[2],                         #
        'k-', lw=0.5, label='path')                                 # DESSIN DU CIRCUIT EN 3D
scale = 0.2 * sPath[-1] / steps                                     # AVEC LES VECTEURS D'ACCÉLÉRATION CURVILIGNE
ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],                    # ET ACCÉLÉRATION NORMALE
          scale * tMarks[0], scale * tMarks[1],                     #
          scale * tMarks[2], color='r', linewidth=0.5, label='As')  #
ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],                    #
          scale * cMarks[0], scale * cMarks[1], scale * cMarks[2],  #
          color='g', linewidth=0.5, label='Gn')                     #

ax.legend()                                                         # affichage de la légende

plt.show()                                                          # affichage du plot
