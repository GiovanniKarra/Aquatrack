import path3d as p3d
import numpy as np
import matplotlib.pyplot as plt


# Paramètres physiques:
g = 9.81  # accélération de gravitation [m/s**2]
b = 0.01  # écart des rails [m]
r = 0.0079375  # rayon de la bille [m]
h = np.sqrt(r ** 2 - b ** 2 / 4)  # hauteur du centre de la bille sur les rails [m]

# e1 = 0.001270193299685568  # coefficient de frottement linéaire [m/(m/s)]
# e1 = 0.000995637359264194
# e1 = 0.00024985386042669975
e1 = 0.001

# chemin de la bille (et autres paramètres)
xyzPoints = np.loadtxt("looping_points.txt", unpack=True)
sPath, xyzPath, tPath, cPath = p3d.path(xyzPoints)

# paramètres pour la simulation:
tEnd = 10  # durée de la simulation [s]
dt = 0.033333333333333326  # pas de la simulation [s]

steps = int(tEnd / dt)  # nombre de pas de la simulation
tSim = np.zeros(steps + 1)  # temps: array[steps+1] * [s]
sSim = np.zeros(steps + 1)  # distance curviligne: array[steps+1] * [m]
VsSim = np.zeros(steps + 1)  # vitesse tangentielle: array[steps+1] * [m/s]

# valeurs initiales:
tSim[0] = 0
sSim[0] = 0
VsSim[0] = 0

# boucle de simulation:
for i in range(steps):
    xyz = p3d.ainterp(sSim[i], sPath, xyzPath)
    # p = 2 * A * x  # pente dz/dx
    # cos_beta = 1 / np.sqrt(1 + p * p)
    # sin_beta = p / np.sqrt(1 + p * p)
    # c = 2 * A / (1 + p * p) ** 1.5  # courbure

    # As = (-g * sin_beta - e1 * VsSim[i] / h * (g * cos_beta + c * VsSim[i] ** 2)) / M

    # VsSim[i + 1] = VsSim[i] + As * dt
    # sSim[i + 1] = sSim[i] + VsSim[i + 1] * dt
    # tSim[i + 1] = tSim[i] + dt
