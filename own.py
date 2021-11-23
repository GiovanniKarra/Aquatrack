import matplotlib.pyplot as plt
import numpy as np

import path3d as p3d

# Physical parameters:
g = 9.8  # Gravity acceleration [m/s**2]
b = 0.009  # Track gap [m]
r = 0.0079375  # B [m]
h = np.sqrt(r ** 2 - b ** 2 / 4)  # hauteur du centre de la bille sur les rails [m]

e = 0.00045  # coefficient de frottement linéaire [m/(m/s)]

xyzPoints = np.loadtxt("points_de_passage.txt", unpack=True)

# Path & vectors
sPath, xyzPath, TPath, CPath = p3d.path(xyzPoints)

# Marks to display on graph
num = 38  # Marks amount
length = sPath[-1]
sMarks = np.linspace(0, length, num)
xyzMarks = np.empty((3, num))  # Coordinates
TMarks = np.empty((3, num))  # Tangent vector
CMarks = np.empty((3, num))  # Curvature vector

# Simulation parameters:
tEnd = 20  # Simulation duration [s]
steps = 38  # Steps amount
dt = tEnd / steps  # Simulation step duration [s]

tSim = np.zeros(steps + 1)  # Time: array[steps+1] * [s]
sSim = np.zeros(steps + 1)  # Curvilinear distance: array[steps+1] * [m]
VsSim = np.zeros(steps + 1)  # Tangential speed: array[steps+1] * [m/s]

M = 1 + 2 / 5 * r ** 2 / h ** 2  # Inertie coefficient

# Initial values:
tSim[0] = 0
sSim[0] = 0
VsSim[0] = 0

# Simulation loop
for i in range(steps):
    Tz = TPath[2, i]

    # Gravity vectors
    g_s = -g * Tz
    gs_vector = g_s * TPath[:, i]
    g_n = np.array([0, 0, -g]) - gs_vector
    G_n = np.linalg.norm(VsSim[i] ** 2 * CPath[:, i] - g_n)

    # Acceleration
    As = (g_s - e * VsSim[i] / h * G_n) / M

    # Register next step speed, curvilinear vector & time
    VsSim[i + 1] = VsSim[i] + As * dt
    sSim[i + 1] = sSim[i] + VsSim[i + 1] * dt
    tSim[i + 1] = tSim[i] + dt

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
scaleC = scale*10
ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],
          scale * TMarks[0], scale * TMarks[1], scale * TMarks[2],
          color='r', linewidth=0.5, label='T')
ax.quiver(xyzMarks[0], xyzMarks[1], xyzMarks[2],
          scaleC * CMarks[0], scaleC * CMarks[1], scaleC * CMarks[2],
          color='g', linewidth=0.5, label='C')
print(xyzMarks)
print(CMarks)
ax.legend()
plt.show()
