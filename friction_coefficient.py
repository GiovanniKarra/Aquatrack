import numpy as np

tEnd = 18.8
dt = 0.033333333333333326
steps = int(tEnd / dt)

L = 0.935
H = 0.375
A = 4 * H / L ** 2

g = 9.81
b = 0.01
r = 0.0079375
h = np.sqrt(r ** 2 - b ** 2 / 4)

xPath = np.linspace(-L / 2, L / 2, steps)
zPath = A * xPath ** 2
dsHalf = np.sqrt(np.diff(xPath) ** 2 + np.diff(zPath) ** 2)
sPath = np.hstack((0, np.cumsum(dsHalf)))

tExp, xExp, VxExp, VyExp, asExp = np.loadtxt('exp_data_example.txt', unpack=True, skiprows=1)

VsExp = np.zeros(steps)
sExp = np.zeros(steps)

for i in range(steps):
    VsExp[i] = np.sqrt(VxExp[i] ** 2 + VyExp[i] ** 2) * (VxExp[i] / np.abs(VxExp[i]))

for i in range(steps - 1):
    sExp[i + 1] = sExp[i] + VsExp[i + 1] * dt

e_array = np.zeros(steps)

for i in range(1, steps - 1):
    x = np.interp(sExp[i], sPath, xPath)
    p = 2 * A * x
    cos_beta = 1 / np.sqrt(1 + p * p)
    sin_beta = p / np.sqrt(1 + p * p)

    c = 2 * A / (1 + p * p) ** 1.5

    e_array[i] = (((2 * r ** 2) / (5 * h)) * asExp[i] + h * (g * sin_beta - asExp[i])) \
                 / (VsExp[i] * (c * VsExp[i] ** 2) + g * cos_beta)

    sExp[i + 1] = sExp[i] + VsExp[i + 1] * dt

    if e_array[i] == np.inf or e_array[i] == -np.inf or e_array[i] == np.nan:
        e_array[i] = 0


# print(e_array)
print(np.average(e_array))
