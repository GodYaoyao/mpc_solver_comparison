#!/usr/bin/env python
import matplotlib.pyplot as plt
import os
from math import *

pwd = os.getcwd()

dt = 0.1
x = 0
y = 0
phi = -0.547976
v = 1
w = 0
refer_point = [(0, 1), (0.433013, 0.75), (0.856179, 0.480982), (1.29914, 0.223054), (1.77071, -0.0495814), (2.27531, -0.311191), (2.76344, -0.574002), (3.19823, -0.834335), (3.61854, -1.11515), (4.05003, -1.41654), (4.49156, -1.70436), (4.89689, -1.97048)]
solver_point = [(0, 0), (-0.600764, -0.502203), (-0.713283, -0.491577), (-0.827468, -0.480851), (-0.941845, -0.470091), (-1.05519, -0.459322), (-0.553027, 1), (-0.551386, 1.2), (-0.546857, 1.4), (-0.540166, 1.6), (-0.531921, 1.8), (-0.52261, 2)]
control_list = [(0, 0), (0.00553763, -2), (0.00145392, -2), (0.000169353, -2), (3.32012e-321, -2), (3.32012e-321, -2), (-2, -2), (-2, -2), (-2, -2), (-2, -2), (-2, -2)]

calculate_point = [(x, y, phi, v, w)]
for c in control_list:
    x += (cos(phi) * dt * v + cos(phi) * dt ** 2 / 2 * c[0])
    y += (sin(phi) * dt * v + sin(phi) * dt ** 2 / 2 * c[0])
    phi += (dt * w + dt ** 2 / 2 * c[1])
    v += dt * c[0]
    w += dt * c[1]
    calculate_point.append((x, y, phi, v, w))

### vehicle origin point
plt.plot(0.0, 0.0, 'xb')
### refer_point
plt.scatter([r[0] for r in refer_point], [r[1] for r in refer_point],
            s=10, color='green', alpha=0.5, label="refer_point")
### solver_point
plt.scatter([c[0] for c in solver_point], [c[1] for c in solver_point],
            s=10, color='blue', alpha=0.5, label="solver_point")
### use control calculate points
plt.scatter([c[0] for c in calculate_point], [c[1] for c in calculate_point],
            color='red', s=10, alpha=0.5, label="control_calculate")
plt.legend()
plt.axis('equal')
plt.grid(True)
plt.show()
