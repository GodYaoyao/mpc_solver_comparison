#!/usr/bin/env python
import matplotlib.pyplot as plt
import os
from math import *

pwd = os.getcwd()

dt = 0.1
x = 0
y = 0
phi = 0
v = 10.1336
w = 0.100727
refer_point = [(0.0253109, -0.0334626), (1.03203, -0.090364), (2.03398, -0.17271), (3.02916, -0.278666), (4.01601, -0.408028), (4.99294, -0.560535), (5.95851, -0.735242), (6.91133, -0.930938), (7.84971, -1.14805), (8.77164, -1.38787), (9.67524, -1.6507), (10.5591, -1.93513), (11.4217, -2.24026), (12.2601, -2.56825), (13.0716, -2.92059)]
solver_point = [(0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0)]
control_list = [(3.55271e-14, -1.00727), (-2, 0), (-2, 0), (-2, 0), (-2, 0), (-2, 0), (-2, 0), (-2, 0), (-2, 0), (-2, 0), (-2, 0), (-2, 0), (-2, 0), (-2, 0)]

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
