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
solver_point = [(0, 0), (1.01052, 0), (2.01549, 7.30617e-05), (3.01501, -0.0197746), (4.0081, -0.0762777), (4.9926, -0.179327), (5.9658, -0.332516), (6.92522, -0.535344), (7.86893, -0.785049), (8.79573, -1.07798), (9.70491, -1.41054), (10.5961, -1.77972), (11.4691, -2.18328), (12.3238, -2.61978), (13.1597, -3.08837)]
control_list = [(-0.568756, -2), (-0.539555, -2), (-0.512365, -1.41057), (-0.49052, -0.684996), (-0.474739, -0.189841), (-0.462847, 0.112949), (-0.450788, 0.264871), (-0.433927, 0.307555), (-0.408082, 0.279724), (-0.370171, 0.214897), (-0.318547, 0.139951), (-0.253154, 0.074215), (-0.17563, 0.0287451), (-0.0894324, 0.00595962)]

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
