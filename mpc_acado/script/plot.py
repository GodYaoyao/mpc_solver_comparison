#!/usr/bin/env python
import matplotlib.pyplot as plt
import os
from math import *

pwd = os.getcwd()

dt = 0.15
x = 0.
y = 0.
theta = -pi/6
v = 6.
w = 0.
refer_point = [(0, 1), (0.649519, 0.625), (1.24658, 0.282862), (1.88371, -0.102789), (2.4752, -0.445936), (3.03536, -0.800703), (3.6101, -1.15855), (4.16354, -1.49049), (4.69827, -1.79146), (5.18763, -2.07799), (5.62428, -2.30602), (6.04539, -2.50867)]
control_list = [(-0.535868, 0.200257), (-0.476821, 0.142657), (-0.426272, 0.097945), (-0.372709, 0.0642854), (-0.320365, 0.0397248), (-0.272996, 0.0222751), (-0.226305, 0.0106159), (-0.180565, 0.00371821), (-0.134631, 0.000212899), (-0.0861327, -0.00104017), (-0.0202673, -0.000300991), (-0.0202673, -0.000300991)]

calc_point = [(x, y, theta, v, w)]
for c in control_list:
    x += (cos(theta) * dt * v + cos(theta) * dt**2 / 2 * c[0])
    y += (sin(theta) * dt * v + sin(theta) * dt**2 / 2 * c[0])
    theta += (dt * w + dt**2 / 2 * c[1])
    v += dt * c[0]
    w += dt * c[1]
    calc_point.append((x, y, theta, v, w))

### vehicle origin point
plt.plot(0.0, 0.0, 'xb')
### use control calculate points
plt.scatter([c[0] for c in calc_point], [c[1] for c in calc_point],
            color='red', s=10, alpha=0.5, label="control_calculate")
### rviz show points
plt.scatter([r[0] for r in refer_point], [r[1] for r in refer_point],
            s=10, color='green', alpha=0.5, label="refer_point")
plt.legend()
plt.axis('equal')
plt.grid(True)
plt.show()
