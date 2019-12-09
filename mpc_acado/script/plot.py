#!/usr/bin/env python
import matplotlib.pyplot as plt
import os
from math import *

pwd = os.getcwd()

dt = 0.1
x = 0.
y = 0.
theta = -pi/6
v = 1.
w = 0.
refer_point = [(0, 1), (0.433013, 0.75), (0.864957, 0.519839), (1.32997, 0.288714), (1.77648, 0.0892567), (2.24396, -0.128873), (2.71131, -0.359726), (3.18967, -0.597368), (3.66598, -0.842069), (4.1907, -1.08231), (4.73851, -1.32419), (5.28693, -1.59021)]
calculate_list = [(0, 0,-0.523599), (0.0952628, -0.055,-0.515363), (0.208378, -0.119071,-0.493385), (0.340488, -0.190112,-0.462524), (0.492626, -0.265967,-0.426467), (0.665608, -0.344562,-0.387803), (0.860014, -0.423975,-0.348168), (1.07621, -0.502445,-0.308444), (1.31442, -0.57834,-0.268972), (1.57471, -0.65009,-0.229773), (1.85338, -0.715271,-0.190738), (2.14266, -0.771128,-0.151759)]
control_list = [(2, 1.6471), (2, 1.10148), (2, 0.675049), (2, 0.364146), (2, 0.157267), (2, 0.0367471), (2, -0.0187255), (2, -0.031685), (2, -0.0230452), (1.23807, -0.00960564), (0.449044, -0.00174126)]

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
### refer points
plt.scatter([r[0] for r in refer_point], [r[1] for r in refer_point],
            s=10, color='green', alpha=0.5, label="refer_point")
### calculate points
plt.scatter([c[0] for c in calculate_list], [c[1] for c in calculate_list],
            s=10, color='blue', alpha=0.5, label="calculate_point")
plt.legend()
plt.axis('equal')
plt.grid(True)
plt.show()
