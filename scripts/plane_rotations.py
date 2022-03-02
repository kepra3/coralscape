#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 24/2/22
@description: TODO
"""

import numpy as np
import matplotlib.pyplot as plt

# create the figure
fig = plt.figure()

# add axes
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel("Reef parallel")
ax.set_ylabel("Reef perpendicular")
ax.set_zlabel("Depth")

xx, yy = np.meshgrid(range(100), range(100))
a = -0.26
b = 0.90
c = -0.34
d = 10.55
# Plane equation
print(f"Plane equation: {a:.2f}x + {b:.2f}y + {c:.2f}z + {d:.2f} = 0")
#ax+by+cz=d
z = (d-(a*xx+b*yy))/c

# plot the plane
ax.plot_surface(xx, yy, z, alpha=0.5)
#angle with x-y plane n1=<0,0,1> n2=<a,b,c>
nxz = [0, 1, 0]  # angle with xz plane
nyz = [1, 0, 0]  # angle with yz plane
nxy = [0, 0, 1]  # angle with the xy plane
n2 = [a, b, c]
magxz = np.linalg.norm(nxz)
magyz = np.linalg.norm(nyz)
magxy = np.linalg.norm(nxy)
mag2 = np.linalg.norm(n2)


cosAx = np.dot(nxz, n2)/(magxz*mag2)
cosAy = np.dot(nyz, n2)/(magyz*mag2)
cosAz = np.dot(nxy, n2)/(magyz*mag2)
print('cosine of theta_xz is', cosAx)
print('cosine of theta_yz is', cosAy)
print('cosine of theta_xy is', cosAz)
Ax = np.arccos(cosAx) * 180./np.pi
print('xz angle is', Ax)
Ay = np.arccos(cosAy) * 180./np.pi
print('yz angle is', Ay)
Az = np.arccos(cosAz) * 180./np.pi
print('xy angle is', Az)
print('180 - xy angle is', 180 - Az)  # the angle probably is 70
plt.show()


