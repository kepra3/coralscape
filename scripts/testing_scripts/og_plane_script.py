#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 24/2/22
@description: plane script from dad
"""
import numpy as np
import matplotlib.pyplot as plt

# create the figure
fig = plt.figure()

# add axes
ax = fig.add_subplot(111,projection='3d')

xx, yy = np.meshgrid(range(100), range(100))
a=-0.25
b= 0.91
c=-0.32
d=10.36
#ax+by+cz=d
z = (d-(a*xx+b*yy) )/c

# plot the plane
ax.plot_surface(xx, yy, z, alpha=0.5)
#angle with x-y plane n1=<0,0,1> n2=<a,b,c>
n1=[0,1,0]
n2=[a,b,c]
mag1=np.linalg.norm(n1)
mag2=np.linalg.norm(n2)


cosA=np.dot(n1,n2)/(mag1*mag2)
print(cosA)
A=np.arccos(cosA)
print(A*180./np.pi)
plt.show()