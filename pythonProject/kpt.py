#!/usr/bin/env python

# Toy graphene model

# Copyright under GNU General Public License 2010, 2012, 2016
# by Sinisa Coh and David Vanderbilt (see gpl-pythtb.txt)

from __future__ import print_function
from pythtb import *
# import TB model class
import numpy as np
import matplotlib.pyplot as plt

#定义坐标轴
fig = plt.figure()
ax1 = plt.axes(projection='3d')
#ax = fig.add_subplot(111,projection='3d')  #这种方法也可以画多个子图

# define lattice vectors
lat=[[1.0,0.0],[0.5,np.sqrt(3.0)/2.0]]
# define coordinates of orbitals
orb=[[1./3.,1./3.],[2./3.,2./3.]]

# make two dimensional tight-binding graphene model
my_model=tb_model(2,2,lat,orb)

# set model parameters
delta=0.0
t=-1.0

# set on-site energies
my_model.set_onsite([-delta,delta])
# set hoppings (one for each connected pair of orbitals)
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_hop(t, 0, 1, [ 0, 0])
my_model.set_hop(t, 1, 0, [ 1, 0])
my_model.set_hop(t, 1, 0, [ 0, 1])

# print tight-binding model
my_model.display()

i = int(3)
j = int(3)
K = 0.7
M = 0.5
w_square=wf_array(my_model,[i,j])
all_kpt=np.zeros((71,51,1))
for i in np.arange(0, K, 0.01):
    for j in np.arange(0, M, 0.01):
        kpt =np.array([i, j])
        all_kpt[i,j,:]=kpt
        (eval, evec) = my_model.solve_one(kpt, eig_vectors=True)
        w_square[i,j]=evec
        plt.plot(i,j,eval,"gray")
        plt.show()

print("eval:",eval)
print(all_kpt)
