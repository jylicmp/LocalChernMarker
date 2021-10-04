# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 16:08:23 2021

@author: Jiayu Li@SUSTech
"""

# Slab model of the magnetic topological insulator MnBi2Te4
# 
# Model function
# MBT_2d_model(n_layers, delta, fm)
# Inputs
#  1. n_layers: (int) number of slabs in z direction
#  2. delta:(float) strength of the exchange coupling
#  3. fm: (int) parameter for different type of the magnetic ordering
#          fm = 0: ferromagnetic ordering
#          fm = 1: antiferromagnetic ordering (available if n_layers > 1)

from pythtb import * # import TB model class
from numpy import linalg as LA
import numpy as np
import matplotlib.pyplot as plt
from LocalChernMarker import *
from MBT_model import *

n_layers = 10
delta = 0.5
mbt_model = MBT_2d_model(n_layers, delta, 0)

n_spin = mbt_model._nspin
n_orb = mbt_model._norb
n_bands = n_orb*n_spin


path = [[.5,.5],[0.,0.],[1./3.,2./3.]]
label = (r'$M$', r'$\Gamma $', r'$K$')
nk = 151

#==============================================================================

# solve bands
(k_vec, k_dist, k_node) = mbt_model.k_path(path, nk, report = False)
evals = mbt_model.solve_all(k_vec)

# plot slab bands
fig, ax = plt.subplots(figsize = [5, 5])
ax.set_xlim(k_node[0] ,k_node[-1])
ax.set_ylim(-3,3)
ax.set_xticks(k_node)
ax.set_xticklabels(label)
for n in range(len(k_node)):
  ax.axvline(x = k_node[n], linewidth = 0.5, color = 'k')
ax.set_title("MBT slab band structure")
ax.set_ylabel("Energy")
for ib in range(n_bands):
    ax.plot(k_dist, evals[ib], color = 'blue')
fig.tight_layout()
fig.savefig("2_slab_band.pdf")

#==============================================================================

# solve orbital-resolved Chern number
n_occ = int(n_bands / 2) # half-filling as default
orbs = int(n_spin * n_orb / n_layers)
cxy = lcm_2d(mbt_model, n_occ)

# collect the layer-resolved Chern number
layers = np.linspace(1, n_layers, num = n_layers)
cxy1 = []
for i in range(n_layers):
    cxy1.append(sum( cxy[i*orbs:(i+1)*orbs:1] ))
    
# plot the layer-resolved Chern number
fig, ax = plt.subplots(figsize = [4, 4])
ax.plot(layers, cxy1, color = 'blue', marker='o', markersize=12)
ax.axhline(0, c = 'k', ls = 'dashed')
ax.set_xlabel("layer")
ax.set_ylabel(r'$C_{xy}(l) $')
fig.tight_layout()
fig.savefig("3_lcm_2d.pdf")

print('Done.\n')