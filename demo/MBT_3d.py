# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 15:47:20 2021

@author: Jiayu Li@SUSTech
"""

# Model of the magnetic topological insulator MnBi2Te4
# 
# Model function
# MBT_3d_model(per, delta, fm)
# Inputs
#  1. per: (int) period in z direction, up to the parity
#          per = 1: a 1*1*1 unit cell
#          per = 2 or 0: a 1*1*2 supercell
#  2. delta:(float) strength of the exchange coupling
#  3. fm: (int) parameter for different type of the magnetic ordering
#          fm = 0: ferromagnetic ordering
#          fm = 1: antiferromagnetic ordering (available if per = 2)

from pythtb import * # import TB model class
from numpy import linalg as LA
import matplotlib.pyplot as plt
from LocalChernMarker import *
from MBT_model import *

per = 1
delta = 0.5
mbt_model = MBT_3d_model(per, delta, 0)

n_spin = mbt_model._nspin
n_orb = mbt_model._norb
n_bands = n_orb*n_spin

path = [[.5,.5,0.],[0.,0.,0.],[0.,0.,.5]]
label = (r'$M$', r'$\Gamma $', r'$Z$')
nk = 151

# Evaluating the bulk band with k_x, k_y, and k_z
(k_vec, k_dist, k_node) = mbt_model.k_path(path, nk, report = False)
evals = mbt_model.solve_all(k_vec)

# Figure for bulk bands
fig, ax = plt.subplots(figsize = [5,5])
ax.set_xlim(k_node[0], k_node[-1])
ax.set_ylim(-3, 3)
ax.set_xticks(k_node)
ax.set_xticklabels(label)
for n in range(len(k_node)):
  ax.axvline(x = k_node[n], linewidth = 0.5, color = 'k')
ax.set_title("MBT bulk band structure")
ax.set_ylabel("Energy")
for ib in range(n_bands):
    ax.plot(k_dist, evals[ib], color = 'blue')
fig.tight_layout()
fig.savefig("1_bulk_band.pdf")

print('Done.\n')