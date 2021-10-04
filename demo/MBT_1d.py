# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 16:45:51 2021

@author: Jiayu Li@SUSTech
"""

# Triangular stick model of the magnetic topological insulator MnBi2Te4
# with translational symmetry along z direction
# 
# Model function
# MBT_tri_1d_model(per, n_len, delta, fm)
# Inputs
#  1. per: (int) period in z direction, up to the parity
#          per = 1: a 1*1*1 unit cell
#          per = 2 or 0: a 1*1*2 supercell
#  2. n_len: (int) number of sites as the length of the equilateral triangle
#  3. delta:(float) strength of the exchange coupling
#  4. fm: (int) parameter for different type of the magnetic ordering
#          fm = 0: ferromagnetic ordering
#          fm = 1: antiferromagnetic ordering (available if per = 2)

from pythtb import * # import TB model class
from numpy import linalg as LA
import numpy as np
import matplotlib.pyplot as plt
from LocalChernMarker import *
from MBT_model import *

per = 1
n_len = 14
delta = 0.5
mbt_model = MBT_tri_1d_model(per, n_len, delta, 0)

n_spin = mbt_model._nspin
n_orb = mbt_model._norb
n_bands = n_orb*n_spin

path = [[-.5],[0.],[0.5]]
label = (r'$-Z$', r'$\Gamma $', r'$Z$')
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
fig.savefig("4_stick_band.pdf")

#==============================================================================

# solve orbital-resolved Chern number
n_occ = int(n_bands / 2) # half-filling as default
n_site = int((1+n_len)*n_len/2)
orbs = int(n_spin * n_orb / n_site)
cxz = lcm_1d(mbt_model, n_occ)

# collect the site-resolved Chern number
cxz1 = []
for i in range(n_site):
    cxz1.append(sum( cxz[i*orbs:(i+1)*orbs:1] ))
# collect site coordinates
rvec = np.zeros((2, n_site), dtype="float")
for i in range(n_site):
    for j in range(2):
        rvec[j,i] = ((mbt_model._orb[2*i])@(mbt_model._lat))[j]

# plot the distribution of Chern number
fig, ax = plt.subplots(figsize=[4, 4])
ax.scatter(rvec[0], rvec[1], c=cxz1, cmap="YlGnBu", s=200)
ax.set_xlabel("x")
ax.set_ylabel("y")
fig.tight_layout()
fig.savefig("5_lcm_distribution.pdf")

# collect the Chern markers along the relevant surface
sites = np.linspace(1, n_len, num = n_len)
surf = []
for i in range(n_len):
    surf.append(cxz1[np.argwhere(rvec[1] < 0.5)[i,0] ])
    
# plot the relevant Chern markers
fig, ax = plt.subplots(figsize = [4, 4])
ax.plot(sites, surf, color = 'blue', marker='o', markersize=12)
ax.axhline(0, c = 'k', ls = 'dashed')
ax.set_xlabel("site")
ax.set_ylabel(r'$C_{xz}(l) $')
fig.tight_layout()
fig.savefig("6_lcm_1d.pdf")

print('Done.\n')
