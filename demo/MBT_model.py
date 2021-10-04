# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 15:47:06 2021

@author: Jiayu Li@SUSTech

Tight-binding model of MnBi2Te4
"""

from pythtb import * # import TB model class
import numpy as np

#==============================================================================

# MnBi2Te4 bulk model
# Inputs
#  1. per: (int) period in z direction, up to the parity
#          per = 1: a 1*1*1 unit cell
#          per = 2 or 0: a 1*1*2 supercell
#  2. delta:(float) strength of the exchange coupling
#  3. fm: (int) parameter for different type of the magnetic ordering
#          fm = 0: ferromagnetic ordering
#          fm = 1: antiferromagnetic ordering (available if per = 2)
def MBT_3d_model(per, delta, fm):
    #model parameters
    acon = 1.0
    ccon = 1.0
    m0 = -2
    m1 = 1
    m2 = 1
    c0 = 0.0
    c1 = 0.3
    c2 = 0.4
    v = 0.5
    vz = 1
    w = 1.5

    m_A = float(delta)
    m_B = float(delta)
    theta_A = 0.0 * np.pi
    phi_A = 0.0 * np.pi
    theta_B = float(fm) * np.pi
    phi_B = 0.0 * np.pi

    c_tilde = c0 + 2*c1 + 4*c2
    m_tilde = m0 + 2*m1 + 4*m2

    #lattice vectors and orbitals
    lat = [[acon, 0.0, 0.0], [acon*0.5, acon*0.5*np.sqrt(3), 0.0], [0.0, 0.0, ccon]]
    orb = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

    my_model = tb_model(3,3,lat,orb,nspin = 2)

    Id = [[1.0, 0.0], [0.0, 1.0]]
    sx = [[0.0, 1.0], [1.0, 0.0]]
    sy = [[0.0, -1j], [1j, 0.0]]
    sz = [[1.0, 0.0], [0.0, -1.0]]

    mm_A = m_A * (np.sin(theta_A)*np.cos(phi_A)*np.array(sx) + np.sin(theta_A)*np.sin(phi_A)*np.array(sy) + np.cos(theta_A)*np.array(sz))
    mm_B = m_B * (np.sin(theta_B)*np.cos(phi_B)*np.array(sx) + np.sin(theta_B)*np.sin(phi_B)*np.array(sy) + np.cos(theta_B)*np.array(sz))

    #onsite energy
    my_model.set_onsite(c_tilde + m_tilde, 0)
    my_model.set_onsite(c_tilde - m_tilde, 1)

    #hopping term
    #a1
    # 0->0
    my_model.set_hop(-2/3*(c2+m2)*np.array(Id), 0, 0, [1,0,0])
    # 1->1
    my_model.set_hop(-2/3*(c2-m2)*np.array(Id), 1, 1, [1,0,0])
    # 0->1
    my_model.set_hop(-v/3*1j*np.array(sx) - w/2*np.array(Id), 1, 0, [1,0,0])
    # 1->0
    my_model.set_hop(-v/3*1j*np.array(sx) + w/2*np.array(Id), 0, 1, [1,0,0])

    #a2
    # 0->0
    my_model.set_hop(-2/3*(c2+m2)*np.array(Id), 0, 0, [0,1,0])
    # 1->1
    my_model.set_hop(-2/3*(c2-m2)*np.array(Id), 1, 1, [0,1,0])
    # 0->1
    my_model.set_hop(-v/6*1j*np.array(sx) - v/2/np.sqrt(3)*1j*np.array(sy) + w/2*np.array(Id), 1, 0, [0,1,0])
    # 1->0
    my_model.set_hop(-v/6*1j*np.array(sx) - v/2/np.sqrt(3)*1j*np.array(sy) - w/2*np.array(Id), 0, 1, [0,1,0])

    #a1-a2
    # 0->0
    my_model.set_hop(-2/3*(c2+m2)*np.array(Id), 0, 0, [1,-1,0])
    # 1->1
    my_model.set_hop(-2/3*(c2-m2)*np.array(Id), 1, 1, [1,-1,0])
    # 0->1
    my_model.set_hop(-v/6*1j*np.array(sx) + v/2/np.sqrt(3)*1j*np.array(sy) + w/2*np.array(Id),1, 0, [1,-1,0])
    # 1->0
    my_model.set_hop(-v/6*1j*np.array(sx) + v/2/np.sqrt(3)*1j*np.array(sy) - w/2*np.array(Id), 0, 1, [1,-1,0])

    #a3
    # 0->0
    my_model.set_hop(-(c1+m1), 0, 0, [0,0,1])
    # 1->1
    my_model.set_hop(-(c1-m1), 1, 1, [0,0,1])
    # 0->1
    my_model.set_hop(-vz*0.5*1j*np.array(sz), 1, 0, [0,0,1])
    # 1->0
    my_model.set_hop(-vz*0.5*1j*np.array(sz), 0, 1, [0,0,1])

    if (per % 2) ==0:
        finite_model = my_model.make_supercell([[1,0,0],[0,1,0],[0,0,2]],to_home=True)
        finite_model.set_onsite(mm_A, 0, mode='add')
        finite_model.set_onsite(mm_A, 1, mode='add')
        finite_model.set_onsite(mm_B, 2, mode='add')
        finite_model.set_onsite(mm_B, 3, mode='add')
    else:
        finite_model = my_model
        finite_model.set_onsite(mm_A, 0, mode='add')
        finite_model.set_onsite(mm_A, 1, mode='add')

    return finite_model

#==============================================================================

# MnBieTe4 slab model
# Inputs
#  1. n_layers: (int) number of slabs in z direction
#  2. delta:(float) strength of the exchange coupling
#  3. fm: (int) parameter for different type of the magnetic ordering
#          fm = 0: ferromagnetic ordering
#          fm = 1: antiferromagnetic ordering (available if n_layers > 1)
def MBT_2d_model(n_layers, delta, fm):
    #model parameters
    acon = 1
    ccon = 1
    m0 = -2
    m1 = 1
    m2 = 1
    c0 = 0.0
    c1 = 0.3
    c2 = 0.4
    v = 0.5
    vz = 1
    w = 1.5

    m_A = float(delta)
    m_B = float(delta)
    theta_A = 0.0 * np.pi
    phi_A = 0.0 * np.pi
    theta_B = float(fm) * np.pi
    phi_B = 0.0 * np.pi

    c_tilde = c0 + 2*c1 + 4*c2
    m_tilde = m0 + 2*m1 + 4*m2

    #lattice vectors and orbitals
    lat = [[acon, 0.0, 0.0], [acon*0.5, acon*0.5*np.sqrt(3), 0.0], [0.0, 0.0, ccon]]
    orb = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

    my_model = tb_model(3,3,lat,orb,nspin = 2)

    Id = [[1.0, 0.0], [0.0, 1.0]]
    sx = [[0.0, 1.0], [1.0, 0.0]]
    sy = [[0.0, -1j], [1j, 0.0]]
    sz = [[1.0, 0.0], [0.0, -1.0]]

    mm_A = m_A * (np.sin(theta_A)*np.cos(phi_A)*np.array(sx) + np.sin(theta_A)*np.sin(phi_A)*np.array(sy) + np.cos(theta_A)*np.array(sz))
    mm_B = m_B * (np.sin(theta_B)*np.cos(phi_B)*np.array(sx) + np.sin(theta_B)*np.sin(phi_B)*np.array(sy) + np.cos(theta_B)*np.array(sz))

    #onsite energy
    my_model.set_onsite(c_tilde + m_tilde, 0)
    my_model.set_onsite(c_tilde - m_tilde, 1)

    #hopping term
    #a1
    # 0->0
    my_model.set_hop(-2/3*(c2+m2)*np.array(Id), 0, 0, [1,0,0])
    # 1->1
    my_model.set_hop(-2/3*(c2-m2)*np.array(Id), 1, 1, [1,0,0])
    # 0->1
    my_model.set_hop(-v/3*1j*np.array(sx) - w/2*np.array(Id), 1, 0, [1,0,0])
    # 1->0
    my_model.set_hop(-v/3*1j*np.array(sx) + w/2*np.array(Id), 0, 1, [1,0,0])

    #a2
    # 0->0
    my_model.set_hop(-2/3*(c2+m2)*np.array(Id), 0, 0, [0,1,0])
    # 1->1
    my_model.set_hop(-2/3*(c2-m2)*np.array(Id), 1, 1, [0,1,0])
    # 0->1
    my_model.set_hop(-v/6*1j*np.array(sx) - v/2/np.sqrt(3)*1j*np.array(sy) + w/2*np.array(Id), 1, 0, [0,1,0])
    # 1->0
    my_model.set_hop(-v/6*1j*np.array(sx) - v/2/np.sqrt(3)*1j*np.array(sy) - w/2*np.array(Id), 0, 1, [0,1,0])

    #a1-a2
    # 0->0
    my_model.set_hop(-2/3*(c2+m2)*np.array(Id), 0, 0, [1,-1,0])
    # 1->1
    my_model.set_hop(-2/3*(c2-m2)*np.array(Id), 1, 1, [1,-1,0])
    # 0->1
    my_model.set_hop(-v/6*1j*np.array(sx) + v/2/np.sqrt(3)*1j*np.array(sy) + w/2*np.array(Id),1, 0, [1,-1,0])
    # 1->0
    my_model.set_hop(-v/6*1j*np.array(sx) + v/2/np.sqrt(3)*1j*np.array(sy) - w/2*np.array(Id), 0, 1, [1,-1,0])

    #a3
    # 0->0
    my_model.set_hop(-(c1+m1), 0, 0, [0,0,1])
    # 1->1
    my_model.set_hop(-(c1-m1), 1, 1, [0,0,1])
    # 0->1
    my_model.set_hop(-vz*0.5*1j*np.array(sz), 1, 0, [0,0,1])
    # 1->0
    my_model.set_hop(-vz*0.5*1j*np.array(sz), 0, 1, [0,0,1])

    finite_model = my_model.cut_piece(n_layers, 2)

    for ni in range(n_layers):
        if (ni % 2) == 0:
            finite_model.set_onsite(mm_A, 2*ni, mode='add')
            finite_model.set_onsite(mm_A, 2*ni+1, mode='add')
        else:
            finite_model.set_onsite(mm_B, 2*ni, mode='add')
            finite_model.set_onsite(mm_B, 2*ni+1, mode='add')

    return finite_model

#==============================================================================

# Triangular MnBi2Te4 model with PBC in z direction
# Inputs
#  1. per: (int) period in z direction, up to the parity
#          per = 1: a 1*1*1 unit cell
#          per = 2 or 0: a 1*1*2 supercell
#  2. n_len: (int) number of sites as the length of the equilateral triangle
#  3. delta:(float) strength of the exchange coupling
#  4. fm: (int) parameter for different type of the magnetic ordering
#          fm = 0: ferromagnetic ordering
#          fm = 1: antiferromagnetic ordering (available if per = 2)

def MBT_tri_1d_model(per, n_len, delta, fm):
    #model parameters
    acon = 1.0
    ccon = 1.0
    m0 = -2
    m1 = 1
    m2 = 1
    c0 = 0.0
    c1 = 0.3
    c2 = 0.4
    v = 0.5
    vz = 1
    w = 1.5

    m_A = float(delta)
    m_B = float(delta)
    theta_A = 0.0 * np.pi
    phi_A = 0.0 * np.pi
    theta_B = float(fm) * np.pi
    phi_B = 0.0 * np.pi

    c_tilde = c0 + 2*c1 + 4*c2
    m_tilde = m0 + 2*m1 + 4*m2

    #lattice vectors and orbitals
    lat = [[acon, 0.0, 0.0], [acon*0.5, acon*0.5*np.sqrt(3), 0.0], [0.0, 0.0, ccon]]
    orb = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

    my_model = tb_model(3,3,lat,orb,nspin=2)

    Id = [[1.0, 0.0], [0.0, 1.0]]
    sx = [[0.0, 1.0], [1.0, 0.0]]
    sy = [[0.0, -1j], [1j, 0.0]]
    sz = [[1.0, 0.0], [0.0, -1.0]]

    mm_A = m_A * (np.sin(theta_A)*np.cos(phi_A)*np.array(sx) + np.sin(theta_A)*np.sin(phi_A)*np.array(sy) + np.cos(theta_A)*np.array(sz))
    mm_B = m_B * (np.sin(theta_B)*np.cos(phi_B)*np.array(sx) + np.sin(theta_B)*np.sin(phi_B)*np.array(sy) + np.cos(theta_B)*np.array(sz))

    #onsite energy
    my_model.set_onsite(c_tilde + m_tilde, 0)
    my_model.set_onsite(c_tilde - m_tilde, 1)

    #hopping term
    #a1
    # 0->0
    my_model.set_hop(-2/3*(c2+m2)*np.array(Id), 0, 0, [1,0,0])
    # 1->1
    my_model.set_hop(-2/3*(c2-m2)*np.array(Id), 1, 1, [1,0,0])
    # 0->1
    my_model.set_hop(-v/3*1j*np.array(sx) - w/2*np.array(Id), 1, 0, [1,0,0])
    # 1->0
    my_model.set_hop(-v/3*1j*np.array(sx) + w/2*np.array(Id), 0, 1, [1,0,0])

    #a2
    # 0->0
    my_model.set_hop(-2/3*(c2+m2)*np.array(Id), 0, 0, [0,1,0])
    # 1->1
    my_model.set_hop(-2/3*(c2-m2)*np.array(Id), 1, 1, [0,1,0])
    # 0->1
    my_model.set_hop(-v/6*1j*np.array(sx) - v/2/np.sqrt(3)*1j*np.array(sy) + w/2*np.array(Id), 1, 0, [0,1,0])
    # 1->0
    my_model.set_hop(-v/6*1j*np.array(sx) - v/2/np.sqrt(3)*1j*np.array(sy) - w/2*np.array(Id), 0, 1, [0,1,0])

    #a1-a2
    # 0->0
    my_model.set_hop(-2/3*(c2+m2)*np.array(Id), 0, 0, [1,-1,0])
    # 1->1
    my_model.set_hop(-2/3*(c2-m2)*np.array(Id), 1, 1, [1,-1,0])
    # 0->1
    my_model.set_hop(-v/6*1j*np.array(sx) + v/2/np.sqrt(3)*1j*np.array(sy) + w/2*np.array(Id),1, 0, [1,-1,0])
    # 1->0
    my_model.set_hop(-v/6*1j*np.array(sx) + v/2/np.sqrt(3)*1j*np.array(sy) - w/2*np.array(Id), 0, 1, [1,-1,0])

    #a3
    # 0->0
    my_model.set_hop(-(c1+m1), 0, 0, [0,0,1])
    # 1->1
    my_model.set_hop(-(c1-m1), 1, 1, [0,0,1])
    # 0->1
    my_model.set_hop(-vz*0.5*1j*np.array(sz), 1, 0, [0,0,1])
    # 1->0
    my_model.set_hop(-vz*0.5*1j*np.array(sz), 0, 1, [0,0,1])

    if (per % 2)==0:
        m0_model = my_model.make_supercell([[1,0,0],[0,1,0],[0,0,2]],to_home=True)
        m0_model.set_onsite(mm_A, 0, mode='add')
        m0_model.set_onsite(mm_A, 1, mode='add')
        m0_model.set_onsite(mm_B, 2, mode='add')
        m0_model.set_onsite(mm_B, 3, mode='add')
    else:
        m0_model = my_model
        m0_model.set_onsite(mm_A, 0, mode='add')
        m0_model.set_onsite(mm_A, 1, mode='add')

    m1_model = m0_model.cut_piece(n_len, 1)
    m2_model = m1_model.cut_piece(n_len, 0)

    remove = []
    for ni in range(n_len-1):
        for nj in range(n_len*2*(per)-2*(ni+1)*(per)):
            remove.append(n_len**2*2*(per)-1-2*n_len*ni*(per)-nj)
    finite_model = m2_model.remove_orb(remove)

    return finite_model

#==============================================================================

# Triangular MnBi2Te4 model with OBC in z direction
# Inputs
#  1. n_layers: (int) number of layers in z direction
#  2. n_len: (int) number of sites as the length of the equilateral triangle
#  3. delta:(float) strength of the exchange coupling
#  4. fm: (int) parameter for different type of the magnetic ordering
#          fm = 0: ferromagnetic ordering
#          fm = 1: antiferromagnetic ordering (available if n_layers > 1)

def MBT_tri_0d_model(n_layers, n_len, delta, fm):
    #model parameters
    acon = 1.0
    ccon = 1.0
    m0 = -2
    m1 = 1
    m2 = 1
    c0 = 0.0
    c1 = 0.3
    c2 = 0.4
    v = 0.5
    vz = 1
    w = 1.5

    m_A = float(delta)
    m_B = float(delta)
    theta_A = 0.0 * np.pi
    phi_A = 0.0 * np.pi
    theta_B = float(fm) * np.pi
    phi_B = 0.0 * np.pi

    c_tilde = c0 + 2*c1 + 4*c2
    m_tilde = m0 + 2*m1 + 4*m2

    #lattice vectors and orbitals
    lat = [[acon, 0.0, 0.0], [acon*0.5, acon*0.5*np.sqrt(3), 0.0], [0.0, 0.0, ccon]]
    orb = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

    my_model = tb_model(3,3,lat,orb,nspin=2)

    Id = [[1.0, 0.0], [0.0, 1.0]]
    sx = [[0.0, 1.0], [1.0, 0.0]]
    sy = [[0.0, -1j], [1j, 0.0]]
    sz = [[1.0, 0.0], [0.0, -1.0]]

    mm_A = m_A * (np.sin(theta_A)*np.cos(phi_A)*np.array(sx) + np.sin(theta_A)*np.sin(phi_A)*np.array(sy) + np.cos(theta_A)*np.array(sz))
    mm_B = m_B * (np.sin(theta_B)*np.cos(phi_B)*np.array(sx) + np.sin(theta_B)*np.sin(phi_B)*np.array(sy) + np.cos(theta_B)*np.array(sz))

    #onsite energy
    my_model.set_onsite(c_tilde + m_tilde, 0)
    my_model.set_onsite(c_tilde - m_tilde, 1)

    #hopping term
    #a1
    # 0->0
    my_model.set_hop(-2/3*(c2+m2)*np.array(Id), 0, 0, [1,0,0])
    # 1->1
    my_model.set_hop(-2/3*(c2-m2)*np.array(Id), 1, 1, [1,0,0])
    # 0->1
    my_model.set_hop(-v/3*1j*np.array(sx) - w/2*np.array(Id), 1, 0, [1,0,0])
    # 1->0
    my_model.set_hop(-v/3*1j*np.array(sx) + w/2*np.array(Id), 0, 1, [1,0,0])

    #a2
    # 0->0
    my_model.set_hop(-2/3*(c2+m2)*np.array(Id), 0, 0, [0,1,0])
    # 1->1
    my_model.set_hop(-2/3*(c2-m2)*np.array(Id), 1, 1, [0,1,0])
    # 0->1
    my_model.set_hop(-v/6*1j*np.array(sx) - v/2/np.sqrt(3)*1j*np.array(sy) + w/2*np.array(Id), 1, 0, [0,1,0])
    # 1->0
    my_model.set_hop(-v/6*1j*np.array(sx) - v/2/np.sqrt(3)*1j*np.array(sy) - w/2*np.array(Id), 0, 1, [0,1,0])

    #a1-a2
    # 0->0
    my_model.set_hop(-2/3*(c2+m2)*np.array(Id), 0, 0, [1,-1,0])
    # 1->1
    my_model.set_hop(-2/3*(c2-m2)*np.array(Id), 1, 1, [1,-1,0])
    # 0->1
    my_model.set_hop(-v/6*1j*np.array(sx) + v/2/np.sqrt(3)*1j*np.array(sy) + w/2*np.array(Id),1, 0, [1,-1,0])
    # 1->0
    my_model.set_hop(-v/6*1j*np.array(sx) + v/2/np.sqrt(3)*1j*np.array(sy) - w/2*np.array(Id), 0, 1, [1,-1,0])

    #a3
    # 0->0
    my_model.set_hop(-(c1+m1), 0, 0, [0,0,1])
    # 1->1
    my_model.set_hop(-(c1-m1), 1, 1, [0,0,1])
    # 0->1
    my_model.set_hop(-vz*0.5*1j*np.array(sz), 1, 0, [0,0,1])
    # 1->0
    my_model.set_hop(-vz*0.5*1j*np.array(sz), 0, 1, [0,0,1])

    m0_model = my_model.cut_piece(n_layers, 2)

    for ni in range(n_layers):
        if (ni % 2) == 0:
            m0_model.set_onsite(mm_A, 2*ni, mode='add')
            m0_model.set_onsite(mm_A, 2*ni+1, mode='add')
        else:
            m0_model.set_onsite(mm_B, 2*ni, mode='add')
            m0_model.set_onsite(mm_B, 2*ni+1, mode='add')

    m1_model = m0_model.cut_piece(n_len, 1)
    m2_model = m1_model.cut_piece(n_len, 0)

    remove = []
    for ni in range(n_len-1):
        for nj in range(n_len*2*n_layers-2*(ni+1)*n_layers):
            remove.append(n_len**2*2*n_layers-1-2*n_len*ni*n_layers-nj)
    finite_model = m2_model.remove_orb(remove)

    return finite_model
