# LocalChernMarker
This is a modified script about the calculation of local Chern marker in a matrix algorithm for 0,1,2 dimensional system. The routine bases on the [partial Chern number package](https://github.com/nicodemosvarnava/pcn) of PythTB ([http://www.physics.rutgers.edu/pythtb/index.html](http://www.physics.rutgers.edu/pythtb/index.html)).

**LocalChernMarker.py**

This script contains three functions: *lcm_2d(my_mode,n_occ)*, *lcm_1d(my_mode,n_occ)*, and *lcm_0d(my_mode,n_occ,i,j)* for calculating the orbital-resolved Chern marker in the spinful tight-binding model with 2,1,0 crystal momenta, respectively.

(1) *lcm_2d(my_model,n_occ)* focuses on a slab model with translational symmetry along in-plane directions (xy), vertically stacking finite slabs along z direction. This function returns Cxy(l) with l the orbital index, including layer, spin, etc.

(2) *lcm_1d(my_model,n_occ)* is prepared for a stick model with open boundaries in x,y directions but translationally invariant along z direction (kz). In default, this function returns Cxz(l).

It is forceful to estimate whether the surface magnetic gap dominates the trivial finite size gap, especially in studying the half-quantized surface response in axion insulators. For more information see

[1] Y. Wan, J. Li, and Q. Liu, Topological Magnetoelectric Response in Ferromagnetic Axion Insulators, [arXiv:2109.14160 (2021)](https://arxiv.org/abs/2109.14160).

(3) *lcm_0d(my_model,n_occ,i,j)* is prepared for a tight-binding model with open boundary condition in all three directions. This function gives Cij(l) with i,j = (0,1,2) for (x,y,z).

**Demo**

Here we provide examples using the script to calculate the local Chern marker in effective MnBi2Te4 models constructed by [PythTB](http://www.physics.rutgers.edu/pythtb/index.html).

*MBT_model.py* contains four effective MnBi2Te4 models, including bulk model *MBT_3d_model(per, delta, fm)* with 3D reciprocal space (kx,ky,kz), slab model *MBT_2d_model(n_layers, delta, fm)* stacking along z direction with 2D reciprocal space (kx,ky), stick model *MBT_tri_1d_model(per, n_len, delta, fm)* with translational symmetry in z direction (kz) and regular triangle geometry in-plane, and finite site model *MBT_tri_0d_model(n_layers, n_len, delta, fm)* of triangular prism geometry without any translational symmetry. Details on the MnBi2Te4 model see Ref. [1] and

[2] R.-X. Zhang, F. Wu, and S. Das Sarma, Möbius Insulator and Higher-Order Topology in MnBi2nTe3n+1, [Phys. Rev. Lett. 124, 136407 (2020)](https://doi.org/10.1103/PhysRevLett.124.136407).

**Remark**

This script is a modified version of [partial Chern number package](https://github.com/nicodemosvarnava/pcn). The calculation is revised in a matrix algorithm that speeds up nearly 10 times than the primary code.

The formula used are followed

[3] Raffaello Bianco and Raffaele Resta, Mapping topological order in coordinate space, [Phys. Rev. B 84, 241106(R) (2011)](https://doi.org/10.1103/PhysRevB.84.241106).

[4] T. Rauch, T. Olsen, D. Vanderbilt, and I. Souza, Geometric and nongeometric contributions to the surface anomalous Hall conductivity, [Phys. Rev. B 98, 115108 (2018)](https://doi.org/10.1103/PhysRevB.98.115108).

[5] N. Varnava and D. Vanderbilt, Surfaces of axion insulators, [Phys. Rev. B 98, 245117 (2018)](https://doi.org/10.1103/PhysRevB.98.245117).
