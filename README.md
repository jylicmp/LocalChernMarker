# LocalChernMarker

**Author:** Jiayu Li

A modified script for computing the **local Chern marker** in 0D, 1D, and 2D systems. It is based on the [partial Chern number package (pcn)](https://github.com/nicodemosvarnava/pcn) of [PythTB](http://www.physics.rutgers.edu/pythtb/index.html), rewritten with a matrix algorithm that achieves roughly **10× speedup** over the original implementation.

> **Compatibility:** This code is only compatible with **PythTB version 1.x**. It is not tested with PythTB 2.x.

---

## Table of Contents

- [LocalChernMarker.py](#localchernmarkerpy)
- [Demo](#demo)
- [Updates](#updates)
- [References](#references)

---

## LocalChernMarker.py

The script provides three main functions for the orbital-resolved Chern marker in tight-binding models with 2, 1, or 0 crystal momenta:

### `lcm_2d(my_model, n_occ)`

For a **slab model** with in-plane (xy) translational symmetry and finite stacking along *z*. Returns **C<sub>xy</sub>(l)** with *l* the orbital index (layer, spin, etc.).

### `lcm_1d(my_model, n_occ)`

For a **stick model** with open boundaries in *x*, *y* and translational symmetry along *z* (k<sub>z</sub>). By default returns **C<sub>xz</sub>(l)**. Useful to assess whether the surface magnetic gap dominates the trivial finite-size gap, e.g. for half-quantized surface response in axion insulators. See Ref. [1].

### `lcm_0d(my_model, n_occ, i, j)`

For a tight-binding model with **open boundaries in all three directions**. Returns **C<sub>ij</sub>(l)** with *i*, *j* ∈ {0, 1, 2} for (x, y, z).

---

## Demo

Examples compute the local Chern marker in effective **MnBi<sub>2</sub>Te<sub>4</sub>** models built with [PythTB](http://www.physics.rutgers.edu/pythtb/index.html).

**`MBT_model.py`** defines four effective MnBi<sub>2</sub>Te<sub>4</sub> models:

| Model | Function | Description |
|-------|----------|-------------|
| Bulk | `MBT_3d_model(per, delta, fm)` | 3D reciprocal space (k<sub>x</sub>, k<sub>y</sub>, k<sub>z</sub>) |
| Slab | `MBT_2d_model(n_layers, delta, fm)` | Stacking along *z*, 2D reciprocal space (k<sub>x</sub>, k<sub>y</sub>) |
| Stick | `MBT_tri_1d_model(per, n_len, delta, fm)` | Translational symmetry along *z* (k<sub>z</sub>), in-plane triangular geometry |
| Finite | `MBT_tri_0d_model(n_layers, n_len, delta, fm)` | Triangular prism, no translational symmetry |

Model details: Ref. [1] and [2].

---

## Updates

| Date | Changes |
|------|---------|
| 2023-10-30 | Support for `nspin=1` systems (suggestion by Rafael Gonzalez-Hernandez); eigenvector handling via NumPy; improved readability; new demo: Haldane model as `nspin=1` case |
| 2021-09-30 | Initial release |

---

## References

[1] Y. Wan, J. Li, and Q. Liu, *Topological magnetoelectric response in ferromagnetic axion insulators*, [Natl. Sci. Rev. **9**, nwac138 (2022)](https://doi.org/10.1093/nsr/nwac138).

[2] R.-X. Zhang, F. Wu, and S. Das Sarma, *Möbius Insulator and Higher-Order Topology in MnBi<sub>2n</sub>Te<sub>3n+1</sub>*, [Phys. Rev. Lett. **124**, 136407 (2020)](https://doi.org/10.1103/PhysRevLett.124.136407).

[3] R. Bianco and R. Resta, *Mapping topological order in coordinate space*, [Phys. Rev. B **84**, 241106(R) (2011)](https://doi.org/10.1103/PhysRevB.84.241106).

[4] T. Rauch, T. Olsen, D. Vanderbilt, and I. Souza, *Geometric and nongeometric contributions to the surface anomalous Hall conductivity*, [Phys. Rev. B **98**, 115108 (2018)](https://doi.org/10.1103/PhysRevB.98.115108).

[5] N. Varnava and D. Vanderbilt, *Surfaces of axion insulators*, [Phys. Rev. B **98**, 245117 (2018)](https://doi.org/10.1103/PhysRevB.98.245117).
