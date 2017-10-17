---
title: Crystal structure
keywords: docs
sidebar: sw_sidebar
permalink: crystal.html
summary: Everything about defining crystal structure
folder: documentation
mathjax: 'true'
---

This section introduces how to define a crystal structure.

## Definition of the lattice

For defining a spin Hamiltonian, first we need a magnetic lattice. SpinW builds magnetic lattices on top of crystal structures, so first we need to define a crystal structure. For simple magnetic lattices this seems to be an unnecessary extra step, however to model real systems with space group symmetry, defining a magnetic lattice on top of a crystal structure make things easier.

Before we define a new crystal structure, we initialize an empty SpinW object:
```
>>model = spinw
>>model.lattice>>
```
The [spinw.lattice] property will store all the information that defines the lattice (without the atoms). The default values give a cubic lattice (angles are stored in radian) with a lattice parameter of 3 \\ang and no symmetry (`spinw.lattice.sym` field is empty).

The [spinw.genlattice] function can be used to modify the lattice:
```
>>>model = spinw
>>model.genlattice('lat_const',[2.71 2.71 13],'angled',[90 90 120],'spgr','R -3 m')
>>model.lattice>>
```
Here we modified the lattice parameters and angles (the `angled` parameter takes angles in units of degree) and defined a space group. Space groups with standard settings can be added simply by using its name, e.g. `'R -3 m'` which corresponds to the labels stored in the `symmetry.dat` file (to see the content of this file type `edit symmetry.dat` into the Matlab Command Window). The `spgr` option also accepts directly space group operators, see section [section.symmetry].


## Generation of crytal structure
## Coordinate systems
## Lattice transformations
## Atomic properties
