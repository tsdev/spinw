---
title: Magnetic structure
keywords: docs
sidebar: sw_sidebar
permalink: magstr.html
summary: Representation of magnetic structures and optimization
folder: documentation
mathjax: true
---

## Complex magnetic representation

If the expectation value of the spin operator is non-zero ($\langle \mathbf{S}\rangle\neq 0$) the spin vectors will make a pattern that can be mathematically best described by its Fourier transform. A concise description of this topic can be found in the article [A. Wills: Magnetic structures and their determination using group theory](https://doi.org/10.1051/jp4:2001906).

To define a magnetic structure first we need to find the correct magnetic unit cell. Although in magnetic structure determination most often the magnetic and the crystallographic unit cells are identical, in SpinW it is possible to define a magnetic supercell, that is the multiples of the crystallographic unit cell. This can help sometimes to optimize the spin wave calculation and gives more flexibility to test new models without changing the crystal structure. The magnetic structure is stored in the [spinw.mag_str] property of the SpinW class with the `nExt` subfield (`int32` type) storing the size of the magnetic supercell in units of the crystallographic unit cell. By default the magnetic supercell is identical to the crystallographic one which is given by `nExt=[1 1 1]`.

Once we defined the magnetic unit cell, the 




## Magnetic lattice

## Single-k representation
