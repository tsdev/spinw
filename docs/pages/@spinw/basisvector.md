---
title: basisvector( )
keywords: sample
summary: "generates basis vectors of the crystal lattice"
sidebar: product1_sidebar
permalink: basisvector.html
folder: @spinw
mathjax: true
---
  generates basis vectors of the crystal lattice
 
  basisVec = BASISVECTOR(obj, {norm})
 
  Input:
 
  obj       spinw class object.
  norm      If true, the basis vectors are normalized to 1, otherwise the
            length is equal to the lattice constants. Default is false.
            Optional.
 
  Output:
 
  basisVec  Stores the three basis vectors in columns, dimensions are 
            [3 3].
 
  The 3x3 basisVec matrix can be used also as a coordinate transformation
  matrix from the relative atomic position to positions in the xyz
  coordinate system in Angstrom units.
 
  Example:
 
  To change coordinate system:
 
  relative atomic positions --> xyz
    r_xyz = basisvector * [ra; rb; rc];
 
  reciprocal lattice units --> Angstrom^-1 (xyz coordinate system)
    Q_xyz =  [h k l] * 2*pi*inv(basisvector);
 
  See also SPINW, SPINW.ABC, SPINW.RL.
 
