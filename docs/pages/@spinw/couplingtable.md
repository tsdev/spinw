---
title: couplingtable( )
keywords: sample
summary: "creates tabulated list of all bonds as stored"
sidebar: product1_sidebar
permalink: couplingtable.html
folder: @spinw
mathjax: true
---
  creates tabulated list of all bonds as stored
 
  bonds = COUPLINGTABLE(obj,{bondIdx})
 
  Input:
 
  obj       spinw class object. 
  bondIdx   List of bond indices, by default all bonds will be output.
            Optional. If a bond index is mutiplied by -1, the table output
            is a matlab built in table type, works only for Matlab R2013b
            or later versions.
 
  Output:
 
  bonds is a struct type data that contains the following fields:
    table   Matrix, where every column defines a bond. The rows are the
            following: (dl_x, dl_y, dl_z, atom1, atom2, idx, mat_idx1,
            mat_idx2, mat_idx3). Where (dl_x, dl_y, dl_z) defines the
            translation vector between the origin of the unit cells of the
            two interacting atom (if they are in the same unit cell, all
            three components are zero) from atom1 to atom2. atom1 and atom2
            are the indices of the atoms in the obj.matom list. idx is the
            index of the bond, where equivalent bonds have identical
            indices, typically index is increasing with bond length. The
            last 3 rows (mat_idx) contains pointers to matrices if they
            are defined, otherwise zeros.
    bondv   Additional information for every bond defined in the .table
            field. The first three rows define the vector pointing from
            atom1 to atom2 in lattice units. The last row define the bond
            length in Angstrom.
    matrix  Contains the coupling matrix for every bond, dimensions are
            [3 3 nCoupling].
 
  Example:
 
  ...
  crystal.gencoupling
  bonds = crystal.couplingtable(-[1 2 3]).table
 
  This will list only the 1st, 2nd and 3rd neighbour bonds in a formatted
  table.
 
  See also SPINW.MATOM, SPINW.INTMATRIX, SPINW.ADDCOUPLING, SPINW.GENCOUPLING.
 
