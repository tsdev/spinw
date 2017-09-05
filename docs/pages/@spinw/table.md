---
title: table( )
keywords: sample
summary: "outputs easy to read tables of internal data"
sidebar: product1_sidebar
permalink: table.html
folder: @spinw
mathjax: true
---
  outputs easy to read tables of internal data
 
  T = SPINW.TABLE(obj,type,{index},{showVal})
 
  The function returns a table in Matlab R2013b or newer while in older
  versions a struct.
 
  For the matrix labels in the list of bonds, the '>>' sign means that the
  matrix value is determined using the symmetry operations.
 
 
  Input:
 
  obj       SpinW object.
  type      String, determines the type of data to show, values:
                'matom'     properties of magnetic atoms in the unit cell
                'matrix'    list of matrices
                'ion'       single ion term in the Hamiltonian
                'bond'      properties of selected bonds
                'mag'       magnetic structure
  index     Indexing into the type of data to show, depending on the option
            type:
                'bond'      indexes the bonds (1 for first neighbors,
                            etc.), if empty all bonds will be shown.
                'mag'       Indexes the propagation vectors, the
                            magnetization of the selected propagation
                            vector will be shown.
            Default value is 1, if empty vector ([]) is given, all
            bonds/propagation vector will be shown.
  showVal   Also show the values of the single ion terms and exchange
            values. The values shown  are the true exchange values after
            the symmetry operations (if there is any). Default is false.
 
  Output:
 
  T         Matlab table object.
 
