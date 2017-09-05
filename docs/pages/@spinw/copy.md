---
title: copy( )
keywords: sample
summary: "clones spinw object with all data"
sidebar: product1_sidebar
permalink: copy.html
folder: @spinw
mathjax: true
---
  clones spinw object with all data
 
  newObj = COPY(obj)
 
  Use this command instead of the '=' sign if you want to
  create an independent duplicate of an spinw class object.
 
  Input:
 
  obj       spinw class object.
 
  Output:
 
  newObj    New spinw class object that contains all the data of
            obj.
 
  Example:
 
  cryst = spinw;
  cryst.addmatrix('label','J1','value',3.1415)
 
  cryst1 = cryst;
  cryst2 = cryst.copy;
 
  cryst.addmatrix('label','J1','value',1)
  J1a = cryst1.matrix.mat;
  J1b = cryst2.matrix.mat;
 
  Where J1a will be a matrix with 1 in the diagonal, while J1b
  has 3.1415 in the diagonal. If cryst is changed, cryst1 will
  be changed as well and viece versa, since they point to the
  same object. However cryst2 is independent of cryst.
 
  See also SPINW, SPINW.STRUCT.
 
