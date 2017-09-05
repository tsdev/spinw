---
title: gm_planard( )
keywords: sample
summary: "planar magnetic structure constraint function"
sidebar: product1_sidebar
permalink: gm_planard.html
folder: swfiles
mathjax: true
---
  planar magnetic structure constraint function 
 
  [M, k, n, name, pname, limit] = GM_PLANAR(M0, x) 
 
  It generates the parameters of arbitrary planar magnetic structure from
  phi angles, ordering wave vector and spin plane normal vector. All angles
  are in degree.
 
  Input:
 
  x         Input parameters in the following order:
            (Phi1, Phi2, ... , kx, ky, kz, nTheta, nPhi).
  M0        Size of magnetic moments: (M1, M2, ...) or scalar if all
            moments are equal.
 
  Output:
 
  M         Array, containing the magnetic moments, dimensions are
            [3 nMagExt]. Every column contain the [Mx; My; Mz] magnetic
            moment components of a magnetic atom in the xyz coordinate
            system.
  k         Magnetic ordering wavevector in r.l.u., dimensions are [1 3].
  n         Normal vector to the plane of the incommensurate spins (if k
            non-zero).
 
  Optional outputs:
  only produced if the output is requested.
 
  name      Name of the function.
  pname     Name of the input parameters in a cell: {'Param1' 'Param2',...}
  limit     Default limits on the input parameters, dimensions are [2 nX].
            Every column contains a lower and upper limit on the parameter.
 
  See also GM_SPHERICAL3D.
 
