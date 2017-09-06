---
{title: gm_spherical3d( ), keywords: sample, summary: magnetic structure constraint function with spherical parameterisation,
  sidebar: sw_sidebar, permalink: swfiles_gm_spherical3d.html, folder: swfiles, mathjax: 'true'}

---
  magnetic structure constraint function with spherical parameterisation
 
  [M, k, n, name, pname, limit] = GM_SPHERICAL3D(M0, x) 
 
  It generates the parameters of magnetic moments and normal vector from
  spherical (theta,phi) coordinates. All angles are in radian.
 
  Input:
 
  x         Input parameters in the following order:
            (Theta1, Phi1, Theta2, Phi2, ... , kx, ky, kz, nTheta, nPhi).
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
 
  See also GM_PLANAR.
 
