---
{title: gm_planar( ), link: gm_planar, summary: planar magnetic structure constraint
    function, keywords: sample, sidebar: sw_sidebar, permalink: gm_planar.html, folder: swfiles,
  mathjax: 'true'}

---
 
:code:`[S, k, n, name, pname, limit] = GM_PLANAR(M0, x)`
 
The function generates the parameters of arbitrary planar magnetic
structure from phi angles (radian), ordering wave vector (rlu) and spin
plane normal vector (xyz).
 
Parameters
----------
 
x:
      Input parameters in the following order: 
      :math:`(\varphi_1, \varphi_2, ... , k_x, k_y, k_z, n_\theta, n_\phi)`.
absS:
      Size of the spins: :math:`(S_1, S_2, ...)` or scalar if all
      moments are equal.
 
Returns
-------
 
S:
      Array, containing the spin orientations with dimensions of [3 nMagExt].
      Every column contain the :math:`[S_x; S_y; S_z]` magnetic moment components of
      a magnetic atom in the xyz coordinate system.
k:
      Magnetic ordering wavevector in rlu units in a row vector.
n:
      Normal vector around which the spins are rotating for non-zero
      k-vector in a row vector.
name:
      String, storing the name of the function. Optional.
pname:
      Name of the input parameters in a cell: {'Phi1_rad', ...}.
      Optional.
limit:
      Limits on the input parameters, dimensions are [2 nParam]. Every
      column contains a lower and upper limit on the corresponding
      parameter. Optional.
 
See also
--------
gm_spherical3d, gm_planard.
 

