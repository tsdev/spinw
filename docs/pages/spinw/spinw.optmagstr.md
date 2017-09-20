---
{title: spinw.optmagstr method, link: spinw.optmagstr, summary: optimises magnetic
    structure by minimizing the energy using non-linear optimization algorithms, keywords: sample,
  sidebar: sw_sidebar, permalink: spinw_optmagstr.html, folder: spinw, mathjax: 'true'}

---

### Syntax

`optm = optmagstr(obj, option1, value1, ...)`

### Description



### Examples

tri = sw_model('triAF',1);
X1 = [0 0 0 0 pi/2 0];
X2 = [0 1/2 1/2 0 pi/2 0];
optRes = tri.optmagstr('func',@gm_planar,'xmin',X1,'xmax',X2);
plot(tri)
The example determined the magnetic structure of the triangular lattice
antiferromagnet assuming planar magnetic structure and constraining the
moments into the [0 y z] plane (nTheta = 90 deg, nPhi = 0 deg or
n = [1 0 0]). Then plots the magnetic structure.

### Input Arguments

`obj`
: [spinw](spinw.html) object.

### Name-Value Pair Arguments

`func`
:unction that produce the magnetic moments, ordering wave
 ector and normal vector from the optimization
 arameters in the following form:
    [M, k, n] = @(x)func(M0, x)
 here M is (3,nMagExt) size matrix. k is the ordering
 ave vector, its size is (1,3). n is the normal vector
 f the spin rotation plane, its size is (1,3). The
 efault value is @gm_spherical3d. For planar magnetic structure
 se @gm_planar.

`xmin`
:inimum limit of the optimisation parameters, optional.

`xmax`
:aximum limit of the optimisation parameters, optional.

`x0`
:tarting value of the optimisation parameters. If empty
 r undefined, then random values are used.

`boundary`
:oundary conditions of the extended unit cell.
    'free'  Free, interactions between extedned unit cells are
            omitted.
    'per'   Periodic, interactions between extended unit cells
            are retained.
 efault is {'per' 'per' 'per'}.

`epsilon`
:he smalles value of incommensurability that is tolerated
 ithout warning. Default is 1e-5.

`nRun`
:umber of runs. If random starting parameters are given, the
 ptimisation process will be rerun nRun times and the best
 esult (lowest ground state energy per spin) will be saved in
 he result.

`title`
:ives a title string to the simulation that is saved in the
 utput.

`Limits`
: on selected prameters:

`Limits`
:be given on any input parameter of the constraint function by

`giving`
:name of the parameter, see the help of the used constraint

`function`
: the following format: optmagstr('ParName',[min max],...).

`For`
: to fix the nTheta value of @gm_planar during the optimisation

`to`
::

`optmagstr(obj,'func',@gm_planar,'nTheta',[0`
:bj,'func',@gm_planar,'nTheta',[0 0]);

`Optimisation`
:n parameters:

`tolx`
:   Minimum change of x when convergence reached, default
    value is 1e-4.

`tolfun`
:   Minimum change of the R value when convergence reached,
    default value is 1e-5.

`maxfunevals`
:   Maximum number of function evaluations, default value
    is 1e7.

`maxiter`
:   Maximum number of iterations, default value is 1e4.

### Output Arguments

'optm' is a struct type variable with the following fields:
obj       spinw object that contains the optimised magnetic structure.
x         Optimised paramters, dimensions are [1 nPar].
fname     Name of the contraint function.
xname     Cell containing the name of the x parameters, dimensions are
          [1 nPar].
e         Energy per spin in the optimised structure.
exitflag  Exit flag of the optimisation code, see fminsearch.
output    Detailed output of the optimisation code, see fminsearch.
param     Input parameters, stored in a struct.

### See Also

[spinw](spinw.html) \| [spinw.anneal](spinw_anneal.html) \| [gm_spherical3d](gm_spherical3d.html) \| [gm_planar](gm_planar.html) \| [fminsearch]

