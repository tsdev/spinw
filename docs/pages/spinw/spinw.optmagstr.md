---
{title: spinw.optmagstr( ), summary: optimises magnetic structure by minimizing the
    energy using non-linear optimization algorithms, keywords: sample, sidebar: sw_sidebar,
  permalink: spinw_optmagstr.html, folder: spinw, mathjax: 'true'}

---
 
optm = OPTMAGSTR(obj, Option1, Value1, ...)
 
Input:
 
obj       spinw class object.
 
Options:
 
func      Function that produce the magnetic moments, ordering wave
          vector and normal vector from the optimization
          parameters in the following form:
              [M, k, n] = @(x)func(M0, x)
          where M is (3,nMagExt) size matrix. k is the ordering
          wave vector, its size is (1,3). n is the normal vector
          of the spin rotation plane, its size is (1,3). The
          default is @gm_spherical3d. For planar magnetic structure
          use @gm_planar.
xmin      Minimum limit of the optimisation parameters, optional.
xmax      Maximum limit of the optimisation parameters, optional.
x0        Starting value of the optimisation parameters. If empty
          or undefined, then random values are used.
boundary  Boundary conditions of the extended unit cell.
              'free'  Free, interactions between extedned unit cells are
                      omitted.
              'per'   Periodic, interactions between extended unit cells
                      are retained.
          Default is {'per' 'per' 'per'}.
epsilon   The smalles value of incommensurability that is tolerated
          without warning. Default is 1e-5.
nRun      Number of runs. If random starting parameters are given, the
          optimisation process will be rerun nRun times and the best
          result (lowest ground state energy per spin) will be saved in
          the result.
title     Gives a title string to the simulation that is saved in the
          output.
 
Limits only on selected prameters:
 
Limits can be given on any input parameter of the constraint function by
giving the name of the parameter, see the help of the used constraint
function in the following format: optmagstr('ParName',[min max],...).
For example to fix the nTheta value of @gm_planar during the optimisation
to zero use:
optmagstr(obj,'func',@gm_planar,'nTheta',[0 0]);
 
Optimisation parameters:
 
tolx          Minimum change of x when convergence reached, default
              value is 1e-4.
tolfun        Minimum change of the R value when convergence reached,
              default value is 1e-5.
maxfunevals   Maximum number of function evaluations, default value
              is 1e7.
maxiter       Maximum number of iterations, default value is 1e4.
 
 
Output:
 
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
 
Example:
 
tri = sw_model('triAF',1);
X1 = [0 0 0 0 pi/2 0];
X2 = [0 1/2 1/2 0 pi/2 0];
optRes = tri.optmagstr('func',@gm_planar,'xmin',X1,'xmax',X2);
plot(tri)
 
The example determined the magnetic structure of the triangular lattice
antiferromagnet assuming planar magnetic structure and constraining the
moments into the [0 y z] plane (nTheta = 90 deg, nPhi = 0 deg or
n = [1 0 0]). Then plots the magnetic structure.
 
See also SPINW, SPINW.ANNEAL, GM_SPHERICAL3D, GM_PLANAR, FMINSEARCH.
 

