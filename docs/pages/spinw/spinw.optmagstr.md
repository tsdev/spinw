---
{title: spinw.optmagstr method, link: spinw.optmagstr, summary: general magnetic structure
    optimizer, keywords: sample, sidebar: sw_sidebar, permalink: spinw_optmagstr,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`optm = optmagstr(obj,Name,Value)`
  
### Description
  
`optm = optmagstr(obj,Name,Value)` is a general magnetic structure
optimizer that as the name suggests is general however usually less
efficient than [spinw.optmagk](spinw_optmagk) or [spinw.optmagsteep](spinw_optmagsteep). However this
function enables the usage of constraint functions to improve the
optimization. This function is most usefull if there is 1-2 parameters
that has to be optimized, such as a canting angle of the spins in
magnetic field. To optimize large number of spin angles
[spinw.optmagsteep](spinw_optmagsteep) might be faster.
  
### Examples
  
The example determines the propagation vector of the ground state of the
triangular lattice antiferromagnet. The magnetic structure is constrained
to be planar in the $$xy$$-plane. The [gm_planard](gm_planard) constraint function is
used where the first 3 parameter determined the propagation vector,
followed by the polar angles of the magnetic moments (here there is only
1 magnetic moment in the unit cell) which is fixed to 0. Finally the last
2 parameters corresponds to the polar angles of the normal to the
spin-plane which is the $$z$$-axis ($$\theta=0$$, $$\varphi=0$$). The optimized
magnetic structure is plotted.
 
```matlab
tri = sw_model('triAF',1)
X1 = [0 0 0 0 0 0]
X2 = [0 1/2 1/2 0 0 0]
optRes = tri.optmagstr('func',@gm_planard,'xmin',X1,'xmax',X2)
km = optRes.x(1:3)
```
*Output*
```
km =
         0    0.3333    0.3333
```
 
```matlab
plot(tri)
```
 
{% include image.html file="generated/spinw_optm_1.png" alt="swplot.zoom(1.5)" %}
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
  
`'func'`
: Function that produces the spin orientations, propagation vector and
  normal vector from the optimization parameters and has the following
  argument list:
  ```matlab
  [M, k, n] = @(x)func(M0, x)
  ```
 here `M` is matrix with dimensions of $$[3\times n_{magExt}]$$, `k` is the
 propagation vector (row vector with 3 elements), `n` is the normal vector
 of the spin rotation plane (row vector with 3 elements). The
 default value is `@gm_spherical3d`. For planar magnetic structures
 use `@gm_planar`.
  
`'xmin'`
: Lower limit of the optimisation parameters.
  
`'xmax'`
: Upper limit of the optimisation parameters.
  
`'x0'`
: Starting value of the optimisation parameters. If empty
  or undefined, then random values are used within the given limits.
  
`'boundary'`
: Boundary conditions of the magnetic cell:
  * `'free'`  Free, interactions between extedned unit cells are
            omitted.
  * `'per'`   Periodic, interactions between extended unit cells
            are retained.
 
  Default value is `{'per' 'per' 'per'}`.
  
`'epsilon'`
: The smallest value of incommensurability that is tolerated
  without warning. Default value is $$10^{-5}$$.
  
`'nRun'`
: Number of runs. If random starting parameters are given, the
  optimisation process will be rerun `nRun` times and the best
  result (lowest ground state energy per spin) will be kept.
  
`'title'`
: Gives a title string to the simulation that is saved in the
  output.
  
#### Limits on selected prameters
 
Limits can be given on any input parameter of the constraint function by
giving the name of the parameter. For parameter names see the help of the
used constraint function. Limits per optimization parameter can be given
in the following format: `optmagstr('ParName',[min max],...)`. For example
to fix the `nTheta` value of [gm_planar](gm_planar) during the optimisation to zero
use: `optmagstr(obj,'func',@gm_planar,'nTheta',[0 0])`.
 
  
#### Optimisation parameters
  
The optimization parameters are identical to the input options of the
Matlab built-in optimizer [fminsearch](https://www.mathworks.com/help/matlab/ref/fminsearch.html).
 
`'tolx'`
: Minimum change of `x` when convergence reached, default
    value is $$10^{-4}$$.
  
`'tolfun'`
: Minimum change of the $$R$$ value when convergence reached,
    default value is $$10^{-5}$$.
  
`'maxfunevals'`
: Maximum number of function evaluations, default value
    is $$10^7$$.
  
`'maxiter'`
: Maximum number of iterations, default value is $$10^4$$.
  
### Output Arguments
  
`optm`
: Struct type variable with the following fields:
  * `obj`       spinw object that contains the optimised magnetic structure.
  * `x`         Optimised paramters in a row vector with $$n_{par}$$ number
                of elements.
  * `fname`     Name of the contraint function.
  * `xname`     Cell containing the name of the $$x$$ parameters with
                  $$n_{par}$$ elements.
  * `e`         Energy per spin in the optimised structure.
  * `exitflag`  Exit flag of the optimisation code, see [fminsearch](https://www.mathworks.com/help/matlab/ref/fminsearch.html).
  * `output`    Detailed output of the optimisation code, see [fminsearch](https://www.mathworks.com/help/matlab/ref/fminsearch.html).
  * `param`     Input parameters, stored in a struct.
  
### See Also
  
[spinw](spinw) \| [spinw.anneal](spinw_anneal) \| [gm_spherical3d](gm_spherical3d) \| [gm_planar](gm_planar) \| [fminsearch]
 

{% include links.html %}
