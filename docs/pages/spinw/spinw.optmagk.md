---
{title: spinw.optmagk method, link: spinw.optmagk, summary: determines the magnetic
    propagation vector, keywords: sample, sidebar: sw_sidebar, permalink: spinw_optmagk.html,
  folder: spinw, mathjax: 'true'}

---

### Syntax

`res = optmagk(obj,Name,Value)`

### Description

The function determines the optimal propagation vector by calculating the
Fourier transform of the Hamiltonian as a function of wave vector and
finding the wave vector that corresponds to the smalles eigenvalue of the
Hamiltonian. It also returns the normal vector that corresponds to the
rotating coordinate system. The optimization is achieved via
Particle-Swarm optimization.
 

### Input Arguments

`obj`
:pinw] object.

`kbase`
:    Provide a set of vectors that define the possible k-vector:
         k = sum_i C(i)*kbase(:,i);
     where the optimiser determines the C(i) values that correspond
     to the lowest ground state energy. kbase is a
     matrix with dimensions [3 nBase], where nBase <=3. The basis
     vectors have to be linearly independent.

### Name-Value Pair Arguments

`'Accepts'`
: all options of ndbase.pso.

### Output Arguments

res       Structure with the following fields:
              k       Value of the optimal k-vector, with values between 0
                      and 1/2.
              n       Normal vector, defines the rotation axis of the
                      rotating coordinate system.
              E       The most negative eigenvalue at the given propagation
                      vector.
              stat    Full output of the ndbase.pso() optimizer.

### See Also

[ndbase.pso]

{% include links.html %}
