---
{title: spinw.optmagk method, link: spinw.optmagk, summary: determines the magnetic
    propagation vector, keywords: sample, sidebar: sw_sidebar, permalink: spinw_optmagk,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`res = optmagk(obj,Name,Value)`
  
### Description
  
`res = optmagk(obj,Name,Value)` determines the optimal propagation vector
using the Luttinger-Tisza method. It calculates the Fourier transform of
the Hamiltonian as a function of wave vector and finds the wave vector
that corresponds to the smalles global eigenvalue of the Hamiltonian. It
also returns the normal vector that corresponds to the rotating
coordinate system. The global optimization is achieved using
Particle-Swarm optimizer.
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
 
`kbase`
: Provides a set of vectors that span the space for possible propagation
  vectors:
 
  $$ \mathbf{k} = \sum_i C(i)\cdot \mathbf{k}_{base}(i);$$
 
  where the optimiser determines the $$C(i)$$ values that correspond
     to the lowest ground state energy. $$\mathbf{k}_{base}$$ is a
     matrix with dimensions $$[3\times n_{base}]$$, where $$n_{base}\leq 3$$. The basis
     vectors have to be linearly independent.
  
The function also accepts all options of [ndbase.pso].
  
### Output Arguments
  
`res`
: Structure with the following fields:
  * `k`       Value of the optimal k-vector, with values between 0
                      and 1/2.
  * `n`       Normal vector, defines the rotation axis of the
                      rotating coordinate system.
  * `E`       The most negative eigenvalue at the given propagation
                      vector.
  * `stat`    Full output of the [ndbase.pso] optimizer.
  
### See Also
  
[ndbase.pso]
 

{% include links.html %}
