---
{title: sw_fstat, link: sw_fstat, summary: calculates thermodynamical averages, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_fstat, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`parOut = sw_fstat(state, parIn, T, E, M, nExt)`
  
### Description
  
`parOut = sw_fstat(state, parIn, T, E, M, nExt)` calculates statistical
properties of different physical variables over several sampled state.
The function is called by [spinw.anneal](spinw_anneal).
  
### Input Arguments
  
`state`
: Defines the task of the function.
  * `1`   Initialize the parOut structure.
  * `2`   Store the parameters of the physical state.
  * `3`   Calculate physical properties from the variable
          statistics.
  
`parIn`
: Same as `parOut`.
  
`T`
: Temperature of the system, row vector with $$n_T$$ number of elements.
  
`E`
: Energy of the system, row vector with $$n_T$$ number of elements.
  
`M`
: Magnetic moment of every atom in a matrix with dimensions of $$[d_{spin}\times n_{magExt}\cdot n_T]$$.
  
`nExt`
: Size of the magnetic supercell, column vector of 3 integers.
  
`kB`
: Boltzmann constant, units of temperature.
  
### Output Arguments
  
`parOut`
: Output parameter structure with the following fields:
  * `nStat`   The number of evaluated states.
  * `M`       $$\langle M\rangle$$ averaged over all magnetic moment stored
              in a matrix with dimensions of $$[d){spin}\times
              n_{magExt}\cdot n_T]$$.
  * `M2`      $$\langle M^2\rangle$$ averaged over all magnetic moment
              stored in a matrix with dimensions of $$[d){spin}\times
              n_{magExt}\cdot n_T]$$.
  * `E`       $$\langle E\rangle$$  summed over all magnetic moment.
  * `E2`      $$\langle E^2\rangle$$  summed over all magnetic moment.
 
For the final execution, the following parameters are calculated:
`parOut`
: Array of struct with $$n_T$$ number of elements:
  * `avgM`    Average components of the magnetisation over $$n_{stat}$$ runs,
              matrix with dimensions of $$[3\times n_{magExt}]$$.
  * `stdM`    Standard deviation of the mgnetisation components over
              $$n_{stat}$$ runs, matrix with dimensions of $$[3\times n_{magExt}]$$.
  * `avgE`    Average system energy per spin over $$n_{stat}$$ runs, scalar.
  * `stdE`    Standard deviation of the system energy per spin over
              $$n_{stat}$$ runs, scalar.
  * `T`       Final temperature of the sample.
  * `Cp`      Heat capacity of the sample: $$(\langle E^2\rangle-\langle E\rangle^2)/k_B/T^2$$.
  * `Chi`     Magnetic susceptibility of the sample: $$(\langle M^2\rangle-\langle M\rangle^2)/k_B/T$$.
  
### See Also
  
[spinw.anneal](spinw_anneal)
 

{% include links.html %}
