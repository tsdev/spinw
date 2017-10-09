---
{title: spinw.moment method, link: spinw.moment, summary: calculates quantum correction
    on ordered moment, keywords: sample, sidebar: sw_sidebar, permalink: spinw_moment,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`M = moment(obj,Name,Value)`
  
### Description
  
`M = moment(obj,Name,Value)` calculates the spin expectation value
including the leading quantum and thermal fluctuations ($$S^{-1}$$ terms).
The magnon poulation is calculated at a given temperature $$T$$ integrated
over the Brillouin zone. To calculate the numerical integral the
Brillouin zone is sampled using Monte Carlo technique.
  
### Example
 
#### Triangular lattice antiferromagnet
 
The example calculates the spin expectation value at zero temperature on
the triangular lattice Heisenberg antiferromagnet. The result can be
compared with the following calculations: [A. V Chubukov, S. Sachdev,
and T. Senthil, J. Phys. Condens. Matter 6, 8891 (1994)](http://iopscience.iop.org/article/10.1088/0953-8984/6/42/019/meta): $$\langle S
\rangle = S - 0.261$$ and 
[S. J. Miyake, J. Phys. Soc. Japan 61, 983 (1992)](http://journals.jps.jp/doi/abs/10.1143/JPSJ.61.983): $$\langle S \rangle = S - 0.2613 +
0.0055/S$$ ($$1/S$$ is a higher order term neglected here).
 
```matlab
tri = sw_model('triAF',1)
M = tri.moment('nRand',1e7)
dS = 1-M.moment
```
*Output*
```
dS =
    0.2612
```
 
 
#### Square lattice antiferromagnet
 
The reduced moment of the Heisenberg square lattice antiferromagnet at
zero temperature can be compared to the published result of 
[D. A. Huse, Phys. Rev. B 37, 2380
(1988)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.37.2380)
$$\langle S \rangle = S - 0.197$$.
 
```matlab
sq = sw_model('squareAF',1)
M = sq.moment('nRand',1e7)
dS = 1-M.moment
```
*Output*
```
dS =
    0.1966    0.1966    0.1966    0.1966
```
 
 
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
  
`'nRand'`
: The number of random $$Q$$ points in the Brillouin-zone,
  default value is 1000.
  
`'T'`
: Temperature, default value is taken from `obj.single_ion.T` and the
  unit is stored in [spinw.unit](spinw_unit) with the default being K.
  
`'tol'`
: Tolerance of the incommensurability of the magnetic
  propagation wavevector. Deviations from integer values of the
  propagation vector smaller than the tolerance are
  considered to be commensurate. Default value is $$10^{-4}$$.
  
`'omega_tol'`
: Tolerance on the energy difference of degenerate modes when
  diagonalising the quadratic form, default value is $$10^{-5}$$.
  
### Output Arguments
  
`M`
: structure, with the following fields:
  * `moment`  Size of the reduced moments in a row vector with
    $$n_{magExt}$$ number of elements.
  * `T`       Temperature.
  * `nRand`   Number of random $$Q$$ points.
  * `obj`     The clone of the input `obj`.
 
### See Also
  
[spinw](spinw) \| [spinw.spinwave](spinw_spinwave) \| [spinw.genmagstr](spinw_genmagstr) \| [spinw.temperature](spinw_temperature)
 

{% include links.html %}
