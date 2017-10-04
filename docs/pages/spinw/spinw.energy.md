---
{title: spinw.energy method, link: spinw.energy, summary: calculates the ground state
    energy, keywords: sample, sidebar: sw_sidebar, permalink: spinw_energy, folder: spinw,
  mathjax: 'true'}

---
  
### Syntax
  
`E = energy(obj,Name,Value)`
  
### Description
  
`E = energy(obj,Name,Value)` calculates the classical ground state energy
per spin. The calculation correctly takes into account the magnetic
supercell. The function gives correct results on single-k magnetic
structures even defined on magnetic supercells. For multi-k magnetic
structures first a definition of a larger supercell is necessary where an
effective $$k=0$$ representation is possible.
  
### Examples
  
After optimising the magnetic structure (by minimizing the ground state
energy), the energy per spin is calculated. This can be compared to
different ground state structures to decide which is the right classical
ground state of the magnetic model in cryst. Here we use the triangular
lattice antiferromagnet where we define the magnetic structure on a
$$[3\times 3]$$ magnetic supercell where the optimal structure (120°
angle between neighboring spins) has a 0 propagation vector. In this case
the exact energy is $$3\cdot 1^2\cdot \cos(120^\circ) = -1.5$$.
 
```matlab
cryst = sw_model('triAF',1)
cryst.genmagstr('mode','random','nExt',[3 3 1])
cryst.optmagsteep('nRun',10)
cryst.energy
```
*Output*
```
Ground state energy: -1.500 meV/spin.
```
 
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
  
`'epsilon'`
: The smallest value of incommensurability that is tolerated 
  without warning. Default is $$10^{-5}$$.
  
### Output Arguments
  
`E`
: Energy per moment (anisotropy + exchange + Zeeman energy).
 
{% include warning.html content=" The calculated energy can be wrong for incommensurate
structures. For example a structure where the spins are rotating in $$XY$$
plane with an incommensurate wavevector of $$(1/3,0,0)$$. The function only
calculates the anisotropy energy in the first unit cell, that is for
single spin $$E_{aniso} = A_{xx}\cdot S_{xx}^2+A_{yy}\cdot S_{yy}^2$$.
While the anisotropy energy in reality is independent of the spin
orientation in the $$XY$$ plane $$E_{aniso}=3S\cdot (A_{xx}+A_{yy})/2$$. Thus
using `spinw.energy` on incommensurate structures together with single
ion anisotropy one has to be carefull! In the triangular case one has to
extend the unit cell to `nExt = [3 3 1]` (in the hexagonal setting), in
this case the energy will be correct." %}
  
### See Also
  
[spinw](spinw) \| [spinw.anneal](spinw_anneal) \| [spinw.newcell](spinw_newcell)
 

{% include links.html %}
