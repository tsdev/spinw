---
{title: spinw.unit property, link: spinw.unit, summary: stores the physical units
    for the Hamiltonian, keywords: sample, sidebar: sw_sidebar, permalink: spinw_unit,
  folder: spinw, mathjax: 'true'}

---
 
Default values are meV, T, Å and K for energy, magnetic
field, length and temperature, respectively.
 
### Sub fields
 
`kB`
: Boltzmann constant, default value is 0.0862 meV/K.
 
`muB`
: Bohr magneton, default values is 0.0579 meV/T.
 
`mu0`
: Vacuum permeability, default value is 201.335431 T$$^2$$Å$$^3$$/meV.
 
`label`
: Labels for distance, energy, magnetic field and temperature
stored in a cell with dimensions $$[1\times 4]$$.
 
`nformula`
: Number of formula units in the unit cell.
 
`qmat`
: Transformation matrix that converts the given $$Q$$ values into
the internal reciprocal lattice. The matrix has dimensions of
$$[3\times 3]$$.
 

{% include links.html %}
