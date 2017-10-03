---
{title: spinw.single_ion property, link: spinw.single_ion, summary: stores single
    ion terms of the Hamiltonian, keywords: sample, sidebar: sw_sidebar, permalink: spinw_single_ion,
  folder: spinw, mathjax: 'true'}

---
 
### Sub fields
 
`aniso`
: Row vector that contains $$n_{magatom}$$ integers, each integer
  assignes one of the matrices from the [spinw.matrix](spinw_matrix) property
  to a magnetic atom in the generated [spinw.matom](spinw_matom) list as a single
  ion anisotropy. Zero value of `aniso` means no single ion
  anisotropy for the corresponding magnetic atom.
 
`g`
: Row vector with $$n_{magatom}$$ integers, each integer
  assignes one of the matrices from the [spinw.matrix](spinw_matrix) property
  to a magnetic atom in the spinw.matom list as a
  g-tensor. Zero value of `g` means a default g-value of 2 for
  the corresponding atoms.
 
`field`
: External magnetic field stored in a row vector with 3 elements,
  unit is defined in [spinw.unit](spinw_unit) (default unit is Tesla).
 
`T`
: Temperature, scalar, unit is defined in [spinw.unit](spinw_unit) (default
  unit is Kelvin).
 

{% include links.html %}
