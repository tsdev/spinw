---
{title: spinw.unit_cell property, link: spinw.unit_cell, summary: stores the atoms
    in the crystallographic unit cell, keywords: sample, sidebar: sw_sidebar, permalink: spinw_unit_cell.html,
  folder: spinw, mathjax: 'true'}

---
 
### Sub fields
 
`r`
: Positions of the atoms in the unit cell, stored in a
  matrix with dimensions of $$[3\times n_{atom}]$$, values are
  in lattice units.
 
`S`
: Spin quantum number of the atoms, stored in a row vector with
  $$n_{atom}$$ number of elements, non-magnetic atoms have `S=0`.
 
`label`
: Label of the atom, strings stored in a $$[1\times n_{atom}]$$
  cell.
 
`color`
: Color of the atom stored in a matrix with dimensions of $$[3\times n_{atom}]$$, where every
  column defines an RGB color with values between 0 and 255.
 
`ox`
: Oxidation number of the atom, stored in a $$[1\times n_{atom}]$$
  matrix.
 
`occ`
: Site occupancy in a $$[1\times n_{atom}]$$ matrix.
 
`b`
: Scattering length of the atoms for neutron and x-ray
  stored in a $$[2\times n_{atom}]$$ matrix, first row is neutron,
  second row is for x-ray.
 
`ff`
: Form factor of the site stored in a $$[2\times 9\times
  n_{atom}]$$ matrix, first row is the magnetic form factor for
  neutrons, the second row is the charge form factor for x-ray
  cross section.
 
`Z`
: Atomic number in a row vector.
 
`A`
: Atomic mass (N+Z) for isotopes and -1 for natural
  distribution of isotopes stored in a row vector.
 
`biso`
: Isotropic displacement factors in units of Å$$^2$$.
  Definition is the same as in
  [FullProf](https://www.ill.eu/sites/fullprof/), defining the
  Debye-Waller factor as $$W(d) = 1/8*b_{iso}/d^2$$ which is
  included in the structure factor as $$\exp(-2W(d))$$.
 

{% include links.html %}
