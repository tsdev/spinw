---
{title: spinw.unit( ), summary: Stores the physical units in the Hamiltonian., keywords: sample,
  sidebar: sw_sidebar, permalink: spinw_unit.html, folder: spinw, mathjax: 'true'}

---
Stores the physical units in the Hamiltonian.
Defaults are meV, Tesla Angstrom and Kelvin.
Sub fields are:
  kB      Boltzmann constant, default is 0.0862 [meV/K]
  muB     Bohr magneton, default is 0.0579 [meV/T]
  mu0     vacuum permeability, 201.335431 [T^2*Angstrom^3/meV]
  qmat    transformation matrix that converts the given q-values
