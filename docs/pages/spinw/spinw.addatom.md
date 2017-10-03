---
{title: spinw.addatom method, link: spinw.addatom, summary: adds new atom, keywords: sample,
  sidebar: sw_sidebar, permalink: spinw_addatom, folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`addatom(obj,Name,Value)`
  
### Description
  
`addatom(obj,Name,Value)` adds a new atom to the list of symmetry
inequivalent sites together with its properties, such as position, spin
quantum number, form factor, etc.
  
### Examples
  
To add a magnetic atom with $$S=1$$ at position $$r=(0,0,0)$$ and a
non-magnetic one at $$r=(1/2 0 0)$$ with red and blue color respectively
use the following command
 
```matlab
crystal = spinw;
crystal.genlattice('lat_const',[4 3 3])
crystal.addatom('r',[0 1/2; 0 0; 0 0],'S',[1 0],'color',{'red' 'blue'})
crystal.plot
```
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
  
`r`
: Atomic positions stored in a matrix with dimensions of $$[3\times
  n_{atom}]$$.
  
`label`
: Names of the atoms in a cell for plotting and form factor
  calculations (see `magion.dat`), e.g. `label={'atom1' 'atom2'
  'atom3'}`.
  Default value is `atomi`, where `i` is the atom index.
  
`S`
: Spin quantum number stored in a row vector with $$n_{atom}$$ elements,
  for non-magnetic atoms set S to zero. If not given the spin quantum
  number is guessed from the given label of the atom. For example if
  `label` is `MCr3+` or `Cr3+` then the $$S=3/2$$ high spin state is
  assumed for Cr$$^{3+}$$. The spin values for every ion is stored in the
  `magion.dat` file. If the atom type is unknown $$S=0$$ is assumed.
  
`color`
: RGB color of the atoms for plotting stored in a matrix with dimensions
  of $$[3\times n_{atom}]$$, where each column describes an RGB color. Each
  value is between 0 and 255. Default value is the color stored in the
  `atom.dat` file. Alternatively a name of the color can be given as a
  string, for example `'White'`, for multiple atoms package it into a
  cell. For the list of colors, see [swplot.color](swplot_color) or the `color.dat`
  file.
  
`ox`
: Oxidation number given as a double or it will be determined
  automatically from label. Default value is 0.
  
`occ`
: Occupancy, given as double. Default value is 1.
  
`formfact`
: Neutron scattering form factor, given as a row vector with 9 numbers,
  for details see [sw_mff](sw_mff). Also string labels can be used from the
  `magion.dat` file.
  
`formfactn`
: Same as the `formfact` option.
  
`formfactx`
: X-ray scattering form factor, given as 9 numbers, for details
  see [sw_cff](sw_cff), also labels can be used from the `xrayion.dat` file.
  
`Z`
: Atomic number, given as integer or determined from the atom label
  automatically. Default value is 113 (Unobtanium).
  
`A`
: Atomic mass, given as integer. Default is -1 for the natural
  mixture of isotopes.
  
`bn`
: Neutron scattering length, given as double. Not implemented yet.
  
`bx`
: X-ray scattering length.
  
`biso`
: Isotropic displacement factors in units of Å$$^2$$.
  Definition is the same as in
  [FullProf](https://www.ill.eu/sites/fullprof/), defining the
  Debye-Waller factor as $$W(d) = 1/8*b_{iso}/d^2$$, which is included in
  the structure factor as $$exp(-2W(d))$$.
  
`update`
: If true, existing atom with the same label and position as a
  new one will be updated. Default is true.
  
### Output Arguments
  
The function modifies the [spinw.unit_cell](spinw_unit_cell) property of the obj
[spinw](spinw) object.
  
### See Also
  
[spinw.genlattice](spinw_genlattice) \| [spinw.addmatrix](spinw_addmatrix) \| [swplot.color](swplot_color) \| [sw_mff](sw_mff) \| [sw_cff](sw_cff)
 

{% include links.html %}
