---
{title: spinw.addcoupling method, link: spinw.addcoupling, summary: assigns an exchange
    matrix to a bond, keywords: sample, sidebar: sw_sidebar, permalink: spinw_addcoupling.html,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`addcoupling(obj,Name,Value)`
  
### Description
  
`addcoupling(obj,Name,Value)` assigns a matrix (will be used as exchange
matrix) to a given bond after bonds are generated using
[spinw.gencoupling](spinw_gencoupling.html).
  
### Examples
  
To add the $$J_1$$ diagonal matrix to all second neighbor bonds
between magnetic atoms use the following:
 
```matlab
cryst = sw_model('squareAF',1)
cryst.addmatrix('label','J1','value',diag([1 0.1 0.1]))
cryst.gencoupling
cryst.addcoupling('mat','J1','bond',2)
plot(cryst,'range',[2 2 1])
```
 
{% include image.html file="generated/spinw_addcoupling_1.png" alt="plot(cryst,'range',[2 2 1])" %}
  
### Input Arguments
  
`obj`
: [spinw](spinw.html) object.
  
### Name-Value Pair Arguments
  
`'mat'`
: Label (string) or index (integer) of the matrix that will be assigned to
  selected bonds, e.g. `'J1'`.
  
`'bond'`
: Integer that selects bonds, e.g. 1 for first neighbor, 2 for second
  neighbor, etc. The given value is compared to the `obj.coupling.idx`
  vector and the exchange matrix will be assigned to matching bonds.
  `'bond'` can be also a row vector to assign matrices to multiple bonds.
  
`'atom'`
: Contains labels of atoms (string) or index of atoms (integer) that is
  compared to [spinw.unit_cell](spinw_unit_cell.html) where all symmetry inequivalent atoms are
  stored. If a single string label or number is given, e.g. `'Cr1'` only
  Cr1-Cr1 bonds will be assigned. If a cell with 2 strings, e.g. `{'Cr1'
  'Cr2'}` only Cr1-Cr2 bonds will be assigned. Default value is `[]`.
  
`'subIdx'`
: If the above options are not enough to select the desired
  bonds, using `subIdx` bonds can be selected one-by-one from
  the list of bonds that fulfill the constraint of `atom` and `bond`.
  
`'type'`
: Type of the coupling with possible values of:
  * `'quadratic'`     Quadratic exchange, default.
  * `'biquadratic'`   Biquadratic exchange.
  
`'sym'`
: If `true`, symmetry operators will be applied on the exchange
  matrices to generate the coupling on symmetry equivalent
  bonds, if `false` all symmetry equivalent bonds will have the
  same exhcange matrix.
  
{% include warning.html content=" Setting `atom` or `subIdx` parameters will remove the symmetry
operations on the selected bonds. This means that assigning any
non-Heisenberg exchange matrix will break the space group defined in
`obj.lattice.sym`. Effectively reducing the symmetry of the given bond to
`P0`" %}
 
### Output Arguments
  
The function adds extra entries to the [spinw.coupling](spinw_coupling.html) property of
`obj`. Specifically it will modify `obj.coupling.mat_idx`,
`obj.coupling.type` and `obj.coupling.sym` matrices.
  
### See Also
  
[spinw](spinw.html) \| [spinw.gencoupling](spinw_gencoupling.html) \| [spinw.addmatrix](spinw_addmatrix.html)

{% include links.html %}
