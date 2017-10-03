---
{title: spinw.unitcell method, link: spinw.unitcell, summary: returns information
    on atoms in the crystallographic unit cell, keywords: sample, sidebar: sw_sidebar,
  permalink: spinw_unitcell, folder: spinw, mathjax: 'true'}

---

### Syntax

`unit_cell_info = unitcell(obj, idx)`

### Description

The function returns information on symmetry inequivalent atoms. 
 

### Examples

...
cryst.unit_cell = unitcell(cryst,[1 3]);
The example keeps only the first and third symmetry inequivalent atoms in
cryst object.
...
cryst.unit_cell = unitcell(cryst,'O');
The example keeps only the Oxygen atoms in cryst object.

### Input Arguments

`obj`
:pinw] object.

`idx`
:    Selects certain atoms. If undefined UNIT_CELL(obj) or
     obj.UNIT_CELL returns information on all atoms. The selection
     can be also done according to the atom labels, in this case
     either a string of the label or cell of strings for several
     labels can be given.

### Output Arguments

'unit_cell_info' is a tructure with that contains all the fields of
unit_cell.

### See Also

[spinw.addtwin](spinw_addtwin) \| [spinw.twinq](spinw_twinq) \| [spinw.unit_cell](spinw_unit_cell)

{% include links.html %}
