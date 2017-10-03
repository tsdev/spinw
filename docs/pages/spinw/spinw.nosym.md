---
{title: spinw.nosym method, link: spinw.nosym, summary: removes the space group symmetry,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_nosym, folder: spinw, mathjax: 'true'}

---

### Syntax

`nosym(obj)`

### Description

The function reduces the crystal symmetry to  0 and but keeps all atomic
positions, that are now all symmetry inequivalent.
 

### Examples

sw_addsym('x+1/2,y+1/2,z;x+1/2,y,z+1/2;x,y+1/2,z+1/2','FCC');
cryst = spinw;
cryst.genlattice('lat_const',[8 8 8],'sym','FCC')
cryst.addatom('r',[0 0 0],'label','Atom1')
cryst.nosym
The example creates four equivalent atomic positions, after the nosym()
command, the cryst.unit_cell.r contains the four generated positions,
that are not symmetry equivalent any more.

### Input Arguments

`obj`
: [spinw](spinw) object.

### Output Arguments

The obj input will have obj.lattice.sym field equal to zero and the
obj.unit_cell field will contain all the generated atomic positions.

### See Also

[spinw](spinw) \| [spinw.newcell](spinw_newcell)

{% include links.html %}
