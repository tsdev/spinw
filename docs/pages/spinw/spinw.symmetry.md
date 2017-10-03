---
{title: spinw.symmetry method, link: spinw.symmetry, summary: true if space group
    is used to generate couplings, keywords: sample, sidebar: sw_sidebar, permalink: spinw_symmetry,
  folder: spinw, mathjax: 'true'}

---

### Syntax

`sym = symmetry(obj)`

### Description

If true, equivalent couplings are generated based on the
crystal space group and all matrices (interaction, anisotropy
and g-tensor) are transformed according to the symmetry
operators. If false, equivalent couplings are generated based
on bond length, equivalent matrices won't be transformed
(all identical).
 
To change it use spinw.gencoupling with the forceNoSym option.
To remove all symmetry operators use spinw.nosym.
 

### See Also

[spinw](spinw) \| [spinw.nosym](spinw_nosym) \| [spinw.gencoupling](spinw_gencoupling)

{% include links.html %}
