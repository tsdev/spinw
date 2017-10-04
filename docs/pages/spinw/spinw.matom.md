---
{title: spinw.matom method, link: spinw.matom, summary: generates magnetic lattice,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_matom, folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`mAtomList = matom(obj)`
  
### Description
  
`mAtomList = matom(obj)` is the same as [spinw.atom](spinw_atom), but only lists the
magnetic atoms, which have non-zero spin. Also this function stores the
generated list in [spinw.cache](spinw_cache).
  
### Output Arguments
  
`mAtomList`
: structure with the following fields:
  * `r`   Position of the magnetic atoms in a matrix with dimensions of 
    $$[3\times n_{magAtom}]$$.
  * `idx` Index in the symmetry inequivalent atom list [spinw.unit_cell](spinw_unit_cell) 
    stored in a row vector with $$n_{magAtom}]$$ number of elements.
  * `S`   Spin of the magnetic atoms stored in a row vectorwith
    $$n_{magAtom}]$$ number of elements.
  
### See Also
  
[spinw.atom](spinw_atom)
 

{% include links.html %}
