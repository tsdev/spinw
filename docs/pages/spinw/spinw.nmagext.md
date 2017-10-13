---
{title: spinw.nmagext method, link: spinw.nmagext, summary: number of magnetic sites,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_nmagext, folder: spinw,
  mathjax: 'true'}

---
 
### Syntax
 
`nMagExt = nmagext(obj)`
 
### Description
 
`nMagExt = nmagext(obj)` returns the number of magnetic sites
in the magnetic supercell. If the magnetic supercell (stored
in `spinw.mag_str.nExt` is identical to the crystal lattice)
the number of magnetic sites is equal to the number of
magnetic atoms in the unit cell. Where the number of magnetic
atoms in the unit cell can be calculated using [spinw.matom](spinw_matom).
 
### See Also
 
[spinw.matom](spinw_matom) \| [spinw.natom](spinw_natom)
 

{% include links.html %}
