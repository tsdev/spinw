---
{title: spinw.abc method, link: spinw.abc, summary: returns lattice parameters and
    angles, keywords: sample, sidebar: sw_sidebar, permalink: spinw_abc.html, folder: spinw,
  mathjax: 'true'}

---
  
### Syntax
  
`latvect = abc(obj)`
  
### Description
  
`latvect = abc(obj)` extracts the lattice vectors and angles from a
[spinw](spinw.html) object.
  
### Input Arguments
  
`obj`
: [spinw](spinw.html) object.
  
### Output Arguments
  
`latVect`
: Vetor with elements `[a, b, c, α, β, γ]`,
  contains the lattice parameters and angles by default in Å and
  degree units respectively (see [spinw.unit](spinw_unit.html) for details).
  
### See Also
  
[spinw.horace](spinw_horace.html)
 

{% include links.html %}
