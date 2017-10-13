---
{title: sw_fsub, link: sw_fsub, summary: simple graph vertex coloring, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_fsub, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`cgraph = sw_fsub(conn, next)`
  
### Description
  
`cgraph = sw_fsub(conn, next)` creates a simple graph vertex coloring,
determines non-connected sublattices for Monte-Carlo calculation.
  
### Input Arguments
  
`conn`
: Contains edge indices which are connected
  `conn(1,idx)-->conn(2,idx)` stored in a matrix with dimensions of $$[2times n_{conn}]$$.
  
`nExt`
: Size of the magnetic supercell in a row vector with 3 integers.
  
### Output Arguments
  
`cGraph`
: Vector, that assigns every magnetic moment to a sublattice.
  
### See Also
  
[spinw.anneal](spinw_anneal)
 

{% include links.html %}
