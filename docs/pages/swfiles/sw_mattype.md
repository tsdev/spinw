---
{title: sw_mattype, link: sw_mattype, summary: classifies square matrices, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_mattype, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`type = sw_mattype(mat)`
  
`type = sw_mattype(mat,epsilon)`
 
### Description
  
`type = sw_mattype(mat)` determines the type of the input matrix `mat`
which stacked $$[3\times 3]$$ matrices. It determines the type of exchnge
interaction that the matrix belongs to.
  
{% include note.html content=" Also works on symbolic matrices, but keep all symbols real for consistent
result!" %}
 
### Input Arguments
  
`mat`
: Matrix with dimensions of $$[3\times 3\times N]$$.
  
`epsilon`
: optional error bar on small matrix elements, default value is $$10^{-5}$$.
  
### Output Arguments
  
`type`
: Row vector with $$N$$ elements each having one of the following value:
  * `1`   Heisenberg exchange,
  * `2`   anisotropic exchange,
  * `3`   DM interaction,
  * `4`   general matrix.
 
*[DM]: Dzyaloshinskii-Moriya
 

{% include links.html %}
