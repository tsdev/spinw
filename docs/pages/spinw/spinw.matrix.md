---
{title: spinw.matrix property, link: spinw.matrix, summary: stores 3x3 matrices for
    using them in the Hailtonian, keywords: sample, sidebar: sw_sidebar, permalink: spinw_matrix,
  folder: spinw, mathjax: 'true'}

---
 
### Sub fields
 
`mat`
: Stores the actual values of 3x3 matrices, in a matrix with
dimensions of $$[3\times 3\times n_{matrix}]$$, if assigned for a 
bond, the unit of energy is stored in [spinw.unit](spinw_unit) (default value 
is meV).
 
`color`
: Color assigned for every matrix, stored in a
  matrix with dimensions of $$[3\times n_{matrix}]$$, with each
  column defining an RGB value.
 
`label`
: Label for every matrix, stored as string in a cell with
  dimensions of $$[1\times n_{matrix}]$$.
 

{% include links.html %}
