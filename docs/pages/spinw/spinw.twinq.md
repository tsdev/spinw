---
{title: spinw.twinq method, link: spinw.twinq, summary: calculates equivalent Q point
    in twins, keywords: sample, sidebar: sw_sidebar, permalink: spinw_twinq, folder: spinw,
  mathjax: 'true'}

---
  
### Syntax
  
`[qTwin, rotQ] = twinq(obj, {Q0})`
  
### Description
  
`[qTwin, rotQ] = twinq(obj, {q0})` calculates the $$Q$$ values in the twin
coordinate systems, in rlu. It also returns the rotation matrices, that
transforms the $$Q$$ point from the original lattice to the selected twin
rlu coordinate system.
  
### Examples
  
This example Calculates the $$[1,0,0]$$ and $$[1,1,0]$$ Bragg reflections
equivalent positions in the twins.
 
```matlab
Q1 = [1 0 0; 1 1 0];
Q2 = cryst.twinq(Q1');
```
  
### Input Arguments
  
`Q0`
: $$Q$$ values in the original crystal in rlu sotred in a matrix with
dimensions of $$[3\times n_Q]$$, optional.
  
### Output Arguments
  
`Qtwin`
: $$Q$$ values in the twin oordinate system in a cell element for
          each twin.
 
`rotQ`
: Rotation matrices with dimensions of $$[3\times 3\times n_{twin}]$$.
  
### See Also
  
[spinw](spinw) \| [spinw.addtwin](spinw_addtwin)
 
*[rlu]: Reciprocal Lattice Unit
 

{% include links.html %}
