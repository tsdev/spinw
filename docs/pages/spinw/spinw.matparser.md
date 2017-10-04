---
{title: spinw.matparser method, link: spinw.matparser, summary: parses parameter vector
    into matrices, keywords: sample, sidebar: sw_sidebar, permalink: spinw_matparser,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`matparser(obj,Name,Value)`
  
### Description
  
`matparser(obj,Name,Value)` modifies the `obj.matrix.mat` matrix,
assigning new values from a given parmeter vector.  
  
### Example
 
To assign a Dzyaloshinskii-Moriya vector to the `'DM'` matrix, the
following input would be sufficient:
 
```matlab
cryst = spinw
cryst.addmatrix('label','DM','value',1)
P = [0.2 0.35 3.14]
M = {'DM' 'DM' 'DM'}
S = cat(3,[0 0 0;0 0 1;0 -1 0],[0 0 -1;0 0 0;1 0 0],[0 1 0;-1 0 0;0 0 0])
cryst.matparser('param',P,'mat',M,'selector',S)
cryst.matrix.mat
```
*Output*
```
    1.0000    3.1400   -0.3500
   -3.1400    1.0000    0.2000
    0.3500   -0.2000    1.0000
```
 
 
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
  
`'param'`
: Input row vector `P` with `nPar` elements that contains the
  new values to be assignd to elements of `obj.matrix.mat`
  matrix.
  
`'mat'`
: Identifies which matrices to be changed according to their
  label or index. To select matrices with given labels use a
  cell of strings with $$n_{par}$$ elements, for example
  `M = {'J1','J2'}`. This will change the diagonal elements of
  matrices $$J_1$$ and $$J_2$$ to a given value that is provided in the
  `param` parameter vector. Alternatively the index of the matrices can
  be given in a vector, such as `[1 2]` (index runs according
  to the order of the previous creation of the matrices using
  [spinw.addmatrix](spinw_addmatrix)).
 
  To assign parameter value only to a selected element of a
  $$[3\times 3]$$ matrix, a bracket notation can be used in any string,
  such as `'D(3,3)'`, in this case only the $$(3,3)$$ element of
  the $$[3\times 3]$$ matrix of `'D'` will be modified, the other elements
  will be unchanged. To modify multiple elements of a matrix
  at once, use the option `selector`.
  
`'selector'`
: Matrix with dimensions of $$[3\times 3\times n_{par}]$$. Each `S(:,:,i)`
  submatrix can contain $$\pm 1$$ and 0. Where `S(:,:,i)` contains
  1, the corresponding matrix elements of
  `spinw.matrix.mat(:,:,M(i))` will be changed to the value
  `P(i)*S(:,:,i)` where `P(i)` is the corresponding parameter
  value. 
  
`'init'`
: Initialize the matrices of `obj.matrix.mat` with zeros for all
  selected labels before assigning parameter values. Default
  is `false`.
  
### Output Arguments
  
The [spinw](spinw) object will contain the modified `obj.matrix.mat` field.
  
### See Also
  
[spinw](spinw) \| [spinw.horace](spinw_horace) \| [spinw.addmatrix](spinw_addmatrix)
 

{% include links.html %}
