---
{title: spinw.addmatrix method, link: spinw.addmatrix, summary: 'adds new [3x3] matrix',
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_addmatrix, folder: spinw,
  mathjax: 'true'}

---
  
### Syntax
  
`addmatrix(obj,Name,Value)`
  
### Description
  
`addmatrix(obj,Name,Value)` adds a new $$[3\times 3]$$ matrix to the
[spinw.matrix](spinw_matrix) field of `obj`. The added matrices can be later assigned
to bonds, single ion anisotropy terms or g-tensors of magnetic atoms. If
the given matrix label already exists in `obj`, instead of adding new
matrix the existing one will be overwritten.
  
### Examples
  
The first example adds a diagonal matrix `eye(3)`, that can describe
Heisenberg interaction if assigned to a bond. The second example adds an
ansisymmetric matrix that can decribe Dzyaloshinskii-Moriya (DM)
interaction if assigned to a bond.
 
```matlab
crystal = spinw
crystal.addmatrix('value',1,'label','J_1')
crystal.matrix.mat
```
*Output*
```
     1     0     0
     0     1     0
     0     0     1
```
 
```matlab
crystal.addmatrix('value',[1 0 0],'label','J_1')
crystal.matrix.mat
```
*Output*
```
     0     0     0
     0     0     1
     0    -1     0
```
 
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
  
`'value'`
: The actual numerical values to be added as a matrix. It can have the
  following shapes:
  * $$[3\times 3]$$ the given values will be stored in [spinw.matrix](spinw_matrix) as
    they are given.
  * $$[1\times 1]$$ the given value will be multiplied with `eye(3)`.
  * `[Mx My Mz]` the given triplet will be used to define an
    antisymmetric matrix `M = [0 M3 -M2;-M3 0 M1;M2 -M1 0]`. 
  
`'label'`
: Label string for plotting default value is `'matI'`, where $$I$$ is the index
  of the matrix.
  
`'color'`
: Color for plotting, either row vector
  that contains color RGB codes (values of 0-255), or a string with the
  name of the color, for possible colors names [swplot.color](swplot_color). Default
  color is a random color.
  
### Output Arguments
  
The `obj` output will contain the additional matrix in the [spinw.matrix](spinw_matrix)
field.
  
### See Also
  
[spinw](spinw) \| [swplot.color](swplot_color)
 
[DM]: Dzyaloshinski-Moriya
[RGB]: Red-Green-Blue
 

{% include links.html %}
