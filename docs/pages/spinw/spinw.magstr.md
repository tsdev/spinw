---
{title: spinw.magstr method, link: spinw.magstr, summary: returns single-k magnetic
    structure representation, keywords: sample, sidebar: sw_sidebar, permalink: spinw_magstr,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`magout = magstr(obj,Name,Value)`
  
### Description
  
`magout = magstr(obj,Name,Value)` converts the internally stored magnetic
structure (general Fourier representation) into a rotating frame
representation with a single propagation vector, real magnetisation
vectors and the normal axis of the rotation. The conversion is not always
possible, in that case the best possible approximation is used, that
might lead sometimes to unexpected magnetic structures. In this case a
warning is triggered.
  
### Example
 
The example shows the equivalent represenation of a simple spin helix in
the $$ab$$-plane using Fourier components of the magnetization and using
the rotating frame. The complex magnetization in the Fourier
representation is converted into a real spin vector and a normal vector
that defines the axis of rotation.
 
```matlab
model = spinw
model.addatom('r',[0 0 0],'S',1)
model.genmagstr('mode','fourier','S',[1i 1 0]','k',[1/3 0 0])
model.mag_str.F
```
*Output*
```
   0.0000 - 1.0000i
   1.0000 + 0.0000i
   0.0000 + 0.0000i
```
 
```matlab
model.magstr
```
*Output*
```
  struct with fields:
        S: [3×1 double]
        k: [0.3333 0 0]
        n: [0 0 1]
    N_ext: [1 1 1]
    exact: 1
```
 
```matlab
model.magstr.S
```
*Output*
```
     0
     1
     0
```
 
 
### Name-Value Pair Arguments
  
`'exact'`
: If `true`, a warning appears in case the conversion is not exact.
  Default is `true`.
  
`'nExt'`
: Size of the magnetic supercell, default value is the value stored in
  the [spinw](spinw) object (on which the Fourier expansion is defined).
  
`'origin'`
: Origin in lattice units, the magnetic structure will be
  calculated relative to this point. Default value is `[0 0 0]`.
  Shifting the origin introduces an overall phase factor.
  
### See Also
  
[spinw.genmagstr](spinw_genmagstr)
 

{% include links.html %}
