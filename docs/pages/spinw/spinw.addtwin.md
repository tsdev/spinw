---
{title: spinw.addtwin method, link: spinw.addtwin, summary: adds crystallographic
    twins, keywords: sample, sidebar: sw_sidebar, permalink: spinw_addtwin.html, folder: spinw,
  mathjax: 'true'}

---
  
### Syntax
  
`addtwin(obj,Name,Value)`
  
### Description
  
`addtwin(obj,Name,Value)` adds crystallographic twins defined by a
rotation matrix and its volume fraction. Using crystallographic twins,
SpinW can simulate imperfect samples and if the relative orientation of
the crystallographic twins are knows, SpinW simulations can be directly
compared to the expeiments on the inperfect sample.
  
### Examples
  
This example shows how to add two extra crystallographic twins to the
crystal.  Together with the original orientation there will be three
twins with equal volumes.
 
```matlab
cryst.addtwin('axis',[0 0 1],'phid',[60 120],'vol',[1 1])
```
  
### Input Arguments
  
`obj`
: [spinw](spinw.html) object.
  
### Name-Value Pair Arguments
  
`'axis'`
: Defines the axis of rotation to generate twins in the xyz
  coordinate system, dimensions are $$[1\times 3]$$.
  
`'phi'`
: Defines the angle of rotation to generate twins in radian
  units. Several twins can be defined parallel if `phi` is a
  row vector. Dimensions are $$[1\times n_{twin}]$$.
  
`'phid'`
: Alternative to `phi` but the unit is °.
  
`'rotC'`
: Rotation matrices, that define crystallographic twins. This is an
  alternative to the `axis`-`phi` parameter pair. Matrix dimensions are 
  $$[3\times 3\times n_{twin}]$$.
  
`'vol'`
: Volume fractions of the twins stored in a row vector with $$n_{twin}$$
  elements. Default value is `ones(1,nTwin)`.
  
`'overwrite'`
: If `true`, the last twin will be overwritten, instead of adding a
  new one. Default is `false`.
  
### Output Arguments
  
The function adds extra entries to the [spinw.twin](spinw_twin.html) property.
  
### See Also
  
[spinw](spinw.html) \| [spinw.twinq](spinw_twinq.html)
 

{% include links.html %}
