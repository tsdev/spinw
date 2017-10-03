---
{title: spinw.gencoupling method, link: spinw.gencoupling, summary: generates bond
    list, keywords: sample, sidebar: sw_sidebar, permalink: spinw_gencoupling, folder: spinw,
  mathjax: 'true'}

---
 
### Syntax
 
`gencoupling(obj,Name,Value)`
 
### Description
 
`gencoupling(obj,Name,Value)` generates all bonds up to a certain length
between magnetic atoms. It also groups bonds based either on crystal
symmetry (is space group is not $$P0$$) or bond length (with `tolDist`
tolerance) is space group is not defined. Sorting bonds based on length
can be forced by setting the `forceNoSym` parameter to true. To check
whether a space group is defined call the [spinw.symmetry](spinw_symmetry) function.
 
{% include warning.html content=" This function has to be used after the crystal structure is defined.
  The [spinw.addcoupling](spinw_addcoupling) function will only work afterwards. " %}
 
### Examples
 
A triangular lattice is generated using `spinw.gencoupling` and
the [spinw.table](spinw_table) function lists the 1st, 2nd and 3rd neighbor bonds:
 
```matlab
cryst = spinw
cryst.genlattice('lat_const',[3 3 5],'angled',[90 90 120])
cryst.addatom('r',[0 0 0],'S',1)
cryst.gencoupling
cryst.table('bond',1:3)
```
 
### Input Arguments
 
`obj`
: [spinw](spinw) object.
 
### Name-Value Pair Arguments
 
`'forceNoSym'`
: If `true`, equivalent bonds are always generated based on
  bond length with `tolDist` length tolerance and effectively reducing
  the bond symmetry to `P0`. If `false` symmetry operators will be used
  if they are given ([spinw.symmetry](spinw_symmetry) returns `true`).
 
`'maxDistance'`
: Maximum bond length that will be stored in the
  [spinw.coupling](spinw_coupling) property in units of Å. Default value is 8.
 
`'maxSym'`
: Maximum bond length until the symmetry equivalent bonds are
  generated. It is usefull if long bonds have to be generated for the
  dipolar interaction, but the symmetry analysis of them is not
  necessary. Default value is equal to `maxDistance`.
 
`'tolDist'`
: Tolerance of distance, within two bonds are considered
  equivalent, default value is $$10^{-3}$$Å. Only used, when no
  space group is defined.
 
`'dMin'`
: Minimum bond length, below which an error is triggered.
  Default value is 0.5 Å.
 
### Output Arguments
 
The [spinw.coupling](spinw_coupling) field of `obj` will store the new bond list, while
overwriting previous bond list. This will also remove any previous
assignment of exchange matrices to bonds.
 
### See Also
 
[spinw](spinw) \| [spinw.symmetry](spinw_symmetry) \| [spinw.nosym](spinw_nosym)
 

{% include links.html %}
