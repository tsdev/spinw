---
{title: spinw.magstr method, link: spinw.magstr, summary: generates magnetic structure
    for the rotating frame, keywords: sample, sidebar: sw_sidebar, permalink: spinw_magstr.html,
  folder: spinw, mathjax: 'true'}

---

### Syntax

`magout = magstr(obj,Name,Value)`

### Description

The function converts the internally stored magnetic structure (general
Fourier representation) into a rotating frame representation with a
single propagation vector, real magnetisation vectors and the normal axis
of the rotation. The conversion is not always possible, in that case the
best possible approximation is used, that might lead sometimes to
unexpected magnetic structures.
 

### Name-Value Pair Arguments

`exact`
: If true, a warning appears in case the conversion is not exact.
  Default is true.

`nExt`
: Size of the magnetic supercell, default is the value stored in
  the SpinW object (on which the Fourier expansion is defined).

`origin`
: Origin in lattice units, the magnetic structure will be
  calculated relative to this point. Default value is [0 0 0].
  Shifting the origin introduces an overall phase factor.

### See Also

[spinw.genmagstr](spinw_genmagstr.html)

