---
{title: sw_qgrid( ), link: sw_qgrid, summary: creates a Q grid, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_qgrid.html, folder: swfiles, mathjax: 'true'}

---

### Syntax

` `

### Description

uses n linearly independent vectors ("lattice vectors") and bin values
(coordinates in "lattice units" or "lu") to generate the points. It works
similarly as the d3d constructor in Horace.
 

### Name-Value Pair Arguments

% `u`
:   Row vector with 3 elements, determining the first axis in 3D

% `space,`
:default value is [1 0 0].

% `v`
:   Second axis, default value is [0 1 0].

% `w`
:   Third axis, default value is [0 0 1].

% `uoffset`
:   Column vector with 3 elements, determining the offset of origin

% `in`
:(fourth element is accepted but discarded).

% `ubin`
:   Bin points along the first axis. Can be a vector with 1, 2 or 3

% ``
:s:
 ]        single value along the u-axis:  B1*u
  B2]     range along the u-axis:         [B1:1/nExt:B2]*u
  dB B2]  range along the u-axis:         [B1:dB:B2]*u

% `vbin`
:   Same as ubin but along the v-axis.

% `wbin`
:   Same as ubin but along the w-axis.

% `nExt`
:   Vector with n-elements that can define fractional bin steps,

% `default`
: is [1 1 1].

% `lab`
:   Cell array of projection axis labels with 3 elements (4th

% `element`
: discarded).

% `n`
:termined by the number of given bins (1<=n<=3), so if only 'ubin'

% `is`
:n, n=1; if 'ubin' and 'vbin' are both given n=2, etc.

### Output Arguments

qGrid     A matrix with dimensions of [3,nAx1,nAx2,...], where nAxi is
the index of points along axis i with 1<=i<=n.

### See Also

[sw_qscan](sw_qscan.html)

