---
{title: spinw.quickham method, link: spinw.quickham, summary: creates magnetic Hamiltonian
    with a single command, keywords: sample, sidebar: sw_sidebar, permalink: spinw_quickham.html,
  folder: spinw, mathjax: 'true'}

---

### Syntax

`quickham(obj, j)`

### Description

The function generates the bonds from the predefined crystal structure
and assigns exchange values to bonds such as J(1) to first neighbor, J(2)
for second neighbor etc. The command will erase all previous bond,
anisotropy, g-tensor and matrix definitions. Even if J(idx) == 0, the
corresponding bond and matrix will be created.
 

### Input Arguments

`obj`
:pinw] object.

`J`
:    Vector containing the Heisenberg exchange values. J(1) for
     first neighbor bonds, etc.

