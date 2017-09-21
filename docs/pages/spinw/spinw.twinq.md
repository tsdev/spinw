---
{title: spinw.twinq method, link: spinw.twinq, summary: calculates equivalent Q point
    in twins, keywords: sample, sidebar: sw_sidebar, permalink: spinw_twinq.html,
  folder: spinw, mathjax: 'true'}

---

### Syntax

`[qtwin, rotq] = twinq(obj, {q0})`

### Description

Qtwin are the Q values in the twin coordinate systems, in r.l.u.
rotQ are the rotation matrices, that transforms Q points (in r.l.u.
units) from the original crystal to the selected twin r.l.u. coordinate
system.
 

### Examples

...
Q1 = [1 0 0; 1 1 0];
Q2 = cryst.twinq(Q1');
This example Calculates the [1 0 0] and [1 1 0] Bragg reflections
equivalent position in the twins.

### Input Arguments

`Q0`
: Q values in the original crystal in r.l.u. Dimensions are
  [3 nQ]. Optional.

### Output Arguments

Qtwin     Q values in the twin oordinate system in a cell element for
          each twin.
rotQ      Rotation matrices, dimensions are [3 3 nTwin].

### See Also

[spinw](spinw.html) \| [spinw.addtwin](spinw_addtwin.html)

{% include links.html %}
