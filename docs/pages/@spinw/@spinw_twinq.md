---
{title: twinq( ), keywords: sample, summary: calculates equivalent Q point in twins,
  sidebar: sw_sidebar, permalink: '@spinw_twinq.html', folder: '@spinw', mathjax: 'true'}

---
  calculates equivalent Q point in twins
 
  [Qtwin, rotQ] = TWINQ(obj, {Q0})
 
  Qtwin are the Q values in the twin coordinate systems, in r.l.u.
  rotQ are the rotation matrices, that transforms Q points (in r.l.u.
  units) from the original crystal to the selected twin r.l.u. coordinate
  system.
 
  Input:
 
  Q0        Q values in the original crystal in r.l.u. Dimensions are
            [3 nQ]. Optional.
 
  Output:
 
  Qtwin     Q values in the twin oordinate system in a cell element for
            each twin.
  rotQ      Rotation matrices, dimensions are [3 3 nTwin].
 
  Example:
 
  ...
  Q1 = [1 0 0; 1 1 0];
  Q2 = cryst.twinq(Q1');
 
  This example Calculates the [1 0 0] and [1 1 0] Bragg reflections
  equivalent position in the twins.
 
  See also SPINW, SPINW.ADDTWIN.
 
