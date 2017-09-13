---
{title: sw_qscan( ), summary: creates linear scans between Q points in 3D, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_qscan.html, folder: swfiles, mathjax: 'true'}

---
creates linear scans between Q points in 3D
 
qOut = SW_QSCAN(qLim)
 
Example:
 
qLim = {[0 1 0] [0 0 0]}
If the last element of qLim is a scalar, it defines the number of point
in each linear scan, by default this value is 100.
qLim = {[0 1 0] [0 0 0] 50}
 
