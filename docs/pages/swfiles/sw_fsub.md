---
{title: sw_fsub( ), summary: simple graph vertex coloring, keywords: sample, sidebar: sw_sidebar,
  permalink: sw_fsub.html, folder: swfiles, mathjax: 'true'}

---
simple graph vertex coloring
 
cGraph = SW_FSUB(conn, nExt)
 
It creates a simple graph vertex coloring, determines non connected
sublattices for Monte-Carlo calculation.
 
Input:
 
conn          Contains edge indices which are connected
              conn(1,idx)-->conn(2,idx), dimensions are [2 nConn].
nExt          Size of the magnetic unit cell in units of cells.
 
Output:
 
cGraph        Vector, that assigns every magnetic moment to a sublattice.
 
See also SPINW.ANNEAL.
 

