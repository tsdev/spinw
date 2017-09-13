---
{title: sw_label( ), summary: returns axis labels for spectrum plot, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_label.html, folder: swfiles, mathjax: 'true'}

---
returns axis labels for spectrum plot
 
[xLabel, xAxis] = SW_LABEL(hkl,hklA) 
 
Input:
 
hkl       Momentum transfer values in r.l.u., dimensions are [3 nQ].
hklA      Momentum transfer values in Angstrom^-1, dimensions are [3 nQ].
lUnit     Length unit, given in a string.
 
Output:
 
It returns the label and axis vector for the x-axis for momentum transfer
scans linear in reciprocal space.
 
See also SW_PLOTSPEC.
 
