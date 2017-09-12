---
{title: sw_plotsf( ), summary: plots the structure factor in the selected Q range
    in 1D or 2D, keywords: sample, sidebar: sw_sidebar, permalink: swfiles_sw_plotsf.html,
  folder: swfiles, mathjax: 'true'}

---
plots the structure factor in the selected Q range in 1D or 2D
 
fHandle = SW_PLOTSF(sFact, 'Option1', Value1, ...)
 
Options:
 
range     Data range in inverse Angstrom, dimensions are [1 2] or [2 2]
          for 1D and 2D plots respectively, default is the full data
          range.
log       Plot 10 based logarithmic intensity, default is false.
colorbar  To show colorbar, default is true.
plotStyle Options to the plot command, default is '-'.
 
See also SPINW, SPINW.STRUCTFACT, SW_INTSF.
 
