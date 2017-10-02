---
{title: sw_plotsf, link: sw_plotsf, summary: plots the structure factor in the selected
    Q range in 1D or 2D, keywords: sample, sidebar: sw_sidebar, permalink: sw_plotsf,
  folder: swfiles, mathjax: 'true'}

---

### Syntax

`fhandle = sw_plotsf(sfact,Name,Value)`

### Description



### Name-Value Pair Arguments

`'range'`
: Data range in inverse Å, dimensions are [1 2] or [2 2]
  for 1D and 2D plots respectively, default value is the full data
  range.

`'log'`
: Plot 10 based logarithmic intensity, default value is false.

`'colorbar'`
: To show colorbar, default value is true.

`'plotStyle'`
: Options to the plot command, default value is '-'.

### See Also

[spinw](spinw) \| [spinw.structfact](spinw_structfact) \| [sw_intsf](sw_intsf)

{% include links.html %}
