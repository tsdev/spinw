---
{title: sw_label, link: sw_label, summary: returns axis labels for spectrum plot,
  keywords: sample, sidebar: sw_sidebar, permalink: sw_label.html, folder: swfiles,
  mathjax: 'true'}

---

### Syntax

`[xlabel, xaxis] = sw_label(hkl,hkla) `

### Description



### Input Arguments

`hkl`
: Momentum transfer values in r.l.u., dimensions are [3 nQ].

`hklA`
: Momentum transfer values in Å$$^{-1}$$, dimensions are [3 nQ].

`lUnit`
: Length unit, given in a string.

### Output Arguments

It returns the label and axis vector for the x-axis for momentum transfer
scans linear in reciprocal space.

### See Also

[sw_plotspec](sw_plotspec.html)

{% include links.html %}
