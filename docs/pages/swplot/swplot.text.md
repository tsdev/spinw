---
{title: swplot.text, link: swplot.text, summary: draws a text at a point in 3D, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_text, folder: swplot, mathjax: 'true'}

---

### Syntax

`htext = swplot.text(r, string, {fontsize})`

### Description

hText = SWPLOT.TEXT(handle,...)
 
Handle of an axes object that selects an axis to plot.
 

### Input Arguments

`handle`
: Handle of an axis object.

`r`
: Coordinate of the center of the text for a single text or
  matrix with dimensions [3 nText] for multiple text.

`string`
: String contains the text or cell of strings to plot multiple
  text.

`fontSize`
: Font size in pt, default value is stored in
  swpref.getpref('fontsize')

### See Also

[text]

{% include links.html %}
