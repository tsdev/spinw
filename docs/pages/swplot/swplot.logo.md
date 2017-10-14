---
{title: swplot.logo, link: swplot.logo, summary: creates the SpinW logo, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_logo, folder: swplot, mathjax: 'true'}

---
  
### Syntax
  
`swplot.logo`
  
`swplot.logo(fName)`
 
### Description
  
`swplot.logo` creates and displays the SpinW logo with credentials. The
logo is using an honest colormap [cm_inferno] and removed the tiles of
the sine wave to symbolize the increase of quality of code (measured as a
number of eliminated for loops) :D The logo is used for SpinW 3.0.
Colormap is expected to change for every major version jump.
  
### Examples
 
This is the logo:
 
```matlab
swplot.logo
```
 
{% include image.html file="generated/swplot_1.png" alt="swplot.logo" %}
 
### Input Arguments
  
`fName`
: File name to save the logo. Optional, if not given the logo
  will be shown in a new figure.
  
### See Also
  
[spinw](spinw)
 

{% include links.html %}
